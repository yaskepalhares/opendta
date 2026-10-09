"""svy: regress, logit, logistic, probit, poisson e svy: tabulate
([SVY] svy estimation, [SVY] svy: tabulate oneway/twoway).

As estimativas pontuais usam os pesos amostrais; a variância vem da
linearização (sanduíche com o "meio" calculado pelo desenho) e os testes
usam o F de Wald ajustado: F = (d - k + 1) W / (d k), gl (k, d - k + 1).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy import special, stats as st

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import match_options
from . import fvars
from .regress import Sample
from .results import Estimates, coef_table, post, rows_from_columns
from .svy import DesignData, _num, _subpop_mask, require_design

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

_TITLES = {"regress": "Linear regression", "logit": "Logistic regression", "logistic": "Logistic regression",
           "probit": "Probit regression", "poisson": "Poisson regression"}


def g(x, fmt):
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def _adjusted_F(b, V, idx, d):
    k = len(idx)
    if not k:
        return SYS, 0, 0
    try:
        W = float(b[idx] @ np.linalg.solve(V[np.ix_(idx, idx)], b[idx]))
    except np.linalg.LinAlgError:
        return SYS, k, d - k + 1
    F = (d - k + 1) / (d * k) * W
    return F, k, d - k + 1


def svy_model(s: "Session", cmd: str, args: str, svyopts: str) -> None:
    d = require_design(s)
    ds = s.data
    so = match_options(svyopts, {"subpop": 3, "vce": 3}) if svyopts.strip() else {}
    smp = Sample(s, args, opts_spec={"or": 2, "irr": 3, "nolog": 4, "coef": 4, "noheader": 6},
                 weights=("pweight",))
    out = s.output
    dd = DesignData(s, d, np.ones(ds.nobs, bool))
    subpop = _subpop_mask(s, str(so["subpop"])) if so.get("subpop") else None
    inside_all = smp.mask & dd.mask & (subpop if subpop is not None else True)
    full = dd.mask
    inside = inside_all[full]
    # colunas avaliadas em todas as observações do desenho (zeros fora da amostra)
    est_mask = smp.mask & dd.mask
    if not est_mask.any():
        raise StataError(2000, "no observations")
    w = dd.w * inside
    out.write(f"\n(running {cmd} on estimation sample)\n", "text")
    cols = fvars.build_design(ds, smp.terms, full, constant=smp.constant)
    notes = fvars.drop_collinear([c for c in cols], smp.constant, None)
    for nm in notes:
        out.write(f"note: {nm} omitted because of collinearity\n", "text")
    est_cols = [c for c in cols if c.values is not None]
    n = int(full.sum())
    Xc = [np.where(inside, c.values, 0.0) for c in est_cols]
    X = np.column_stack(Xc + ([np.ones(n)] if smp.constant else [])) if (Xc or smp.constant) else np.zeros((n, 0))
    yv = fvars.comp_values(ds, smp.depcomp)[full]
    y = np.where(inside & (yv < SYS), yv, 0.0)
    k = X.shape[1]
    if cmd == "regress":
        XtWX = (X * w[:, None]).T @ X
        A = np.linalg.inv(XtWX)
        b = A @ (X.T @ (w * y))
        e = (y - X @ b) * inside
        Z = (X * (w * e)[:, None]) @ A
        ybar = float(w @ y) / w.sum()
        r2 = 1 - float(w @ e ** 2) / float(w @ ((y - ybar) ** 2 * inside))
        extra = {"r2": r2}
    else:
        from .ml import maximize
        from .models import _logit, _poisson, _probit
        if cmd in ("logit", "logistic", "probit"):
            yb = (y != 0).astype(float)
            fun, scores = (_probit if cmd == "probit" else _logit)(yb, X, w)
            p0 = float(w @ yb) / w.sum()
            b0 = np.zeros(k)
            if smp.constant:
                b0[-1] = np.log(p0 / (1 - p0)) if cmd != "probit" else st.norm.ppf(p0)
        else:
            fun, scores = _poisson(y, X, w, np.zeros(n))
            b0 = np.zeros(k)
            if smp.constant:
                b0[-1] = np.log(float(w @ y) / w.sum())
        res = maximize(fun, b0)
        if not smp.o.get("nolog"):
            out.write("\n", "text")
            from .ml import print_log
            print_log(out, res, "log pseudolikelihood")
        b = res.b
        A = np.linalg.inv(-res.H)
        Z = (X * (w * scores(b) * inside)[:, None]) @ A
        extra = {}
    V = dd.variance(Z)
    V = (V + V.T) / 2
    df = dd.df
    slopes = list(range(len(est_cols)))
    F, df1, df2 = _adjusted_F(b, V, slopes, df)
    names = [c.name if not c.omitted else fvars.omitted_name(c.name) for c in cols] + (["_cons"] if smp.constant else [])
    kk = len(names)
    pos = [i for i, c in enumerate(cols) if c.values is not None] + ([kk - 1] if smp.constant else [])
    bf = np.zeros(kk)
    Vf = np.zeros((kk, kk))
    bf[pos] = b
    Vf[np.ix_(pos, pos)] = V
    bmap = {c.name: bf[i] for i, c in enumerate(cols)}
    semap = {c.name: float(np.sqrt(max(Vf[i, i], 0))) for i, c in enumerate(cols)}
    if smp.constant:
        bmap["_cons"], semap["_cons"] = bf[-1], float(np.sqrt(max(Vf[-1, -1], 0)))
    est = Estimates(cmd if cmd != "logistic" else "logit", smp.depname, names, bf, Vf, int(inside.sum()), "t", df)
    est.rows = rows_from_columns(cols, bmap, semap, ds, constant=smp.constant)
    est.cols, est.constant, est.terms, est.depcomp = cols, smp.constant, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = "Linearized"
    est.ml_cmd = cmd if cmd != "regress" else "regress"
    N = int(est_mask.sum())
    pop = float(dd.w[est_mask[full]].sum())
    est.extra = dict(N=N, pop=pop, df=df, F=F, df1=df1, df2=df2, strata=dd.n_strata, psu=dd.n_psu,
                     sub=(int(inside.sum()), float(w.sum())) if subpop is not None else None, **extra)
    est.display = _display
    eform = ""
    if cmd == "logistic" and not smp.o.get("coef"):
        eform = "Odds Ratio"
    elif smp.o.get("or"):
        eform = "Odds Ratio"
    elif smp.o.get("irr"):
        eform = "IRR"
    est.eform_title = eform
    est.svytitle = _TITLES[cmd]
    from .models import predict_ml
    est.predict = predict_ml if cmd != "regress" else _predict_xb
    scalars = [("N", float(N)), ("N_pop", pop), ("N_strata", float(dd.n_strata)), ("N_psu", float(dd.n_psu)),
               ("df_r", float(df)), ("df_m", float(df1)), ("F", F)]
    if subpop is not None:
        scalars += [("N_sub", float(inside.sum())), ("N_subpop", float(w.sum()))]
    if "r2" in extra:
        scalars.append(("r2", extra["r2"]))
    macros = [("prefix", "svy"), ("cmdname", cmd), ("cmd", cmd if cmd != "logistic" else "logistic"),
              ("cmdline", f"svy: {cmd} {args.strip()}"), ("vce", "linearized"), ("vcetype", "Linearized"),
              ("depvar", smp.depname), ("properties", "b V"), ("title", f"Survey: {_TITLES[cmd]}")]
    post(s, est, scalars, macros, est_mask)
    _display(s, est)


def _predict_xb(s, est, stat, mask):
    from .postest import xb_all
    xb, ok = xb_all(s, est)
    if stat in ("", "xb"):
        return np.where(ok, xb, SYS), ("(option xb assumed; fitted values)" if not stat else "")
    raise StataError(198, f"option {stat} not allowed")


def _display(s: "Session", est: Estimates, *, header: bool = True, table: bool = True, level=None, **_kw):
    out = s.output
    x = est.extra
    if header:
        out.write("\n", "text")
        out.write(f"Survey: {est.svytitle}\n\n", "text")
        right = [("Number of obs", f"{x['N']:,}"), ("Population size", _num(x["pop"]))]
        if x.get("sub"):
            right += [("Subpop. no. obs", f"{x['sub'][0]:,}"), ("Subpop. size", _num(x["sub"][1]))]
        right += [("Design df", f"{x['df']:,}"), (f"F({x['df1']:>4},{x['df2']:>7})", g(x["F"], "%10.2f")),
                  ("Prob > F", g(st.f.sf(x["F"], x["df1"], x["df2"]) if x["F"] < SYS else SYS, "%10.4f"))]
        if "r2" in x:
            right.append(("R-squared", g(x["r2"], "%10.4f")))
        left = [f"Number of strata   = {x['strata']:>9,}", f"Number of PSUs     = {x['psu']:>9,}"]
        for i, (rl, rv) in enumerate(right):
            lft = left[i] if i < len(left) else ""
            out.write(f"{lft:<46}{rl:<18}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
    if table:
        title = getattr(est, "eform_title", "") or "Coef."
        coef_table(s, est.depvar, est.rows, stat="t", df=est.df_r, level=level or est.level,
                   vcetype="Linearized", coef_title=title, eform=bool(getattr(est, "eform_title", "")))


# ---------------------------------------------------------------------------
# svy: tabulate (uma e duas vias)
# ---------------------------------------------------------------------------

def svy_tabulate(s: "Session", args: str, svyopts: str) -> None:
    from ..commands._util import touse
    from ..core.varlist import expand
    from ..lang.syntax import parse_standard
    d = require_design(s)
    ds = s.data
    p = parse_standard(args)
    o = match_options(p.options, {"count": 3, "cell": 2, "row": 1, "column": 3, "se": 2, "ci": 2,
                                  "percent": 3, "format": 3, "missing": 4, "pearson": 4, "obs": 3,
                                  "proportion": 4}) if p.options.strip() else {}
    so = match_options(svyopts, {"subpop": 3}) if svyopts.strip() else {}
    names = expand(ds, p.varlist, allow_empty=False)
    if len(names) > 2:
        raise StataError(103, "too many variables specified")
    dd = DesignData(s, d, np.ones(ds.nobs, bool))
    full = dd.mask
    mask = touse(s, p)
    for n in names:
        mask &= ds.get(n).data < SYS
    if so.get("subpop"):
        mask &= _subpop_mask(s, str(so["subpop"]))
    inside = mask[full]
    w = dd.w * inside
    out = s.output
    W = float(w.sum())
    pct = bool(o.get("percent"))
    mult = 100.0 if pct else 1.0
    vars_ = [ds.get(n) for n in names]
    data = [v.data[full] for v in vars_]
    levels = [sorted(set(x[inside].tolist())) for x in data]

    def lab(v, lv):
        t = ds.value_labels.get(v.value_label, {}).get(int(lv)) if v.value_label else None
        return t if t else (f"{int(lv)}" if lv == int(lv) else f"{lv:g}")
    N = int(inside.sum())
    out.write("\n", "text")
    right = [("Number of obs", f"{N:,}"), ("Population size", _num(W)), ("Design df", f"{dd.df:,}")]
    left = [f"Number of strata   = {dd.n_strata:>9,}", f"Number of PSUs     = {dd.n_psu:>9,}"]
    for i, (rl, rv) in enumerate(right):
        lft = left[i] if i < len(left) else ""
        out.write(f"{lft:<46}{rl:<18}= ", "text")
        out.write(f"{rv:>10}\n", "result")
    out.write("\n", "text")
    stat_name = "percentage" if pct else "proportion"
    if len(names) == 1:
        v0 = vars_[0]
        cells = [(lv, float(w @ (data[0] == lv)) / W) for lv in levels[0]]
        Z = np.column_stack([w * ((data[0] == lv) - p_) / W for lv, p_ in cells])
        V = dd.variance(Z)
        se = np.sqrt(np.maximum(np.diag(V), 0))
        wl = max(9, max(len(lab(v0, lv)) for lv, _ in cells), len(names[0]))
        colw = max(len(stat_name), 10)
        extra_se = bool(o.get("se"))
        line = "-" * (wl + 1) + "-" + "-" * (colw + 2) + ("-" * (colw + 2) if extra_se else "")
        out.write(line + "\n", "text")
        out.write(f"{names[0]:<{wl}} | {stat_name:>{colw}}" + (f" {'se':>{colw}}" if extra_se else "") + "\n", "text")
        out.write("-" * wl + "-+-" + "-" * colw + ("-" * (colw + 1) if extra_se else "") + "\n", "text")
        for i, (lv, p_) in enumerate(cells):
            out.write(f"{lab(v0, lv):>{wl}} | ", "text")
            out.write(f"{_sig(p_ * mult):>{colw}}"
                      + (f" {_sig(se[i] * mult):>{colw}}" if extra_se else "")
                      + "\n", "result")
        out.write(f"{'':>{wl}} |\n", "text")
        out.write(f"{'Total':>{wl}} | ", "text")
        out.write(f"{_sig(mult):>{colw}}\n", "result")
        out.write(line + "\n", "text")
        out.write(f"  Key:  {stat_name}  =  cell {stat_name}\n", "text")   # VERIFICAR
        return
    # duas vias: proporções nas células e teste de Pearson com correção de Rao-Scott
    v0, v1 = vars_
    R, C = len(levels[0]), len(levels[1])
    ind = []
    for a in levels[0]:
        for b in levels[1]:
            ind.append(((data[0] == a) & (data[1] == b)).astype(float))
    P = np.array([float(w @ i_) / W for i_ in ind]).reshape(R, C)
    Z = np.column_stack([w * (i_ - P.ravel()[j]) / W for j, i_ in enumerate(ind)])
    V = dd.variance(Z)
    wl = max(9, max(len(lab(v0, a)) for a in levels[0]))
    cw = max(7, max(len(lab(v1, b)) for b in levels[1]))
    head = f"{names[0]:<{wl}} | " + "".join(f"{lab(v1, b):>{cw}} " for b in levels[1]) + f"{'Total':>{cw}}"
    out.write("-" * len(head) + "\n", "text")
    out.write(f"{'':<{wl}} | {names[1]:<{(cw + 1) * C}}\n", "text")
    out.write(head + "\n", "text")
    out.write("-" * wl + "-+-" + "-" * (len(head) - wl - 3) + "\n", "text")
    for i, a in enumerate(levels[0]):
        out.write(f"{lab(v0, a):>{wl}} | ", "text")
        out.write("".join(f"{_sig(P[i, j] * mult):>{cw}} " for j in range(C))
                  + f"{_sig(P[i].sum() * mult):>{cw}}\n", "result")
    out.write(f"{'':>{wl}} |\n", "text")
    out.write(f"{'Total':>{wl}} | ", "text")
    out.write("".join(f"{_sig(P[:, j].sum() * mult):>{cw}} " for j in range(C))
              + f"{_sig(mult):>{cw}}\n", "result")
    out.write("-" * len(head) + "\n", "text")
    out.write(f"  Key:  {stat_name}s in cells\n", "text")   # VERIFICAR
    # Pearson: X² = n Σ (p - pr pc)² / (pr pc); Rao-Scott de 2ª ordem
    pr, pc = P.sum(1), P.sum(0)
    E = np.outer(pr, pc)
    X2 = N * float(np.sum((P - E) ** 2 / np.where(E > 0, E, 1)))
    df0 = (R - 1) * (C - 1)
    # contrastes de interação (efeito de desenho generalizado, Rao & Scott 1984)
    Cm = _interaction_contrasts(R, C)
    Vsrs = (np.diag(E.ravel()) - np.outer(E.ravel(), E.ravel())) / N
    try:
        A = np.linalg.solve(Cm.T @ Vsrs @ Cm, Cm.T @ V @ Cm)
        tr, tr2 = float(np.trace(A)), float(np.trace(A @ A))
        dnum = tr ** 2 / tr2
        dden = dnum * dd.df
        Fv = X2 / tr
        pv = float(st.f.sf(Fv, dnum, dden))
    except np.linalg.LinAlgError:
        Fv, dnum, dden, pv = SYS, SYS, SYS, SYS
    out.write("\n  Pearson:\n", "text")
    out.write(f"    Uncorrected   chi2({df0})         = {X2:>9.4f}\n", "result")
    out.write(f"    Design-based  F({dnum:.2f}, {dden:.2f}) = {Fv:>9.4f}     P = {pv:.4f}\n", "result")
    s.r = {"N": float(N), "r": float(R), "c": float(C), "cun_Pear": X2, "F_Pear": Fv, "df1_Pear": dnum,
           "df2_Pear": dden, "p_Pear": pv}


def _sig(x: float, digits: int = 4) -> str:
    """%6.4g do Stata: 4 algarismos significativos, sem o zero à esquerda."""
    if not np.isfinite(x):
        return "."
    t = f"{x:.{digits}g}"
    if "e" in t:
        return format_value(x, "%9.0g", pad=False).strip()
    return t[1:] if t.startswith("0.") else ("-" + t[2:] if t.startswith("-0.") else t)


def _interaction_contrasts(R: int, C: int) -> np.ndarray:
    """Matriz RC x (R-1)(C-1) de contrastes de interação (efeitos centrados)."""
    def contr(k):
        m = np.zeros((k, k - 1))
        for j in range(k - 1):
            m[j, j] = 1
            m[k - 1, j] = -1
        return m
    return np.kron(contr(R), contr(C))

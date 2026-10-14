"""svy: regress, logit, logistic, probit, poisson and svy: tabulate
([SVY] svy estimation, [SVY] svy: tabulate oneway/twoway).

Point estimates use the sampling weights; the variance comes from
linearization (a sandwich with the "meat" computed from the design) and the tests
use the adjusted Wald F: F = (d - k + 1) W / (d k), df (k, d - k + 1).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy import special, stats as st

from ..core import missing as M
from ..core.dataset import is_string_type
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
    from .svy import parse_svy_opts
    so = parse_svy_opts(svyopts)
    smp = Sample(s, args, opts_spec={"or": 2, "irr": 3, "nolog": 4, "coef": 4, "noheader": 6},
                 weights=("pweight",))
    out = s.output
    dd = DesignData(s, d, np.ones(ds.nobs, bool), vce=so.get("vce"), mse=so.get("mse"),
                    hadamard=so.get("hadamard"))
    subpop = _subpop_mask(s, str(so["subpop"])) if so.get("subpop") else None
    inside_all = smp.mask & dd.mask & (subpop if subpop is not None else True)
    full = dd.mask
    inside = inside_all[full]
    # columns evaluated on all design observations (zeros outside the sample)
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

    def fit(wv):
        """Estimates and linearized scores (n x k) with weights wv."""
        if cmd == "regress":
            XtWX = (X * wv[:, None]).T @ X
            A = np.linalg.inv(XtWX)
            b_ = A @ (X.T @ (wv * y))
            e = (y - X @ b_) * inside
            Z_ = (X * (wv * e)[:, None]) @ A
            ybar = float(wv @ y) / wv.sum()
            r2_ = 1 - float(wv @ e ** 2) / float(wv @ ((y - ybar) ** 2 * inside))
            return b_, Z_, {"r2": r2_}
        from .ml import maximize
        from .models import _logit, _poisson, _probit
        if cmd in ("logit", "logistic", "probit"):
            yb = (y != 0).astype(float)
            fun, scores = (_probit if cmd == "probit" else _logit)(yb, X, wv)
            p0 = float(wv @ yb) / wv.sum()
            b0 = np.zeros(k)
            if smp.constant:
                b0[-1] = np.log(p0 / (1 - p0)) if cmd != "probit" else st.norm.ppf(p0)
        else:
            fun, scores = _poisson(y, X, wv, np.zeros(n))
            b0 = np.zeros(k)
            if smp.constant:
                b0[-1] = np.log(float(wv @ y) / wv.sum())
        # svy: does not show the iteration log (compat 0601)
        res = maximize(fun, b0)
        A = np.linalg.inv(-res.H)
        Z_ = (X * (wv * scores(res.b) * inside)[:, None]) @ A
        return res.b, Z_, {}

    b, Z, extra = fit(w)
    n_rep = None
    if dd.replicated:
        V, n_rep = dd.rep_variance(lambda ww: fit(ww * inside)[0], b)
    else:
        V = dd.variance(Z)
    V = (V + V.T) / 2
    from .svy import _srs_variances
    srs = _srs_variances(Z, w, dd, int((smp.mask & dd.mask).sum()))
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
    Swor, Swr = np.zeros((kk, kk)), np.zeros((kk, kk))
    Swor[np.ix_(pos, pos)] = srs[0]
    Swr[np.ix_(pos, pos)] = srs[1]
    bmap = {c.name: bf[i] for i, c in enumerate(cols)}
    semap = {c.name: float(np.sqrt(max(Vf[i, i], 0))) for i, c in enumerate(cols)}
    if smp.constant:
        bmap["_cons"], semap["_cons"] = bf[-1], float(np.sqrt(max(Vf[-1, -1], 0)))
    est = Estimates(cmd if cmd != "logistic" else "logit", smp.depname, names, bf, Vf, int(inside.sum()), "t", df)
    est.rows = rows_from_columns(cols, bmap, semap, ds, constant=smp.constant)
    est.cols, est.constant, est.terms, est.depcomp = cols, smp.constant, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = dd.vcetype
    est.ml_cmd = cmd if cmd != "regress" else "regress"
    N = int(est_mask.sum())
    pop = float(dd.w[est_mask[full]].sum())
    est.extra = dict(N=N, pop=pop, df=df, F=F, df1=df1, df2=df2, strata=dd.n_strata, psu=dd.n_psu, n_rep=n_rep,
                     sub=(int(inside.sum()), float(w.sum())) if subpop is not None else None, **extra)
    est.display = _display
    est.svy_effects = dict(V=Vf, srs=(Swor, Swr), labels=[(nm, "") for nm in names], vals=bf, title="Coef.",
                           vcetype=dd.vcetype, over="", level=smp.level, depvar=smp.depname)
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
              ("cmdline", f"svy: {cmd} {args.strip()}"), ("vce", dd.vce), ("vcetype", dd.vcetype),
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
        if x.get("n_rep"):
            right.append(("Replications", f"{x['n_rep']:,}"))   # VERIFY position
        right += [("Design df", f"{x['df']:,}"), (f"F({x['df1']:>4},{x['df2']:>7})", g(x["F"], "%10.2f")),
                  ("Prob > F", g(st.f.sf(x["F"], x["df1"], x["df2"]) if x["F"] < SYS else SYS, "%10.4f"))]
        if "r2" in x:
            right.append(("R-squared", g(x["r2"], "%10.4f")))
        left = [f"Number of strata   = {x['strata']:>9,}", f"Number of PSUs     = {x['psu']:>9,}"]
        for i, (rl, rv) in enumerate(right):
            lft = left[i] if i < len(left) else ""
            out.write(f"{lft:<48}{rl:<18}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
    if table:
        title = getattr(est, "eform_title", "") or "Coef."
        coef_table(s, est.depvar, est.rows, stat="t", df=est.df_r, level=level or est.level,
                   vcetype=est.vcetype, coef_title=title, eform=bool(getattr(est, "eform_title", "")))


# ---------------------------------------------------------------------------
# svy: tabulate (one- and two-way)
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
    out.write("\n(running tabulate on estimation sample)\n\n", "text")
    right = [("Number of obs", f"{N:,}"), ("Population size", _num(W)), ("Design df", f"{dd.df:,}")]
    left = [f"Number of strata   = {dd.n_strata:>9,}", f"Number of PSUs     = {dd.n_psu:>9,}"]
    for i, (rl, rv) in enumerate(right):
        lft = left[i] if i < len(left) else ""
        out.write(f"{lft:<48}{rl:<18}= ", "text")
        out.write(f"{rv:>10}\n", "result")
    out.write("\n", "text")
    stat_name = "percentage" if pct else "proportion"
    if len(names) == 1:
        v0 = vars_[0]
        cells = [(lv, float(w @ (data[0] == lv)) / W) for lv in levels[0]]
        Z = np.column_stack([w * ((data[0] == lv) - p_) / W for lv, p_ in cells])
        V = dd.variance(Z)
        se = np.sqrt(np.maximum(np.diag(V), 0))
        # Stata 14 layout (compat 0601): variable name right-aligned in 9 columns, a column of 10
        # to 11 characters
        wl = max(9, max(len(lab(v0, lv)) for lv, _ in cells), len(names[0]))
        colw = max(len(stat_name), 10)
        extra_se = bool(o.get("se"))
        line = "-" * (wl + 1) + "-" + "-" * (colw + 1) + ("-" * (colw + 1) if extra_se else "")
        out.write(line + "\n", "text")
        out.write(f"{names[0]:>{wl}} | {stat_name:>{colw}}" + (f" {'se':>{colw}}" if extra_se else "") + "\n", "text")
        out.write("-" * wl + "-+-" + "-" * (colw - 1) + ("-" * (colw + 1) if extra_se else "") + "-\n", "text")
        for i, (lv, p_) in enumerate(cells):
            out.write(f"{lab(v0, lv):>{wl}} | ", "text")
            out.write(f"{_sig(p_ * mult):>{colw}}"
                      + (f" {_sig(se[i] * mult):>{colw}}" if extra_se else "")
                      + "\n", "result")
        out.write(f"{'':>{wl}} | \n", "text")
        out.write(f"{'Total':>{wl}} | ", "text")
        out.write(f"{_sig(mult):>{colw}}\n", "result")
        out.write(line + "\n", "text")
        out.write(f"  Key:  {stat_name}  =  cell {stat_name}\n", "text")   # VERIFY
        return
    # two-way: cell proportions and Pearson test with the Rao-Scott correction
    v0, v1 = vars_
    R, C = len(levels[0]), len(levels[1])
    ind = []
    for a in levels[0]:
        for b in levels[1]:
            ind.append(((data[0] == a) & (data[1] == b)).astype(float))
    P = np.array([float(w @ i_) / W for i_ in ind]).reshape(R, C)
    Z = np.column_stack([w * (i_ - P.ravel()[j]) / W for j, i_ in enumerate(ind)])
    V = dd.variance(Z)
    # layout (compat 0601): columns as wide as the longest label or value, separated by two spaces;
    # column name centered over them
    wl = max(9, max(len(lab(v0, a)) for a in levels[0]), len(names[0]))
    texts = [[_sig(P[i, j] * mult) for j in range(C)] + [_sig(P[i].sum() * mult)] for i in range(R)]
    tot = [_sig(P[:, j].sum() * mult) for j in range(C)] + [_sig(mult)]
    heads = [lab(v1, b) for b in levels[1]] + ["Total"]
    cw = max([len(h) for h in heads] + [len(x) for row in texts for x in row] + [len(x) for x in tot])
    body = len(heads) * cw + 2 * (len(heads) - 1) + 1
    out.write("-" * (wl + 2 + body) + "\n", "text")
    out.write(f"{'':>{wl}} |" + names[1].center(body) + "\n", "text")
    out.write(f"{names[0]:>{wl}} | " + "  ".join(f"{h:>{cw}}" for h in heads) + "\n", "text")
    out.write("-" * wl + "-+" + "-" * body + "\n", "text")
    for i, a in enumerate(levels[0]):
        out.write(f"{lab(v0, a):>{wl}} | ", "text")
        out.write("  ".join(f"{x:>{cw}}" for x in texts[i]) + "\n", "result")
    out.write(f"{'':>{wl}} | \n", "text")
    out.write(f"{'Total':>{wl}} | ", "text")
    out.write("  ".join(f"{x:>{cw}}" for x in tot) + "\n", "result")
    out.write("-" * (wl + 2 + body) + "\n", "text")
    out.write(f"  Key:  cell {stat_name}\n", "text")
    # Pearson: X² = n Σ (p - pr pc)² / (pr pc); 2nd-order Rao-Scott
    pr, pc = P.sum(1), P.sum(0)
    E = np.outer(pr, pc)
    X2 = N * float(np.sum((P - E) ** 2 / np.where(E > 0, E, 1)))
    df0 = (R - 1) * (C - 1)
    # interaction contrasts (generalized design effect, Rao & Scott 1984) [SVY] svy: tabulate
    # twoway: Δ = (C'D⁻¹V_srs D⁻¹C)⁻¹ (C'D⁻¹ V D⁻¹C), with D = diag(p̂) and V_srs of the observed
    # proportions (compat 0601: F = 4.9245)
    Cm = _interaction_contrasts(R, C)
    pv_ = P.ravel()
    Vsrs = (np.diag(pv_) - np.outer(pv_, pv_)) / N
    Di = np.diag(np.where(pv_ > 0, 1 / np.where(pv_ > 0, pv_, 1), 0.0))
    try:
        A = np.linalg.solve(Cm.T @ Di @ Vsrs @ Di @ Cm, Cm.T @ Di @ V @ Di @ Cm)
        tr, tr2 = float(np.trace(A)), float(np.trace(A @ A))
        dnum = tr ** 2 / tr2
        dden = dnum * dd.df
        Fv = X2 / tr
        pv = float(st.f.sf(Fv, dnum, dden))
    except np.linalg.LinAlgError:
        Fv, dnum, dden, pv = SYS, SYS, SYS, SYS
    out.write("\n  Pearson:\n", "text")
    def _dfmt(x):
        return f"{x:.0f}" if abs(x - round(x)) < 1e-9 else f"{x:.2f}"   # VERIFY non-integer format
    ftxt = f"F({_dfmt(dnum)}, {_dfmt(dden)})"
    out.write(f"    Uncorrected   {'chi2(' + str(df0) + ')':<16}= {X2:>9.4f}\n", "result")
    out.write(f"    Design-based  {ftxt:<16}= {Fv:>9.4f}     P = {pv:.4f}\n", "result")
    s.r = {"N": float(N), "r": float(R), "c": float(C), "cun_Pear": X2, "F_Pear": Fv, "df1_Pear": dnum,
           "df2_Pear": dden, "p_Pear": pv}


def _sig(x: float, digits: int = 4) -> str:
    """Stata's %6.4g: 4 significant digits, without the leading zero."""
    if not np.isfinite(x):
        return "."
    t = f"{x:.{digits}g}"
    if "e" in t:
        return format_value(x, "%9.0g", pad=False).strip()
    return t[1:] if t.startswith("0.") else ("-" + t[2:] if t.startswith("-0.") else t)


def _interaction_contrasts(R: int, C: int) -> np.ndarray:
    """RC x (R-1)(C-1) matrix of interaction contrasts (centered effects)."""
    def contr(k):
        m = np.zeros((k, k - 1))
        for j in range(k - 1):
            m[j, j] = 1
            m[k - 1, j] = -1
        return m
    return np.kron(contr(R), contr(C))


# ---------------------------------------------------------------------------
# svy: with the other models (ologit, oprobit, mlogit, nbreg, glm, tobit)
# ---------------------------------------------------------------------------

_GENERIC_TITLES = {"ologit": "Ordered logistic regression", "oprobit": "Ordered probit regression",
                   "mlogit": "Multinomial logistic regression", "nbreg": "Negative binomial regression",
                   "glm": "Generalized linear models", "tobit": "Tobit regression"}


def svy_generic(s: "Session", cmd: str, args: str, svyopts: str) -> None:
    """Runs the command with the design weights and replaces the variance with the
    design one: _Prep.vce receives the scores and returns the linearized V (or, with
    a replication vce, the command runs again with each replicate weight)."""
    from ..core.dataset import Variable
    from ..lang.syntax import parse_standard
    from .svy import parse_svy_opts
    d = require_design(s)
    ds = s.data
    so = parse_svy_opts(svyopts)
    p = parse_standard(args)
    if p.weight:
        raise StataError(101, "weights not allowed with the svy prefix")   # VERIFY
    if "vce(" in p.options.replace(" ", ""):
        raise StataError(198, "option vce() not allowed")   # VERIFY
    dd = DesignData(s, d, np.ones(ds.nobs, bool), vce=so.get("vce"), mse=so.get("mse"),
                    hadamard=so.get("hadamard"))
    subpop = _subpop_mask(s, str(so["subpop"])) if so.get("subpop") else None
    full = dd.mask
    inside_all = full & (subpop if subpop is not None else True)
    pos = np.full(ds.nobs, -1)
    pos[np.flatnonzero(full)] = np.arange(dd.n)
    out = s.output
    out.write(f"\n(running {cmd} on estimation sample)\n", "text")
    wname, iname = "__svy_w", "__svy_in"
    for nm in (wname, iname):
        if ds.has(nm):
            ds.drop_vars([nm])
    wfull = np.zeros(ds.nobs)
    wfull[full] = dd.w
    ds.add(Variable(wname, "double", wfull))
    ds.add(Variable(iname, "byte", inside_all.astype(float)))
    cond = f"{iname}" + (f" & ({p.if_})" if p.if_ else "")
    rng = f" in {p.in_}" if p.in_ else ""
    opts = (p.options + " " if p.options else "") + "vce(robust)"
    if cmd != "tobit":
        opts += " nolog"
    text = f"{cmd} {p.varlist} [pweight={wname}] if {cond}{rng}, {opts}"
    ctx = {"dd": dd, "pos": pos}
    s._svy_ctx = ctx
    out.quiet_depth += 1
    try:
        s.interp._execute_expanded(text)
    finally:
        out.quiet_depth -= 1
        s._svy_ctx = None
    est = s._estimates
    b = np.asarray(est.b, dtype=float)
    n_rep = None
    try:
        if dd.replicated:
            def estimator(ww):
                wv = np.zeros(ds.nobs)
                wv[full] = ww
                ds.get(wname).data = wv
                out.quiet_depth += 1
                try:
                    s.interp._execute_expanded(text.replace(", vce(robust)", ",").replace("vce(robust)", ""))
                finally:
                    out.quiet_depth -= 1
                return np.asarray(s._estimates.b, dtype=float)
            saved = (s.e, s._estimates, ds.esample)
            V, n_rep = dd.rep_variance(estimator, b)
            s.e, s._estimates, ds.esample = saved
            est = s._estimates
            est.V = (V + V.T) / 2
            _refresh_rows(est)
    finally:
        ds.drop_vars([wname, iname])
    V = np.asarray(est.V, dtype=float)
    df = dd.df
    names = list(est.names)
    eqs = list(est.eqnames) if est.eqnames else ["_"] * len(names)
    aux_eq = {"lnalpha", "sigma", "/"}
    slopes = [j for j, (nm, eq) in enumerate(zip(names, eqs))
              if nm != "_cons" and V[j, j] > 0 and not eq.startswith("cut") and eq not in aux_eq
              and "b." not in nm and "o." not in nm]
    F, df1, df2 = _adjusted_F(b, V, slopes, df)
    N = int(np.count_nonzero(ds.esample)) if ds.esample is not None else int(inside_all.sum())
    pop = float(dd.w[(ds.esample if ds.esample is not None else inside_all)[full]].sum())
    if subpop is not None:
        # with subpop(), Number of obs and Population size count the whole design sample (complete
        # cases in the model variables); compat 0602
        import re as _re
        ok = full.copy()
        for tok in set(_re.findall(r"[A-Za-z_][A-Za-z0-9_]*", _re.sub(r"\b[A-Za-z]+\.", " ", p.varlist))):
            if ds.has(tok):
                v = ds.get(tok)
                if not is_string_type(v.vtype):
                    ok &= ~(v.data.astype(float) >= M.MISSING_THRESHOLD)
        N = int(ok.sum())
        pop = float(dd.w[ok[full]].sum())
    est.vcetype = dd.vcetype
    est.stat = "t"
    est.df_r = df
    est.extra = dict(getattr(est, "extra", {}) or {})
    est.extra.update(dict(N=N, pop=pop, df=df, F=F, df1=df1, df2=df2, strata=dd.n_strata, psu=dd.n_psu,
                          n_rep=n_rep, sub=(int(inside_all.sum()), float(dd.w[inside_all[full]].sum()))
                          if subpop is not None else None))
    est.svytitle = _GENERIC_TITLES.get(cmd, cmd)
    est.display = _display_generic
    est.extra["no_F"] = cmd == "glm"            # svy: glm does not show the F test (compat 0602)
    est.extra["dispersion"] = getattr(est, "dispersion", "") if cmd == "nbreg" else ""
    # under svy the cut points get t and P>|t| (compat 0602)
    for r in est.rows:
        if r.kind == "aux" and r.label.startswith("/cut") and not r.eq:
            r.kind = "coe"
    # exponentiated form requested in the command's options
    from ..lang.syntax import parse_options
    words = {n.lower() for n, _ in parse_options(p.options)} if p.options.strip() else set()
    title_ef = ""
    if "eform" in words or "ef" in words:
        title_ef = "exp(b)" if cmd == "glm" else getattr(est, "eform_title_opt", "exp(b)")
    elif "or" in words:
        title_ef = "Odds Ratio"
    elif "irr" in words:
        title_ef = "IRR"
    elif "rrr" in words:
        title_ef = "RRR"
    est.coef_title, est.eform_on = (title_ef, True) if title_ef else ("Coef.", False)
    Z = ctx.get("Z")
    if Z is not None:
        from .svy import _srs_variances
        srs = _srs_variances(Z, dd.w * inside_all[full], dd, N)
        labs = [("/" + eq if (eq.startswith("cut") or eq in aux_eq) and nm == "_cons" else nm, "")
                for nm, eq in zip(names, eqs)]
        est.svy_effects = dict(V=V, srs=srs, labels=labs, vals=b, title="Coef.",
                               vcetype=dd.vcetype, over="", level=est.level,
                               depvar=getattr(est, "depvar", ""))
    e = s.e
    e.update({"N_strata": float(dd.n_strata), "N_psu": float(dd.n_psu), "N_pop": pop, "df_r": float(df),
              "F": F, "df_m": float(df1)})
    e["prefix"] = "svy"
    e["vce"] = dd.vce
    e["vcetype"] = dd.vcetype
    e["cmdline"] = f"svy: {cmd} {args.strip()}"
    from ..commands.matrix import Matrix
    if isinstance(e.get("V"), Matrix):
        e["V"] = Matrix(V.copy(), list(e["V"].rownames), list(e["V"].colnames))
        e["V"].roweq = eqs
        e["V"].coleq = eqs
    _display_generic(s, est)


def _refresh_rows(est) -> None:
    """Standard errors of the table rows from est.V (same order as e(b))."""
    se = np.sqrt(np.maximum(np.diag(est.V), 0))
    k = 0
    for r in est.rows:
        if r.kind in ("coef",):
            while k < len(est.b) and np.asarray(est.V)[k, k] == 0 and est.b[k] == 0:
                k += 1
            if k < len(se):
                r.se = float(se[k])
            k += 1
        elif r.kind == "aux" and not r.eq:
            if k < len(se) and r.se < SYS:     # a missing standard error stays missing
                r.se = float(se[k])
            k += 1


def _display_generic(s: "Session", est, *, header: bool = True, table: bool = True, level=None, **_kw):
    out = s.output
    x = est.extra
    if header:
        out.write("\n", "text")
        out.write(f"Survey: {est.svytitle}\n\n", "text")
        right = [("Number of obs", f"{x['N']:,}"), ("Population size", _num(x["pop"]))]
        if x.get("sub"):
            right += [("Subpop. no. obs", f"{x['sub'][0]:,}"), ("Subpop. size", _num(x["sub"][1]))]
        if x.get("n_rep"):
            right.append(("Replications", f"{x['n_rep']:,}"))
        right += [("Design df", f"{x['df']:,}")]
        if not x.get("no_F"):
            right += [(f"F({x['df1']:>4},{x['df2']:>7})", g(x["F"], "%10.2f")),
                      ("Prob > F", g(st.f.sf(x["F"], x["df1"], x["df2"]) if x["F"] < SYS else SYS, "%10.4f"))]
        left = [f"Number of strata   = {x['strata']:>9,}", f"Number of PSUs     = {x['psu']:>9,}"]
        if x.get("dispersion"):
            # "Dispersion         = mean" on the Prob > F line (compat 0602)
            left += [""] * (len(right) - 1 - len(left)) + [f"{'Dispersion':<19}= {x['dispersion']}"]
        for i, (rl, rv) in enumerate(right):
            lft = left[i] if i < len(left) else ""
            out.write(f"{lft:<48}{rl:<18}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="t", df=est.df_r, level=level or est.level,
                   vcetype=est.vcetype, coef_title=getattr(est, "coef_title", "Coef."),
                   eform=getattr(est, "eform_on", False))

"""Hotelling's test, multivariate regression, MANOVA, canonical correlations and mvtest
([MV] hotelling, [MV] mvreg, [MV] manova, [MV] canon, [MV] mvtest).

The statistics are in core/mvtests.py. Layouts follow the manuals' examples; what they do not
show carries a VERIFY mark and is settled by the reference log of compat case 1906."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import multivariate as mv
from ..core import mvtests as T
from ..core.errors import StataError
from ..lang.syntax import match_options, parse_standard
from .mvcommon import g as gnum, lower_matrix, matrix, numeric_names, post_e, sample
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _fp(p: float) -> str:
    return f"{p:.4f}" if np.isfinite(p) else "."


def _summary_table(out, ds, names, X, w=None, wtype="") -> None:
    """The table of `summarize`: observations, mean, standard deviation, minimum and maximum."""
    from ..core import stats as S
    from .summarize import g9
    out.write("\n    Variable |        Obs        Mean    Std. Dev.       Min        Max\n", "text")
    out.write("-------------+" + "-" * 57 + "\n", "text")
    for j, n_ in enumerate(names):
        m = S.moments(X[:, j], w, wtype)
        # the mean is accumulated in extended precision, so that a mean that is exactly -.86921875 shows as Stata does
        xj = X[:, j].astype(np.longdouble)
        wj = np.ones(len(xj), dtype=np.longdouble) if w is None else np.asarray(w, dtype=np.longdouble)
        mean = float((wj * xj).sum() / wj.sum())
        out.write(f"{n_[:12]:>12} | {m.N:>10,.0f}   {g9(mean):>9}   {g9(m.sd):>9}  {g9(m.min):>9}  {g9(m.max):>9}\n", "text")


# ---------------------------------------------------------------------------
# hotelling
# ---------------------------------------------------------------------------

@command("hotelling", "hot")
def cmd_hotelling(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"by": 2, "notable": 3}) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 1)
    by = str(o["by"]).strip() if isinstance(o.get("by"), str) else None
    X, w, wtype, mask = sample(s, p, names, extra=[by] if by else [])
    out = s.output
    k = len(names)
    ds = s.data
    if by is None:
        mean, S, R, sw = mv.weighted_moments(X, w)
        n = sw
        if not o.get("notable"):
            _summary_table(out, ds, names, X, w, wtype)
        r = T.hotelling_one(mean, S, n)
        t2, f_ = gnum(r["T2"], "%10.0g"), gnum(r["F"], "%10.0g")
        out.write(f"\n1-group Hotelling's T-squared = {t2}\n", "text")
        out.write(f"F test statistic: (({int(n)}-{k})/({int(n)}-1)({k})) x {t2} = {f_}\n", "text")
        out.write("\nH0: Vector of means is equal to a vector of zeros\n", "text")
        out.write(f"{'F(' + str(k) + ',' + gnum(r['df2']) + ') =':>24}{r['F']:>10.4f}\n", "text")
        out.write(f"{'Prob > F(' + str(k) + ',' + gnum(r['df2']) + ') =':>24}{r['p']:>10.4f}\n", "text")
        s.r = {"N": float(n), "k": float(k), "df": r["df2"], "T2": r["T2"]}
        return
    g = s.data.get(by).data[mask]
    levels = sorted(set(g.tolist()))
    if len(levels) != 2:
        raise StataError(420, f"{by} has {len(levels)} groups; must be 2")                 # VERIFY
    parts = []
    for lv in levels:
        sel = g == lv
        ww = None if w is None else w[sel]
        mean, S, R, sw = mv.weighted_moments(X[sel], ww)
        parts.append((mean, S, sw, sel))
        if not o.get("notable"):
            from .inspect import cell_text
            lab = cell_text(ds, ds.get(by), int(np.flatnonzero(mask)[np.flatnonzero(sel)[0]]))
            out.write("\n" + "-" * int(s.settings.get("linesize", 80)) + "\n", "text")
            out.write(f"-> {by} = {lab}\n", "text")
            _summary_table(out, ds, names, X[sel], ww, wtype)
    (m1, S1, n1, _), (m2, S2, n2, _) = parts
    r = T.hotelling_two(m1, S1, n1, m2, S2, n2)
    t2, f_ = gnum(r["T2"], "%10.0g"), gnum(r["F"], "%10.0g")
    out.write(f"\n2-group Hotelling's T-squared = {t2}\n", "text")
    out.write(f"F test statistic: (({int(n1 + n2)}-{k}-1)/({int(n1 + n2)}-2)({k})) x {t2} = {f_}\n", "text")
    out.write("\nH0: Vectors of means are equal for the two groups\n", "text")
    d2 = gnum(r["df2"])
    out.write(f"{'F(' + str(k) + ',' + d2 + ') =':>24}{r['F']:>10.4f}\n", "text")
    out.write(f"{'Prob > F(' + str(k) + ',' + d2 + ') =':>24}{r['p']:>10.4f}\n", "text")
    s.r = {"N": float(n1 + n2), "k": float(k), "df": r["df2"], "T2": r["T2"]}


# ---------------------------------------------------------------------------
# mvreg
# ---------------------------------------------------------------------------

@command("mvreg", "mvreg")
def cmd_mvreg(s: "Session", args: str) -> None:
    if not args.strip() or args.strip().startswith(","):
        from ..estimation.postest import replay
        replay(s, "mvreg", args)
        return
    from ..estimation import fvars
    from ..estimation.results import CoefRow, Estimates, post, rows_from_columns
    from ._util import touse
    from .summarize import weights
    p = parse_standard(args)
    o = match_options(p.options, {"noconstant": 3, "corr": 4, "level": 1}) if p.options.strip() else {}
    if p.exp is None:
        raise StataError(198, "invalid syntax")
    ds = s.data
    deps = numeric_names(s, p.varlist, 1)
    constant = not o.get("noconstant")
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    with fvars.caching():
        terms = fvars.expand_fv(ds, p.exp)
        mask = touse(s, p)
        for d_ in deps:
            mask &= ds.get(d_).data < SYS
        mask &= ~fvars.term_missing(ds, terms)
        w, wtype, mask = weights(s, p, mask, ("fweight",))
        if not mask.any():
            raise StataError(2000, "no observations")
        cols = fvars.build_design(ds, terms, mask, constant=constant)
    n = int(mask.sum())
    wv = np.ones(n) if w is None else w[mask].astype(float)
    N = float(wv.sum())
    for nm in fvars.drop_collinear(cols, constant, wv):
        s.output.write(f"note: {nm} omitted because of collinearity\n", "text")
    est_cols = [c for c in cols if c.values is not None]
    X = np.column_stack([c.values for c in est_cols] + ([np.ones(n)] if constant else []))
    Y = np.column_stack([ds.get(d_).data[mask] for d_ in deps]).astype(float)
    k = X.shape[1]
    XtWXi = np.linalg.inv((X * wv[:, None]).T @ X)
    B = XtWXi @ (X * wv[:, None]).T @ Y
    E = Y - X @ B
    df_r = N - k
    if df_r <= 0:
        raise StataError(2001, "insufficient observations")
    SSCP = (E * wv[:, None]).T @ E
    Sigma = SSCP / df_r
    q = len(deps)
    df_m = k - (1 if constant else 0)
    out = s.output
    names_one = [c.name if not c.omitted else fvars.omitted_name(c.name) for c in cols] + (["_cons"] if constant else [])
    kk = len(names_one)
    pos = [i for i, c in enumerate(cols) if c.values is not None] + ([kk - 1] if constant else [])
    eq_stats = []
    b_all, V_all, names_all, eq_all, rows = [], np.zeros((kk * q, kk * q)), [], [], []
    for e_i, d_ in enumerate(deps):
        bfull = np.zeros(kk)
        bfull[pos] = B[:, e_i]
        bmap = {c.name: bfull[i] for i, c in enumerate(cols)}
        semap = {}
        for i, c in enumerate(cols):
            semap[c.name] = math.sqrt(Sigma[e_i, e_i] * XtWXi[pos.index(i), pos.index(i)]) if i in pos else 0.0
        if constant:
            bmap["_cons"] = bfull[-1]
            semap["_cons"] = math.sqrt(Sigma[e_i, e_i] * XtWXi[-1, -1])
        if e_i:
            rows.append(CoefRow("sep"))
        rows.append(CoefRow("header", d_, left=True))
        rows += rows_from_columns(cols, bmap, semap, ds, constant=constant)
        b_all.append(bfull)
        names_all += names_one
        eq_all += [d_] * kk
        y = Y[:, e_i]
        ybar = float(np.average(y, weights=wv)) if constant else 0.0
        tss = float((wv * (y - ybar) ** 2).sum())
        rss = float(SSCP[e_i, e_i])
        r2 = 1 - rss / tss if tss > 0 else float("nan")
        F = ((tss - rss) / df_m) / (rss / df_r) if df_m > 0 else float("nan")
        pF = float(T.st.f.sf(F, df_m, df_r)) if df_m > 0 and np.isfinite(F) else float("nan")
        eq_stats.append((d_, N, k, math.sqrt(Sigma[e_i, e_i]), r2, F, pF))
    for a in range(q):
        for b_ in range(q):
            blk = Sigma[a, b_] * XtWXi
            full = np.zeros((kk, kk))
            full[np.ix_(pos, pos)] = blk
            V_all[a * kk:(a + 1) * kk, b_ * kk:(b_ + 1) * kk] = full
    b_vec = np.concatenate(b_all)
    est = Estimates("mvreg", deps[0], names_all, b_vec, V_all, int(round(N)), "t", df_r)
    est.eqnames = eq_all
    est.rows = rows
    est.cols = cols
    est.constant = constant
    est.terms = terms
    est.level = level
    est.vcetype = ""
    est.wtype = wtype
    est.resid_df = df_r
    est.extra = {"N": N, "df_m": df_m, "df_r": df_r, "Sigma": Sigma, "B": B, "deps": deps,
                 "rank": k}
    est.predict = _mvreg_predict
    est.predict_opts = {}
    scalars = [("k", kk), ("df_r", df_r), ("k_eq", q), ("N", N), ("rank", kk * q)]
    if o.get("corr") and q > 1:
        Rr = Sigma / np.sqrt(np.outer(np.diag(Sigma), np.diag(Sigma)))
        scalars += [("df_chi2", q * (q - 1) / 2.0), ("chi2", float(N * sum(Rr[i, j] ** 2 for i in range(q) for j in range(i))))]
    strs = lambda vals, fmt: "".join(" " + fmt(v) for v in vals)
    g7 = lambda v: gnum(v, "%9.0g")
    macros = [("cmd", "mvreg"), ("cmdline", "mvreg " + args.strip()), ("estat_cmd", "mvreg_estat"),
              ("marginsdefault", " ".join(f"predict(xb equation({d_}))" for d_ in deps)), ("marginsok", "XB default"),
              ("marginsnotok", "stdp Residuals stddp"), ("predict", "reg3_p"), ("depvar", " ".join(deps)),
              ("eqnames", " ".join(deps)), ("F", strs([x[5] for x in eq_stats], g7)),
              ("rmse", strs([x[3] for x in eq_stats], g7)), ("p_F", strs([x[6] for x in eq_stats], lambda v: f"{v:.4f}")),
              ("r2", strs([x[4] for x in eq_stats], lambda v: f"{v:.4f}")), ("properties", "b V"), ("wtype", wtype)]
    post(s, est, scalars, macros, mask, extra_mats=[("Sigma", matrix(Sigma, deps, deps))])
    est.no_b = False
    _display_mvreg(s, est, eq_stats, o, deps, N)


def _display_mvreg(s, est, eq_stats, o, deps, N) -> None:
    from ..estimation.results import coef_table
    out = s.output
    out.write("\n", "text")
    out.write(f"{'Equation':<21}{'Obs':>3}{'Parms':>8}{'RMSE':>12}{'\"R-sq\"':>10}{'F':>11}{'P':>9}\n", "text")
    out.write("-" * 74 + "\n", "text")
    for d_, ob, parms, rmse, r2, F, pF in eq_stats:
        out.write(f"{d_[:21]:<21}{int(ob):>3}{parms:>8}{gnum(rmse, '%9.0g'):>12}{r2:>10.4f}{gnum(F, '%9.0g'):>11}{pF:>9.4f}\n", "text")
    out.write("\n", "text")
    coef_table(s, "", est.rows, stat="t", df=est.df_r, level=est.level)
    if o.get("corr") and len(deps) > 1:
        q = len(deps)
        Sg = est.extra["Sigma"]
        R = Sg / np.sqrt(np.outer(np.diag(Sg), np.diag(Sg)))
        chi2 = float(N * sum(R[i, j] ** 2 for i in range(q) for j in range(i)))
        df = q * (q - 1) / 2.0
        out.write("\nCorrelation matrix of residuals:\n\n", "text")
        lower_matrix(out, deps, R)
        out.write(f"\nBreusch-Pagan test of independence: chi2({int(df)}) = {chi2:9.3f}, Pr = {T.st.chi2.sf(chi2, df):.4f}\n", "text")


def _mvreg_predict(s, est, stat, mask):
    from ..estimation import fvars
    from ..estimation.postest import _columns
    o = getattr(est, "_pred_opts", {}) or {}
    deps = est.extra["deps"]
    eq = str(o["equation"]).strip().lstrip("#") if isinstance(o.get("equation"), str) else deps[0]
    if eq.isdigit():
        eq = deps[int(eq) - 1]
    if eq not in deps:
        raise StataError(111, f"equation {eq} not found")                                  # VERIFY
    i = deps.index(eq)
    kk = len(est.names) // len(deps)
    b = est.b[i * kk:(i + 1) * kk]
    with fvars.caching():
        ok = ~fvars.term_missing(s.data, est.terms)
        xb = np.zeros(s.data.nobs)
        for j, col in _columns(s, est):
            xb += b[j] * col
    xb = np.where(ok, xb, SYS)
    note = ""
    if stat in ("", "xb"):
        if not stat:
            note = "(option xb assumed; fitted values)"
        return xb, note
    if stat == "residuals":
        y = s.data.get(eq).data
        return np.where(ok & (y < SYS), y - xb, SYS), ""
    raise StataError(198, f"option {stat} not allowed")


# ---------------------------------------------------------------------------
# manova
# ---------------------------------------------------------------------------

def _resid_sscp(X: np.ndarray, Y: np.ndarray, w: np.ndarray):
    sw = np.sqrt(w)
    Xw, Yw = X * sw[:, None], Y * sw[:, None]
    beta, *_ = np.linalg.lstsq(Xw, Yw, rcond=None)
    R = Yw - Xw @ beta
    return R.T @ R, int(np.linalg.matrix_rank(Xw))


def _test_mats(lam, stats, q: int, dfh: float, dfe: float) -> tuple:
    """(eigvals, aux, stat) of one test: the s = min(q, df) nonzero eigenvalues, (s, m, n) and the
    4 x 6 table of the statistics (type 1 = exact, 2 = approximate, 3 = upper bound)."""
    sdim = int(min(q, dfh))
    mpar, npar = (abs(q - dfh) - 1.0) / 2.0, (dfe - q - 1.0) / 2.0
    code = {"e": 1.0, "a": 2.0, "u": 3.0}
    stat = None
    if stats:
        stat = matrix(np.array([[r_["stat"], r_["df1"], r_["df2"], r_["F"], r_["p"], code[r_["flag"]]]
                                for r_ in stats]), ["W", "P", "L", "R"], ["statistic", "df1", "df2", "F", "p", "type"])
    aux = matrix(np.array([[sdim], [mpar], [npar]], dtype=float), ["s", "m", "n"], ["aux"])
    return matrix(np.asarray(lam[:sdim], dtype=float).reshape(1, -1)), aux, stat


def _stat_rows(out, label: str, stats: list[dict], df_h: float) -> None:
    pieces = _wrap_label(label)
    for i, r in enumerate(stats):
        src = pieces[i] if i < len(pieces) else ""
        dfc = f"{int(df_h):>10}" if i == 0 else f"{'':>10}"
        out.write(f"{src:>20} |{r['name']}{r['stat']:>9.4f}{dfc}{r['df1']:>9.1f}{r['df2']:>9.1f}"
                  f"{r['F']:>9.2f} {r['p']:.4f} {r['flag']}\n", "text")


@command("manova", "manova")
def cmd_manova(s: "Session", args: str) -> None:
    from ..estimation import anova as A
    from ._util import touse
    from .summarize import weights
    p = parse_standard(args)
    o = match_options(p.options, {"sequential": 3, "continuous": 4, "noconstant": 4, "partial": 4,
                                  "dropemptycells": 4}) if p.options.strip() else {}
    if p.exp is None:
        raise StataError(198, "invalid syntax")
    deps = numeric_names(s, p.varlist, 2)
    continuous = set(s.expand_varlist(o["continuous"])) if isinstance(o.get("continuous"), str) else set()
    segments = A._parse_terms(s, p.exp, continuous)
    terms = [t_ for seg in segments for t_ in seg]
    if not terms:
        raise StataError(102, "too few variables specified")
    allvars = sorted({n for t_ in terms for n in t_.variables})
    ds = s.data
    mask = touse(s, p)
    for n_ in deps + allvars:
        if ds.get(n_).is_string:
            raise StataError(109, "type mismatch")
        mask &= ds.get(n_).data < SYS
    w, wtype, mask = weights(s, p, mask, ("aweight", "fweight"))
    if not mask.any():
        raise StataError(2000, "no observations")
    Y = np.column_stack([ds.get(d_).data[mask] for d_ in deps]).astype(float)
    ww = np.ones(Y.shape[0]) if w is None else w[mask].astype(float)
    if wtype == "aweight":
        ww = ww * len(ww) / ww.sum()
    data = {}
    for n_ in allvars:
        x = ds.get(n_).data[mask].astype(float)
        data[n_] = (x, sorted(set(x.tolist())))
    blocks = [A._term_matrix(t_, data) for t_ in terms]
    sequential = bool(o.get("sequential"))
    one = np.ones((Y.shape[0], 1))
    Xfull = np.column_stack([one] + blocks)
    E, rank = _resid_sscp(Xfull, Y, ww)
    N = float(ww.sum())
    ybar = (ww[:, None] * Y).sum(0) / N
    Yc = Y - ybar
    Ttot = (Yc * ww[:, None]).T @ Yc
    Hm = Ttot - E
    q = len(deps)
    df_m, df_r = rank - 1, N - rank
    out = s.output
    out.write(f"\n{'':>23}Number of obs = {int(N):>10}\n\n", "text")
    out.write(f"{'':>23}W = Wilks' lambda      L = Lawley-Hotelling trace\n", "text")
    out.write(f"{'':>23}P = Pillai's trace     R = Roy's largest root\n\n", "text")
    out.write(f"{'Source':>20} | Statistic        df    F(df1,     df2) =   F   Prob>F\n", "text")
    rule = f"{'':>10}{'-' * 11}+{'-' * 55}"
    sep = f"{'':>21}|{'-' * 55}"
    out.write(rule + "\n", "text")
    res = {}
    lam_m = T.eig_ratio(E, Hm)
    res["model"] = (T.mv_statistics(lam_m, q, df_m, df_r), df_m, Hm)
    rows = [("Model", res["model"])]
    prev_E, prev_rank = _resid_sscp(one, Y, ww)
    for k_, t_ in enumerate(terms):
        if sequential:
            Ek, rk = _resid_sscp(np.column_stack([one] + blocks[:k_ + 1]), Y, ww)
            H, dfh = prev_E - Ek, rk - prev_rank
            prev_E, prev_rank = Ek, rk
        else:
            Er, rk = _resid_sscp(np.column_stack([one] + [b for j, b in enumerate(blocks) if j != k_]), Y, ww)
            H, dfh = Er - E, rank - rk
        lam = T.eig_ratio(E, H) if dfh > 0 else np.zeros(0)
        rows.append((t_.name, (T.mv_statistics(lam, q, dfh, df_r) if dfh > 0 else [], dfh, H)))
    index = {id(t_): k_ for k_, t_ in enumerate(terms)}
    groups = _error_groups(segments, index)
    resid = f"{'Residual':>20} |{int(df_r):>20}\n"
    if len(terms) > 1:
        _stat_rows(out, "Model", rows[0][1][0], rows[0][1][1])
        out.write(sep + "\n" + resid + rule + "\n", "text")
    tested = {}
    for members, err in groups:
        for k_ in members:
            label, (_, dfh, H) = rows[k_ + 1]
            if err is None:
                stats = rows[k_ + 1][1][0]
            else:
                He, dfe = rows[err + 1][1][2], rows[err + 1][1][1]
                stats = T.mv_statistics(T.eig_ratio(He, H), q, dfh, dfe) if dfh > 0 and dfe > 0 else []
            tested[k_] = (stats, err)
            _stat_rows(out, label, stats, dfh)
            out.write(sep + "\n", "text")
        if err is None:
            out.write(resid + rule + "\n", "text")
        else:
            _label_df_row(out, terms[err].name, rows[err + 1][1][1])
            out.write(rule + "\n", "text")
    out.write(f"{'Total':>20} |{int(N - 1):>20}\n{'':>10}{'-' * 67}\n", "text")
    out.write(f"{'':>23}e = exact, a = approximate, u = upper bound on F\n", "text")
    # e() in the order of Stata 14's ereturn list (compat 1929)
    Xo = np.column_stack([_overparam_block(t_, data) for t_ in terms] + [one])
    Xw = Xo * ww[:, None]
    G = np.linalg.pinv(Xo.T @ Xw)
    B = G @ Xw.T @ Y
    fv = _design_fvnames(terms, data)
    bm = matrix(B.T.reshape(1, -1), ["y1"], fv * q)
    bm.coleq = [d_ for d_ in deps for _ in fv]
    Vm = matrix(np.kron(E / df_r, G), fv * q, fv * q)
    Vm.roweq = Vm.coleq = bm.coleq
    mats = [("b", bm), ("V", Vm)]
    for j in range(len(terms), 0, -1):
        stats, err = tested.get(j - 1, (rows[j][1][0], None))
        dfh, H = rows[j][1][1], rows[j][1][2]
        He, dfe = (E, df_r) if err is None else (rows[err + 1][1][2], rows[err + 1][1][1])
        lam = T.eig_ratio(He, H) if dfh > 0 and dfe > 0 else np.zeros(0)
        eig, aux, stat = _test_mats(lam, stats, q, dfh, dfe)
        mats += [(f"eigvals_{j}", eig), (f"aux_{j}", aux), (f"stat_{j}", stat), (f"H_{j}", matrix(H, deps, deps))]
    eig, aux, stat = _test_mats(lam_m, res["model"][0], q, df_m, df_r)
    mats += [("eigvals_m", eig), ("aux_m", aux), ("stat_m", stat), ("H_m", matrix(Hm, deps, deps)),
             ("E", matrix(E, deps, deps)), ("B", matrix(B, fv, deps)), ("xpxinv", matrix(G, fv, fv))]
    mats = [(k_, v_) for k_, v_ in mats if v_ is not None]
    scalars = [("N", N), ("k", float(rank)), ("rank", float(rank * q)), ("df_m", df_m), ("df_r", df_r),
               ("k_eq", float(q))] + [(f"df_{j}", float(rows[j][1][1])) for j in range(1, len(terms) + 1)]
    # the univariate regression of each dependent variable (VERIFY: the number formats)
    rss, tss = np.diag(E), np.diag(Ttot)
    F_ = (tss - rss) / df_m / (rss / df_r)
    pF = [float(T.st.f.sf(f_, df_m, df_r)) for f_ in F_]
    indep = list(dict.fromkeys(n_ for t_ in terms for n_ in t_.variables))
    macros = [("cmdline", "manova " + args.strip()), ("marginsnotok", "stdp Residuals STDDp"),
              ("marginsdefault", " ".join(f"predict(xb equation({d_}))" for d_ in deps)),
              ("p_F", " ".join(f"{x:.4f}" for x in pF)), ("F", " ".join(gnum(x).strip() for x in F_)),
              ("r2", " ".join(f"{1 - r_ / t_:.4f}" for r_, t_ in zip(rss, tss))),
              ("rmse", " ".join(gnum(math.sqrt(r_ / df_r)).strip() for r_ in rss)),
              ("cmd", "manova"), ("predict", "manova_p"), ("estat_cmd", "manova_estat"), ("properties", "b V"),
              ("depvar", " ".join(deps)), ("indepvars", " ".join(indep)), ("wtype", wtype)]
    tm = []
    for j, t_ in enumerate(terms, start=1):
        tm.append((f"term_{j}", t_.name))
        err = tested.get(j - 1, (None, None))[1]
        if err is not None:
            tm.append((f"errorterm_{j}", terms[err].name))
    macros += tm[::-1]
    post_e(s, "manova", scalars, macros, mats, mask, int(N),
           extra={"manova": {"Y": Y, "w": ww, "blocks": blocks, "terms": [t_.name for t_ in terms],
                             "order": _design_order(terms, data),
                             "fvnames": _design_fvnames(terms, data),
                             "Xo": Xo,
                             "ref_blocks": [_reference_block(t_, data) for t_ in terms],
                             "nested": [bool(t_.nest) for t_ in terms],
                             "E": E, "rank": rank, "df_r": df_r, "deps": deps}})


def _error_groups(segments, index) -> list[tuple[list[int], int | None]]:
    """Which terms are tested against which error term ([MV] manova, nested designs): the terms
    of a segment are tested against the first term of the next segment (the residual when the
    next segment is empty, a trailing slash). A term that follows a slash is an error term and is
    tested only when it is itself followed by a slash, that is, when it is alone in its segment
    (VERIFY beyond the manual's two examples). Returns [(term indices, error index or None)]."""
    groups = []
    last = len(segments) - 1
    for i, seg in enumerate(segments):
        ids = [index[id(t_)] for t_ in seg]
        if i == last:
            members = ids if i == 0 else ids[1:]
            err = None
        else:
            members = ids if i == 0 or len(ids) == 1 else ids[1:]
            nxt = segments[i + 1]
            err = index[id(nxt[0])] if nxt else None
        if members:
            groups.append((members, err))
    return groups


def _wrap_label(label: str, width: int = 20) -> list[str]:
    """The lines of a source label: terms joined by spaces while they fit in the column ("a b",
    compat 1917), a term longer than the column broken after a | or # (the manual's
    "associate|store|" + "video"; VERIFY for other shapes)."""
    pieces = []
    for term in label.split() or [""]:
        while len(term) > width:
            cut = max(term.rfind("|", 0, width), term.rfind("#", 0, width))
            if cut <= 0:
                cut = width - 1
            pieces.append((term[:cut + 1], True))
            term = term[cut + 1:]
        pieces.append((term, False))
    lines: list[str] = []
    joinable = False
    for piece, broken in pieces:
        if lines and joinable and len(lines[-1]) + 1 + len(piece) <= width:
            lines[-1] += " " + piece
        else:
            lines.append(piece)
        joinable = not broken
    return lines


def _label_df_row(out, label: str, df: float) -> None:
    pieces = _wrap_label(label)
    out.write(f"{pieces[0]:>20} |{int(df):>20}\n", "text")
    for piece in pieces[1:]:
        out.write(f"{piece:>20} |\n", "text")


def _overparam_block(t_, data: dict) -> np.ndarray:
    """The columns of a term in the overparameterized design (an indicator for every level,
    continuous variables as they are, crossed for interactions), in manovatest, showorder's
    order."""
    mats = []
    for name, cont in list(t_.comps) + [(v, False) for v in t_.nest]:
        x, levels = data[name]
        mats.append(x[:, None] if cont else np.column_stack([(x == lev).astype(float) for lev in levels]))
    out = mats[0]
    for m_ in mats[1:]:
        out = np.einsum("ij,ik->ijk", out, m_).reshape(len(out), -1)
    return out


def _design_fvnames(terms, data: dict) -> list[str]:
    """Factor-variable names of the overparameterized columns, as the test() constraints are
    written: 1.group, 1.a#2.b, x, _cons ([MV] manova postestimation, example 2)."""
    out = []
    for t_ in terms:
        labs = [""]
        for name, cont in list(t_.comps) + [(v, False) for v in t_.nest]:
            parts = [name] if cont else [f"{lev:g}.{name}" for lev in data[name][1]]
            labs = [f"{a}#{b}" if a else b for a in labs for b in parts]
        out += labs
    return out + ["_cons"]


def _design_order(terms, data: dict) -> list[str]:
    """Labels of the columns of the overparameterized design matrix, as manovatest, showorder
    lists them: (a==0), ..., (a==0)*(b==0), ..., _cons (compat 1917; nested terms VERIFY)."""
    out = []
    for t_ in terms:
        labs = [""]
        for name, cont in list(t_.comps) + [(v, False) for v in t_.nest]:
            parts = [name] if cont else [f"({name}=={lev:g})" for lev in data[name][1]]
            labs = [f"{a}*{b}" if a else b for a in labs for b in parts]
        out += labs
    return out + ["_cons"]


def _reference_block(t_, data: dict) -> np.ndarray:
    """The columns of a term in the base-level coding (indicators of every level but the
    first; continuous variables as they are), crossed for interactions."""
    mats = []
    for name, cont in list(t_.comps) + [(v, False) for v in t_.nest]:
        x, levels = data[name]
        mats.append(x[:, None] if cont else np.column_stack([(x == lev).astype(float) for lev in levels[1:]]))
    out = mats[0]
    for m_ in mats[1:]:
        out = np.einsum("ij,ik->ijk", out, m_).reshape(len(out), -1)
    return out


def _manova_terms(text: str, names: list[str]) -> list[int]:
    """Indices of the terms of the last manova named in text (as written in e(term_#))."""
    out = []
    norm = {n_.replace(" ", ""): i for i, n_ in enumerate(names)}
    for w_ in text.split():
        key = w_.replace(" ", "")
        if key not in norm:
            raise StataError(111, f"manova variable {w_} not found")
        out.append(norm[key])
    return out


def _joint_hypothesis(Y, ww, m, idx):
    """SSCP matrix and degrees of freedom of the test of the terms idx: one term is the partial
    test (the model without the term, sum-to-zero coding); several terms drop their columns
    together from the model in the base-level coding, which is what Stata's joint test
    reproduces (compat 1917; VERIFY beyond two main effects next to their interaction)."""
    # nested terms keep their within-cell coding: the base-level coding crosses the nesting
    # variables as if they were factors (the pooled error of [MV] manova postestimation, example 4)
    blocks = m["blocks"] if len(idx) == 1 else [b if nested else r for b, r, nested in
                                                 zip(m["blocks"], m["ref_blocks"], m["nested"])]
    one = np.ones((Y.shape[0], 1))
    E, rank = _resid_sscp(np.column_stack([one] + blocks), Y, ww)
    Er, rk = _resid_sscp(np.column_stack([one] + [b for k, b in enumerate(blocks) if k not in idx]), Y, ww)
    return Er - E, rank - rk


def _test_hypothesis(Y, ww, Xo, C):
    """SSCP matrix and degrees of freedom of the hypothesis C beta = 0 on the overparameterized
    design ([MV] manova, Methods and formulas): H = (CB)'{C (X'X)^- C'}^-1 (CB), B = (X'X)^- X'Y."""
    if C.shape[1] != Xo.shape[1]:
        raise StataError(503, "conformability error")
    Xw = Xo * ww[:, None]
    G = np.linalg.pinv(Xo.T @ Xw)
    B = G @ Xw.T @ Y
    CB = C @ B
    M_ = C @ G @ C.T
    return CB.T @ np.linalg.pinv(M_) @ CB, int(np.linalg.matrix_rank(M_))


def _ytransform_lines(Mt, deps) -> list[str]:
    """The rows of ytransform() written as combinations of the dependent variables (VERIFY for
    coefficients other than 1 and -1)."""
    lines = []
    for row in Mt:
        text = ""
        for c, d in zip(row, deps):
            if c == 0:
                continue
            mag = "" if abs(c) == 1 else f"{abs(c):g}*"
            if not text:
                text = ("-" if c < 0 else "") + mag + d
            else:
                text += (" - " if c < 0 else " + ") + mag + d
        lines.append(text)
    return lines


@command("manovatest", "manovatest")
def cmd_manovatest(s: "Session", args: str) -> None:
    """manovatest term [term ...] [/ term [term ...]] [, ytransform(matname) showorder] after
    manova ([MV] manova postestimation): the multivariate tests of the terms, with the residual
    or the terms after the slash as the error. test() is not available yet."""
    from .matrix import get_matrix
    est = getattr(s, "_estimates", None)
    if est is None or est.cmd != "manova":
        raise StataError(301, "last estimates not found")
    c = _find_comma(args)
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    o = match_options(opt, {"ytransform": 2, "test": 4, "showorder": 5}) if opt.strip() else {}
    m = est.extra["manova"]
    out = s.output
    if o.get("showorder"):
        out.write("\n Order of columns in the design matrix\n", "text")
        for i, lab in enumerate(m["order"], start=1):
            out.write(f"{i:>7}: {lab}\n", "text")
        if not body.strip():
            return
    num, _, den = body.partition("/")
    Y, ww = m["Y"], m["w"]
    q = len(m["deps"])
    Cm = None
    if isinstance(o.get("test"), str):
        if body.strip():
            raise StataError(198, "invalid syntax")
        Cm = get_matrix(s, o["test"].strip()).data
        H, dfh = _test_hypothesis(Y, ww, m["Xo"], Cm)
        num, den = "manovatest", ""
    elif not num.strip():
        raise StataError(198, "invalid syntax")
    else:
        H, dfh = _joint_hypothesis(Y, ww, m, _manova_terms(num, m["terms"]))
    if den.strip():
        Eh, dfe = _joint_hypothesis(Y, ww, m, _manova_terms(den, m["terms"]))
        elab = " ".join(den.split())
    else:
        Eh, dfe = m["E"], m["df_r"]
        elab = "Residual"
    Mt = None
    if isinstance(o.get("ytransform"), str):
        Mt = get_matrix(s, o["ytransform"].strip()).data
        if Mt.shape[1] != q:
            raise StataError(503, "conformability error")
        H, Eh = Mt @ H @ Mt.T, Mt @ Eh @ Mt.T
        q = Mt.shape[0]
    if dfe < q:
        # an error term with fewer degrees of freedom than dependent variables (compat 1917)
        out.write("matrix not positive definite\n", "text")
        raise StataError(459, "insufficient denominator degrees of freedom")               # VERIFY rc
    if np.linalg.matrix_rank(Eh) < Eh.shape[0]:
        raise StataError(506, "error matrix is singular")                                   # VERIFY
    lam = T.eig_ratio(Eh, H) if dfh > 0 else np.zeros(0)
    stats = T.mv_statistics(lam, q, dfh, dfe) if dfh > 0 else []
    if Mt is not None:
        lines = _ytransform_lines(Mt, m["deps"])
        out.write(f"\n Transformation{'s' if len(lines) > 1 else ''} of the dependent variables\n", "text")
        for i, line in enumerate(lines, start=1):
            out.write(f" ({i})    {line}\n", "text")
    if Cm is not None:
        lines = _ytransform_lines(Cm, m["fvnames"])
        out.write(f"\n Test constraint{'s' if len(lines) > 1 else ''}\n", "text")
        for i, line in enumerate(lines, start=1):
            out.write(f" ({i})    {line} = 0\n", "text")
    out.write(f"\n{'':>23}W = Wilks' lambda      L = Lawley-Hotelling trace\n", "text")
    out.write(f"{'':>23}P = Pillai's trace     R = Roy's largest root\n\n", "text")
    out.write(f"{'Source':>20} | Statistic        df    F(df1,     df2) =   F   Prob>F\n", "text")
    rule = f"{'':>10}{'-' * 11}+{'-' * 55}"
    out.write(rule + "\n", "text")
    _stat_rows(out, " ".join(num.split()), stats, dfh)
    out.write(f"{'':>21}|{'-' * 55}\n", "text")
    _label_df_row(out, elab, dfe)
    out.write(f"{'':>10}{'-' * 67}\n", "text")
    out.write(f"{'':>23}e = exact, a = approximate, u = upper bound on F\n", "text")
    # r() as Stata lists it (return list shows the last stored first): df, df_r, then eigvals,
    # aux, stat, E, H; with an error term after the slash aux, stat, H, E, eigvals (compat 1917, 1929)
    eig, aux, stat = _test_mats(lam, stats, q, dfh, dfe)
    if den.strip():
        s.r = {"eigvals": eig, "E": matrix(Eh), "H": matrix(H), "stat": stat, "aux": aux}
    else:
        s.r = {"H": matrix(H), "E": matrix(Eh), "stat": stat, "aux": aux, "eigvals": eig}
    s.r.update({"df_r": float(dfe), "df": float(dfh)})
    s.r = {k_: v_ for k_, v_ in s.r.items() if v_ is not None}


def _find_comma(text: str) -> int:
    depth = 0
    for i, ch in enumerate(text):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == "," and depth == 0:
            return i
    return -1

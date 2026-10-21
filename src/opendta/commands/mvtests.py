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
from .mvcommon import chi2_p, g as gnum, matrix, numeric_names, post_e, sample
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _fp(p: float) -> str:
    return f"{p:.4f}" if np.isfinite(p) else "."


def _summary_table(out, ds, names, X, w=None) -> None:
    out.write("\n    Variable |       Obs        Mean    Std. Dev.       Min        Max\n", "text")
    out.write("-------------+" + "-" * 56 + "\n", "text")
    for j, n_ in enumerate(names):
        x = X[:, j]
        n = len(x) if w is None else float(w.sum())
        m = float(np.average(x, weights=w))
        sd = math.sqrt(float(np.average((x - m) ** 2, weights=w)) * n / (n - 1)) if n > 1 else float("nan")
        out.write(f"{n_[:12]:>12} |{int(n):>10}{gnum(m, '%9.0g'):>12}{gnum(sd, '%9.0g'):>12}"
                  f"{gnum(float(x.min()), '%9.0g'):>11}{gnum(float(x.max()), '%9.0g'):>11}\n", "text")


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
            _summary_table(out, ds, names, X, w)
        r = T.hotelling_one(mean, S, n)
        out.write(f"\n1-group Hotelling's T-squared = {r['T2']:.4f}\n", "text")
        out.write(f"F test statistic: (({int(n)}-{k})/({int(n)}-1)({k})) x {r['T2']:.4f} = {r['F']:.4f}\n", "text")
        out.write("\nH0: Vector of means is equal to a vector of zeros\n", "text")
        out.write(f"{'F(' + str(k) + ', ' + gnum(r['df2']) + ') =':>25}{r['F']:>10.4f}\n", "text")
        out.write(f"{'Prob > F(' + str(k) + ', ' + gnum(r['df2']) + ') =':>25}{r['p']:>10.4f}\n", "text")
        s.r = {"T2": r["T2"], "F": r["F"], "df": float(k), "df_r": r["df2"], "p": r["p"]}
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
            out.write(f"\n-> {by} = {lab}\n", "text")
            _summary_table(out, ds, names, X[sel], ww)
    (m1, S1, n1, _), (m2, S2, n2, _) = parts
    r = T.hotelling_two(m1, S1, n1, m2, S2, n2)
    out.write(f"\n2-group Hotelling's T-squared = {r['T2']:.4f}\n", "text")
    out.write(f"F test statistic: (({int(n1)}+{int(n2)}-{k}-1)/({int(n1)}+{int(n2)}-2)({k})) x {r['T2']:.4f} = {r['F']:.4f}\n", "text")
    out.write("\nH0: Vectors of means are equal for the two groups\n", "text")
    d2 = gnum(r["df2"])
    out.write(f"{'F(' + str(k) + ', ' + d2 + ') =':>25}{r['F']:>10.4f}\n", "text")
    out.write(f"{'Prob > F(' + str(k) + ', ' + d2 + ') =':>25}{r['p']:>10.4f}\n", "text")
    s.r = {"T2": r["T2"], "F": r["F"], "df": float(k), "df_r": r["df2"], "p": r["p"]}


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
    from ..estimation.results import CoefRow, Estimates, coef_table, post, rows_from_columns
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
        rows.append(CoefRow("header", d_))
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
    scalars = [("N", N), ("k", kk * q), ("k_eq", q), ("k_dv", q), ("df_r", df_r), ("chi2", SYS)]
    for i, (d_, ob, parms, rmse, r2, F, pF) in enumerate(eq_stats, start=1):
        scalars += [(f"rmse_{i}", rmse), (f"r2_{i}", r2), (f"F_{i}", F), (f"df_m{i}", float(df_m))]
    macros = [("cmdline", "mvreg " + args.strip()), ("cmd", "mvreg"), ("depvar", " ".join(deps)),
              ("eqnames", " ".join(deps)), ("title", "Multivariate regression"), ("predict", "mvreg_p"),
              ("marginsnotok", "_ALL"), ("properties", "b V"), ("wtype", wtype)]
    post(s, est, scalars, macros, mask, extra_mats=[("Sigma", matrix(Sigma, deps, deps))])
    est.no_b = False
    _display_mvreg(s, est, eq_stats, o, deps, N)


def _display_mvreg(s, est, eq_stats, o, deps, N) -> None:
    from ..estimation.results import coef_table
    out = s.output
    out.write("\n", "text")
    out.write('Equation          Obs  Parms        RMSE    "R-sq"          F        P\n', "text")
    out.write("-" * 70 + "\n", "text")
    for d_, ob, parms, rmse, r2, F, pF in eq_stats:
        out.write(f"{d_[:14]:<14}{int(ob):>7}{parms:>7}{rmse:>12.7g}{r2:>10.4f}{F:>11.2f}{_fp(pF):>9}\n", "text")
    out.write("\n", "text")
    if o.get("corr"):
        q = len(deps)
        Sg = est.extra["Sigma"]
        R = Sg / np.sqrt(np.outer(np.diag(Sg), np.diag(Sg)))
        chi2 = float(N * sum(R[i, j] ** 2 for i in range(q) for j in range(i)))
        df = q * (q - 1) / 2.0
        out.write("Correlation matrix of residuals:\n\n", "text")
        out.write(f"{'':>13}" + "".join(f"{d_[:8]:>10}" for d_ in deps) + "\n", "text")
        for i, d_ in enumerate(deps):
            out.write(f"{d_[:12]:>12} " + "".join(f"{R[i, j]:>10.4f}" for j in range(i + 1)) + "\n", "text")
        out.write(f"\nBreusch-Pagan test of independence: chi2({int(df)}) = {chi2:8.3f}, Pr = {T.st.chi2.sf(chi2, df):.4f}\n\n"
                  if q > 1 else "", "text")
    coef_table(s, "", est.rows, stat="t", df=est.df_r, level=est.level)


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


def _stat_rows(out, label: str, stats: list[dict], df_h: float) -> None:
    for i, r in enumerate(stats):
        src = label if i == 0 else ""
        dfc = f"{int(df_h):>5}" if i == 0 else f"{'':>5}"
        out.write(f"{src:>20} | {r['name']} {r['stat']:>8.4f} {dfc} {r['df1']:>7.1f} {r['df2']:>7.1f} "
                  f"{r['F']:>8.2f} {_fp(r['p']):>7} {r['flag']}\n", "text")


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
    if len(segments) > 1:
        raise StataError(198, "error terms (/) are not available yet in OpenDTA")
    terms = segments[0]
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
    out.write(f"\n{'Number of obs =':>40}{int(N):>9}\n\n", "text")
    out.write(f"{'':>26}W = Wilks' lambda      L = Lawley-Hotelling trace\n", "text")
    out.write(f"{'':>26}P = Pillai's trace     R = Roy's largest root\n\n", "text")
    out.write(f"{'Source':>20} | Statistic     df   F(df1,    df2) =   F   Prob>F\n", "text")
    rule = f"{'-' * 21}+{'-' * 52}"
    out.write(rule + "\n", "text")
    res = {}
    lam_m = T.eig_ratio(E, Hm)
    res["model"] = (T.mv_statistics(lam_m, q, df_m, df_r), df_m, Hm)
    rows = [("model", res["model"])]
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
    for label, (stats, dfh, H) in rows:
        _stat_rows(out, label, stats, dfh)
        out.write(rule + "\n", "text")
    out.write(f"{'Residual':>20} | {'':>9} {int(df_r):>5}\n{rule}\n", "text")
    out.write(f"{'Total':>20} | {'':>9} {int(N - 1):>5}\n{'-' * 74}\n", "text")
    out.write(f"{'':>20}   e = exact, a = approximate, u = upper bound on F\n", "text")
    names_cols = ["statistic", "df1", "df2", "F", "p"]
    mats = [("E", matrix(E, deps, deps)), ("H_m", matrix(Hm, deps, deps)),
            ("eigvals_m", matrix(lam_m.reshape(1, -1), ["eigvals"], [f"c{i + 1}" for i in range(len(lam_m))])),
            ("stat_m", matrix(np.array([[r["stat"], r["df1"], r["df2"], r["F"], r["p"]] for r in res["model"][0]]),
                              ["W", "P", "L", "R"], names_cols))]
    for j, (label, (stats, dfh, H)) in enumerate(rows[1:], start=1):
        mats.append((f"H_{j}", matrix(H, deps, deps)))
    scalars = [("N", N), ("k_dv", q), ("df_m", df_m), ("df_r", df_r), ("rank", rank)]
    macros = [("cmdline", "manova " + args.strip()), ("cmd", "manova"), ("depvar", " ".join(deps)),
              ("title", "Multivariate model"), ("model", "ols"), ("properties", "nob noV"), ("wtype", wtype)]
    macros += [(f"term_{j}", t_.name) for j, t_ in enumerate(terms, start=1)]
    post_e(s, "manova", scalars, macros, mats, mask, int(N))


# ---------------------------------------------------------------------------
# canon
# ---------------------------------------------------------------------------

def _split_sets(text: str) -> tuple[str, str, str]:
    t = text.strip()
    if not t.startswith("("):
        raise StataError(198, "invalid syntax")
    depth, parts, cur, i = 0, [], "", 0
    while i < len(t):
        ch = t[i]
        if ch == "(":
            depth += 1
            if depth == 1:
                cur = ""
                i += 1
                continue
        elif ch == ")":
            depth -= 1
            if depth == 0:
                parts.append(cur)
                i += 1
                continue
        if depth >= 1:
            cur += ch
        elif depth == 0 and ch not in " \t" and ch != ",":
            break
        i += 1
    rest = t[i:]
    if len(parts) != 2:
        raise StataError(198, "invalid syntax")
    return parts[0], parts[1], rest


@command("canon", "canon")
def cmd_canon(s: "Session", args: str) -> None:
    c = _find_comma(args)
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    set1, set2, rest = _split_sets(body)
    p = parse_standard("x " + rest + ("," + opt if opt else ""))
    o = match_options(p.options, {"lc": 2, "noconstant": 3, "stdcoef": 3, "test": 3, "level": 1}) if p.options.strip() else {}
    n1 = numeric_names(s, set1, 1)
    n2 = numeric_names(s, set2, 1)
    X, w, wtype, mask = sample(s, p, n1 + n2)
    p1, p2 = len(n1), len(n2)
    noconstant = bool(o.get("noconstant"))
    if noconstant:
        n_ = X.shape[0]
        sw = float(n_ if w is None else w.sum())
        W = np.ones(n_) if w is None else w
        C = (X * W[:, None]).T @ X / (sw - 1.0)
        means = np.zeros(p1 + p2)
    else:
        means, C, _, sw = mv.weighted_moments(X, w)
    N = int(round(sw))
    S11, S22, S12 = C[:p1, :p1], C[p1:, p1:], C[:p1, p1:]
    cc = T.canonical(S11, S22, S12)
    r, A1, A2 = cc["r"], cc["A1"], cc["A2"]
    # sign convention: the first element of each coefficient column is positive (VERIFY)
    for j in range(len(r)):
        if A1[0, j] < 0:
            A1[:, j] *= -1
            A2[:, j] *= -1
    k = len(r)
    if o.get("stdcoef"):
        A1 = A1 * np.sqrt(np.diag(S11))[:, None]
        A2 = A2 * np.sqrt(np.diag(S22))[:, None]
    out = s.output
    out.write(f"\nCanonical correlation analysis{'':>20}Number of obs    ={N:>11}\n", "text")
    cn = [str(j + 1) for j in range(k)]
    for title, nm, A_ in (("first", n1, A1), ("second", n2, A2)):
        out.write(f"\n{'Standardized' if o.get('stdcoef') else 'Raw'} coefficients for the {title} variable set\n\n", "text")
        out.write(f"{'':>16} |" + "".join(f"{c_:>10}" for c_ in cn) + "\n", "text")
        out.write("-" * 17 + "+" + "-" * (10 * k + 1) + "\n", "text")
        for i, n_ in enumerate(nm):
            out.write(f"{n_[:15]:>16} |" + "".join(f"{A_[i, j]:>10.4f}" for j in range(k)) + "\n", "text")
    out.write("\n" + "-" * 76 + "\nCanonical correlations:\n  " + "  ".join(f"{x:.4f}" for x in r) + "\n", "text")
    out.write("-" * 76 + "\nTests of significance of all canonical correlations\n\n", "text")
    out.write(f"{'':>25}{'Statistic':>10}{'df1':>8}{'df2':>10}{'F':>10}{'Prob>F':>9}\n", "text")
    lam = r ** 2 / (1.0 - r ** 2)
    v = N - 1.0 - p2 if not noconstant else N - p2
    stats = T.mv_statistics(lam, p1, p2, v)
    names_full = {"W": "Wilks' lambda", "P": "Pillai's trace", "L": "Lawley-Hotelling trace", "R": "Roy's largest root"}
    for r_ in stats:
        out.write(f"{names_full[r_['name']]:>24} {r_['stat']:>10.4f}{r_['df1']:>8.0f}{r_['df2']:>10.1f}"
                  f"{r_['F']:>10.4f}{_fp(r_['p']):>9} {r_['flag']}\n", "text")
    out.write(f"{'':>45}e = exact, a = approximate, u = upper bound on F\n", "text")
    sc = [("N", N), ("n_lc", k), ("df_m", p1 * p2), ("df_r", v)]
    mats = [("ccorr", matrix(r.reshape(1, -1), ["ccorr"], cn)), ("coef1", matrix(A1, n1, cn)), ("coef2", matrix(A2, n2, cn))]
    macros = [("cmdline", "canon " + args.strip()), ("cmd", "canon"), ("varlist1", " ".join(n1)), ("varlist2", " ".join(n2)),
              ("title", "Canonical correlation analysis"), ("properties", "nob noV"), ("wtype", wtype)]
    post_e(s, "canon", sc, macros, mats, mask, N)


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


# ---------------------------------------------------------------------------
# mvtest
# ---------------------------------------------------------------------------

@command("mvtest", "mvtest")
def cmd_mvtest(s: "Session", args: str) -> None:
    words = args.strip().split(None, 1)
    if not words:
        raise StataError(198, "invalid syntax")
    sub = words[0].rstrip(",")
    rest = args.strip()[len(words[0]):] if len(words) > 1 else ""
    for name, mn, fn in (("means", 1, _mvtest_means), ("covariances", 3, _mvtest_cov), ("correlations", 4, _mvtest_corr),
                         ("normality", 1, _mvtest_norm)):
        if len(sub) >= mn and name.startswith(sub):
            fn(s, rest)
            return
    raise StataError(198, "invalid syntax")


def _groups(s, p, names, o, allowed=("fweight", "aweight")):
    by = str(o["by"]).strip() if isinstance(o.get("by"), str) else None
    X, w, wtype, mask = sample(s, p, names, allowed, extra=[by] if by else [])
    if by is None:
        return X, w, None, mask, None
    g = s.data.get(by).data[mask]
    return X, w, g, mask, by


def _chi_lines(out, title: str, chi2: float, df: float, p: float) -> None:
    out.write(f"\n{title}\n\n", "text")
    out.write(f"{'chi2(' + gnum(df) + ') =':>25}{chi2:>10.4f}\n{'Prob > chi2 =':>25}{p:>10.4f}\n", "text")


def _mvtest_means(s, rest: str) -> None:
    p = parse_standard(rest)
    o = match_options(p.options, {"by": 2, "zero": 4, "equals": 3, "homogeneous": 4, "heterogeneous": 6}) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 1)
    X, w, g, mask, by = _groups(s, p, names, o)
    out = s.output
    k = len(names)
    if g is None:
        mean, S, R, n = mv.weighted_moments(X, w)
        if o.get("zero") or isinstance(o.get("equals"), str):
            mu0 = np.zeros(k)
            if isinstance(o.get("equals"), str):
                from .matrix import get_matrix
                mu0 = get_matrix(s, str(o["equals"]).strip()).data.ravel()
            r = T.hotelling_one(mean, S, n, mu0)
            what = "Test that all means are equal to the hypothesized values" if not o.get("zero") else "Test that all means are equal to 0"
            df1 = k
        else:
            if k < 2:
                raise StataError(102, "too few variables specified")
            Cm = np.hstack([np.eye(k - 1), -np.ones((k - 1, 1))])
            r = T.hotelling_one(Cm @ mean, Cm @ S @ Cm.T, n)
            what = "Test that all means are equal"
            df1 = k - 1
        out.write(f"\n{what}\n\n", "text")
        out.write(f"{'Hotelling T2 =':>25}{r['T2']:>12.4f}\n", "text")
        out.write(f"{'F(' + str(df1) + ',' + gnum(r['df2']) + ') =':>25}{r['F']:>12.4f}\n", "text")
        out.write(f"{'Prob > F(' + str(df1) + ',' + gnum(r['df2']) + ') =':>25}{r['p']:>12.4f}\n", "text")
        s.r = {"T2": r["T2"], "F": r["F"], "df_m": float(df1), "df_r": r["df2"], "p_F": r["p"]}
        return
    levels = sorted(set(g.tolist()))
    if len(levels) < 2:
        raise StataError(198, "by() must define at least 2 groups")                         # VERIFY
    E = np.zeros((k, k))
    gm = np.average(X, axis=0, weights=w)
    H = np.zeros((k, k))
    N = float(X.shape[0] if w is None else w.sum())
    for lv in levels:
        sel = g == lv
        ww = np.ones(int(sel.sum())) if w is None else w[sel]
        xm = (ww[:, None] * X[sel]).sum(0) / ww.sum()
        Xc = X[sel] - xm
        E += (Xc * ww[:, None]).T @ Xc
        d = xm - gm
        H += ww.sum() * np.outer(d, d)
    dfh, dfe = len(levels) - 1, N - len(levels)
    stats = T.mv_statistics(T.eig_ratio(E, H), k, dfh, dfe)
    out.write(f"\nTest of equality of means across the groups of {by}\n\n", "text")
    out.write(f"{'':>26}{'Statistic':>10}{'F(df1,':>9}{'df2)':>9}{'=':>3}{'F':>9}{'Prob>F':>9}\n", "text")
    nm = {"W": "Wilks' lambda", "P": "Pillai's trace", "L": "Lawley-Hotelling trace", "R": "Roy's largest root"}
    for r_ in stats:
        out.write(f"{nm[r_['name']]:>25} {r_['stat']:>9.4f}{r_['df1']:>9.1f}{r_['df2']:>9.1f}   {r_['F']:>9.4f}{_fp(r_['p']):>9} {r_['flag']}\n", "text")
    s.r = {"p_F": stats[0]["p"], "F": stats[0]["F"], "df_r": stats[0]["df2"], "df_m": stats[0]["df1"]}


def _mvtest_cov(s, rest: str) -> None:
    p = parse_standard(rest)
    o = match_options(p.options, {"by": 2, "equals": 3, "diagonal": 4, "spherical": 4, "compound": 4, "block": 4}) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 2)
    for opt in ("compound", "block", "equals"):
        if o.get(opt):
            raise StataError(198, f"option {opt} is not available yet in OpenDTA")
    X, w, g, mask, by = _groups(s, p, names, o)
    out = s.output
    if g is not None:
        levels = sorted(set(g.tolist()))
        Ss, ns = [], []
        for lv in levels:
            sel = g == lv
            _, S_, _, sw = mv.weighted_moments(X[sel], None if w is None else w[sel])
            Ss.append(S_)
            ns.append(sw)
        r = T.box_m(Ss, ns)
        _chi_lines(out, f"Test of equality of covariance matrices across {len(levels)} samples", r["chi2"], r["df"], r["p"])
        s.r = {"chi2": r["chi2"], "df": r["df"], "p_chi2": r["p"]}
        return
    mean, S, R, n = mv.weighted_moments(X, w)
    if o.get("spherical"):
        r = T.lr_spherical(S, n)
        _chi_lines(out, "Test that the covariance matrix is spherical", r["chi2"], r["df"], r["p"])
    else:
        r = T.lr_diagonal(R, n)
        _chi_lines(out, "Test that the covariance matrix is diagonal", r["chi2"], r["df"], r["p"])
    s.r = {"chi2": r["chi2"], "df": r["df"], "p_chi2": r["p"]}


def _mvtest_corr(s, rest: str) -> None:
    p = parse_standard(rest)
    o = match_options(p.options, {"by": 2, "zero": 4, "equals": 3, "pairwise": 4}) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 2)
    for opt in ("equals", "pairwise"):
        if o.get(opt):
            raise StataError(198, f"option {opt} is not available yet in OpenDTA")
    X, w, g, mask, by = _groups(s, p, names, o)
    out = s.output
    if g is not None:
        levels = sorted(set(g.tolist()))
        Rs, ns = [], []
        for lv in levels:
            sel = g == lv
            _, _, R_, sw = mv.weighted_moments(X[sel], None if w is None else w[sel])
            Rs.append(R_)
            ns.append(sw)
        r = T.lr_equal_correlations(Rs, ns)
        _chi_lines(out, f"Test of equality of correlation matrices across {len(levels)} samples", r["chi2"], r["df"], r["p"])
    else:
        mean, S, R, n = mv.weighted_moments(X, w)
        r = T.lr_diagonal(R, n)
        _chi_lines(out, "Test that all correlations are zero", r["chi2"], r["df"], r["p"])
    s.r = {"chi2": r["chi2"], "df": r["df"], "p_chi2": r["p"]}


def _mvtest_norm(s, rest: str) -> None:
    p = parse_standard(rest)
    o = match_options(p.options, {"stats": 2, "univariate": 3, "bivariate": 3, "by": 2}) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 2)
    for opt in ("univariate", "bivariate", "by"):
        if o.get(opt):
            raise StataError(198, f"option {opt} is not available yet in OpenDTA")
    want = str(o.get("stats", "mardia hz")).split() if isinstance(o.get("stats"), str) else ["mardia", "hz"]
    if want == ["all"]:
        want = ["mardia", "hz"]
    if any(w_.startswith("dh") or w_.startswith("doornik") for w_ in want):
        raise StataError(198, "stats(dhansen) is not available yet in OpenDTA")
    X, w, wtype, mask = sample(s, p, names, ("fweight",))
    out = s.output
    out.write("\nTest for multivariate normality\n\n", "text")
    sres = {}
    if any(x.startswith("mar") for x in want):
        m = T.mardia(X)
        out.write(f"{'Mardia mSkewness =':>25}{m['b1']:>12.6f}   chi2({gnum(m['df_skew'])}) = {m['chi2_skew']:>9.3f}   Prob>chi2 = {m['p_skew']:.4f}\n", "text")
        out.write(f"{'Mardia mKurtosis =':>25}{m['b2']:>12.6f}   chi2(1) = {m['z_kurt'] ** 2:>11.3f}   Prob>chi2 = {m['p_kurt']:.4f}\n", "text")
        sres.update({"p_skew": m["p_skew"], "p_kurt": m["p_kurt"]})
    if any(x.startswith("hz") or x.startswith("henze") for x in want):
        h = T.henze_zirkler(X)
        out.write(f"{'Henze-Zirkler =':>25}{h['hz']:>12.6f}   chi2(1) = {h['z'] ** 2:>11.3f}   Prob>chi2 = {h['p']:.4f}\n", "text")
        sres["p_hz"] = h["p"]
    s.r = sres

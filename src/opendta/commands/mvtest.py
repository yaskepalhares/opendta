"""mvtest: tests of means, covariances, correlations and normality ([MV] mvtest).

The statistics are in core/mvtests.py."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import multivariate as mv
from ..core import mvtests as T
from ..core.errors import StataError
from ..lang.syntax import match_options, parse_standard
from .mvcommon import g as gnum, matrix, numeric_names, sample
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

# ---------------------------------------------------------------------------
# mvtest
# ---------------------------------------------------------------------------

@command("mvtest", "mvtest", byable=True)
def cmd_mvtest(s: "Session", args: str) -> None:
    words = args.strip().split(None, 1)
    if not words:
        raise StataError(198, "invalid syntax")
    sub = words[0].rstrip(",")
    rest = args.strip()[len(words[0]):] if len(words) > 1 else ""
    for name, mn, fn in (("means", 1, _mvtest_means), ("covariances", 3, _mvtest_cov), ("correlations", 4, _mvtest_corr),
                         ("normality", 1, _mvtest_norm)):
        if len(sub) >= mn and name.startswith(sub):
            from .data import _by_header
            from .summarize import by_groups
            for group, first in by_groups(s):
                if first is not None:
                    _by_header(s, first)
                s._by_group_mask = group if first is not None else None
                try:
                    fn(s, rest)
                finally:
                    s._by_group_mask = None
            return
    raise StataError(198, "invalid syntax")


def _groups(s, p, names, o, allowed=("fweight", "aweight")):
    """The sample and, with by(groupvars), the group of each observation (one code per distinct
    combination of values). Observations with missing groupvars are dropped unless `missing`."""
    by = str(o["by"]).strip() if isinstance(o.get("by"), str) else None
    if by is None:
        X, w, wtype, mask = sample(s, p, names, allowed)
        return X, w, None, mask, None
    byvars = s.expand_varlist(by)
    X, w, wtype, mask = sample(s, p, names, allowed, extra=[] if o.get("missing") else byvars)
    cols = np.column_stack([s.data.get(v).data[mask].astype(float) for v in byvars])
    _, g = np.unique(cols, axis=0, return_inverse=True)
    return X, w, g.ravel(), mask, by


def _chi_lines(out, title: str, label: str, chi2: float, df: float, p: float) -> None:
    out.write(f"\n{title}\n\n", "text")
    out.write(f"{label + '(' + gnum(df) + ') =':>24}{chi2:>10.2f}\n{'Prob > chi2 =':>24}{p:>10.4f}\n", "text")


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
            what = "Test that all means are the same"
            df1 = k - 1
        out.write(f"\n{what}\n\n", "text")
        out.write(f"{'Hotelling T2 =':>24}{r['T2']:>10.2f}\n", "text")
        out.write(f"{'Hotelling F(' + str(df1) + ',' + gnum(r['df2']) + ') =':>24}{r['F']:>10.2f}\n", "text")
        out.write(f"{'Prob > F =':>24}{r['p']:>10.4f}\n", "text")
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
    out.write(f"\nTest for equality of {len(levels)} group means, assuming homogeneity\n\n", "text")
    out.write(f"{'':>28}| Statistic   F(df1,     df2)   =  F      Prob>F\n", "text")
    out.write(f"{'':>5}{'-' * 23}+{'-' * 47}\n", "text")
    nm = {"W": "Wilks' lambda", "P": "Pillai's trace", "L": "Lawley-Hotelling trace", "R": "Roy's largest root"}
    for r_ in stats:
        out.write(f"{nm[r_['name']]:>27} |{r_['stat']:>10.4f}{r_['df1']:>8.1f}{r_['df2']:>10.1f}{r_['F']:>10.2f}"
                  f"{r_['p']:>9.4f} {r_['flag']}\n", "text")
    out.write(f"{'':>5}{'-' * 71}\n{'':>28}e = exact, a = approximate, u = upper bound on F\n", "text")
    s.r = {"p_F": stats[0]["p"], "F": stats[0]["F"], "df_r": stats[0]["df2"], "df_m": stats[0]["df1"]}


def _mvtest_cov(s, rest: str) -> None:
    p = parse_standard(rest)
    o = match_options(p.options, {"by": 2, "missing": 4, "equals": 3, "diagonal": 4, "spherical": 4, "compound": 4,
                                  "block": 4}) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 2)
    hyp = [h for h in ("diagonal", "spherical", "compound", "equals", "block") if o.get(h)]
    if len(hyp) > 1:
        raise StataError(198, f"options {hyp[0]} and {hyp[1]} may not be combined")       # VERIFY
    if o.get("by") and hyp:
        raise StataError(198, f"option {hyp[0]} not allowed with by()")                   # VERIFY
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
        out.write(f"\nTest of equality of covariance matrices across {len(levels)} samples\n\n", "text")
        out.write(f"{'Modified LR chi2 =':>24}{r['M']:>10.5f}\n", "text")
        out.write(f"{'Box F(' + gnum(r['df']) + ',' + gnum(r['df2']) + ') =':>24}{r['F']:>10.2f}"
                  f"{'Prob > F =':>16}{r['p_F']:>8.4f}\n", "text")
        out.write(f"{'Box chi2(' + gnum(r['df']) + ') =':>24}{r['chi2']:>10.2f}{'Prob > chi2 =':>16}{r['p']:>8.4f}\n", "text")
        # names as the manual's Stored results list them; the order of return list is VERIFY
        s.r = {"chi2type": "Box", "p_F_Box": r["p_F"], "df_r_Box": r["df2"], "df_m_Box": r["df"], "F_Box": r["F"],
               "p_chi2": r["p"], "df": r["df"], "chi2": r["chi2"]}
        return
    mean, S, R, n = mv.weighted_moments(X, w)
    SN = S * (n - 1.0) / n
    if o.get("spherical"):
        r = T.lr_spherical(S, n)
        title = "Test that covariance matrix is spherical"
    elif o.get("compound"):
        r = T.lr_compound_cov(SN, n)
        title = "Test that covariance matrix is compound symmetric"
    elif isinstance(o.get("equals"), str):
        from .matrix import get_matrix
        mname = str(o["equals"]).strip()
        S0 = get_matrix(s, mname).data
        if S0.shape != S.shape:
            raise StataError(503, "conformability error")
        r = T.lr_equals_cov(SN, S0, n)
        title = f"Test that covariance matrix equals matrix {mname}"
    elif isinstance(o.get("block"), str):
        blocks = [[names.index(v) for v in s.expand_varlist(part)] for part in str(o["block"]).split("||")]
        used = [i for b_ in blocks for i in b_]
        if len(set(used)) != len(used):
            raise StataError(198, "variables may appear in only one block")                    # VERIFY
        rest_ = [i for i in range(len(names)) if i not in used]
        if rest_:
            blocks.append(rest_)            # VERIFY: the variables not listed form the last block
        r = T.lr_block_cov(SN, blocks, n)
        title = "Test that covariance matrix is block diagonal"
    else:
        r = T.lr_diagonal(R, n)
        title = "Test that covariance matrix is diagonal"
    _chi_lines(out, title, "Adjusted LR chi2", r["chi2"], r["df"], r["p"])
    s.r = {"chi2type": "LR", "p_chi2": r["p"], "df": r["df"], "chi2": r["chi2"]}


def _mvtest_corr(s, rest: str) -> None:
    p = parse_standard(rest)
    o = match_options(p.options, {"by": 2, "missing": 4, "compound": 4, "equals": 3}) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 2)
    if o.get("by") and (o.get("compound") or o.get("equals")):
        raise StataError(198, "options compound and equals() not allowed with by()")      # VERIFY
    if o.get("compound") and o.get("equals"):
        raise StataError(198, "options compound and equals() may not be combined")        # VERIFY
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
        r = T.jennrich_equal(Rs, ns)
        _chi_lines(out, "Test of equality of correlation matrices across samples", "Jennrich chi2", r["chi2"], r["df"], r["p"])
        typ = "Wald"
    else:
        mean, S, R, n = mv.weighted_moments(X, w)
        if isinstance(o.get("equals"), str):
            from .matrix import get_matrix
            mname = str(o["equals"]).strip()
            C = get_matrix(s, mname).data
            if C.shape != R.shape:
                raise StataError(503, "conformability error")
            sd = np.sqrt(np.diag(C))
            r = T.jennrich_one(R, C / np.outer(sd, sd), n)
            _chi_lines(out, f"Test that correlation matrix equals specified pattern {mname}", "Jennrich chi2",
                       r["chi2"], r["df"], r["p"])
            typ = "Wald"
        else:
            r = T.lawley(R, n)
            _chi_lines(out, "Test that correlation matrix is compound symmetric (all correlations equal)",
                       "Lawley chi2", r["chi2"], r["df"], r["p"])
            typ = "Lawley"                                                                     # VERIFY
    s.r = {"chi2type": typ, "p_chi2": r["p"], "df": r["df"], "chi2": r["chi2"]}


_NORM_STATS = ("dhansen", "hzirkler", "kurtosis", "skewness")


def _mvtest_norm(s, rest: str) -> None:
    from .nptests import _royston_adjust, _sk_z
    p = parse_standard(rest)
    o = match_options(p.options, {"stats": 2, "univariate": 3, "bivariate": 3}) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 2)
    req = str(o["stats"]).split() if isinstance(o.get("stats"), str) else ["dhansen"]
    want = set()
    for w_ in req:
        if w_ == "all":
            want.update(_NORM_STATS)
            continue
        full = next((n_ for n_ in _NORM_STATS if len(w_) >= 2 and n_.startswith(w_)), None)
        if full is None:
            raise StataError(198, f"stats({w_}) invalid")
        want.add(full)
    X, w, wtype, mask = sample(s, p, names, ("fweight",))
    if w is not None:
        X = np.repeat(X, w.astype(int), axis=0)                        # fweights: replicate the observations
    out = s.output
    res: dict = {}
    if o.get("univariate"):
        out.write("\nTest for univariate normality\n\n", "text")
        out.write(f"{'':>46}------- joint ------\n", "text")
        out.write("    Variable | Pr(Skewness)   Pr(Kurtosis)  adj chi2(2)    Prob>chi2\n", "text")
        out.write(f"{'-' * 13}+{'-' * 58}\n", "text")
        U = []
        for j, nm in enumerate(names):
            z1, z2 = _sk_z(X[:, j])
            ps, pk = 2 * float(T.st.norm.sf(abs(z1))), 2 * float(T.st.norm.sf(abs(z2)))
            K2, pj = _royston_adjust(z1 * z1 + z2 * z2, X.shape[0])
            U.append([ps, pk, K2, pj])
            out.write(f"{nm[:12]:>12} |{ps:>11.4f}{pk:>15.4f}{K2:>13.2f}{pj:>13.4f}\n", "text")
        out.write(f"{'-' * 72}\n", "text")
        res["Utest"] = matrix(np.array(U), names, ["Pr_skew", "Pr_kurt", "adj_chi2", "Prob_chi2"])   # VERIFY: names
    if o.get("bivariate"):
        out.write("\nDoornik-Hansen test for bivariate normality\n\n", "text")
        out.write(f"{'-' * 51}\n", "text")
        out.write("   Pair of variables |      chi2    df   Prob>chi2\n", "text")
        out.write(f"{'-' * 21}+{'-' * 29}\n", "text")
        B = []
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                d = T.doornik_hansen(X[:, [i, j]])
                B.append([i + 1, j + 1, d["chi2"], d["df"], d["p"]])
                left = names[i][:9] if j == i + 1 else ""
                out.write(f"{left:>9} {names[j][:9]:<10} |{d['chi2']:>10.2f}{gnum(d['df']):>6}{d['p']:>12.4f}\n", "text")
        out.write(f"{'-' * 51}\n", "text")
        res["Btest"] = matrix(np.array(B), cols=["var1", "var2", "chi2", "df", "p"])          # VERIFY: names
    out.write("\nTest for multivariate normality\n\n", "text")
    rank = float(np.linalg.matrix_rank(np.cov(X, rowvar=False)))
    if "skewness" in want or "kurtosis" in want:
        m = T.mardia(X)
    if "skewness" in want:
        out.write(f"{'Mardia mSkewness =':>25}{gnum(m['b1']):>10}   chi2({gnum(m['df_skew'])}) = {m['chi2_skew']:>8.3f}"
                  f"   Prob>chi2 = {m['p_skew']:.4f}\n", "text")
        res.update({"mskew": m["b1"], "chi2_mskew": m["chi2_skew"], "df_mskew": m["df_skew"], "p_mskew": m["p_skew"],
                    "rank_mskew": rank})
    if "kurtosis" in want:
        out.write(f"{'Mardia mKurtosis =':>25}{gnum(m['b2']):>10}   chi2(1) = {m['z_kurt'] ** 2:>8.3f}"
                  f"   Prob>chi2 = {m['p_kurt']:.4f}\n", "text")
        res.update({"mkurt": m["b2"], "chi2_mkurt": m["z_kurt"] ** 2, "z_mkurt": m["z_kurt"], "p_mkurt": m["p_kurt"],
                    "rank_mkurt": rank})
    if "hzirkler" in want:
        h = T.henze_zirkler(X)
        out.write(f"{'Henze-Zirkler =':>25}{gnum(h['hz']):>10}   chi2(1) = {h['z'] ** 2:>8.3f}"
                  f"   Prob>chi2 = {h['p']:.4f}\n", "text")
        res.update({k_: h[v_] for k_, v_ in (("hz", "hz"), ("E_hz", "E"), ("V_hz", "V"), ("z_hz", "z"), ("p_hz", "p"))
                    if v_ in h})
        res["rank_hz"] = rank
    if "dhansen" in want:
        d = T.doornik_hansen(X)
        out.write(f"{'Doornik-Hansen':>18}{'':>19}chi2({gnum(d['df'])}) = {d['chi2']:>8.3f}   Prob>chi2 = {d['p']:>7.4f}\n", "text")
        res.update({"chi2_dh": d["chi2"], "df_dh": d["df"], "p_dh": d["p"]})
        res["U_dh"] = matrix(d["U"], names, ["b1", "b2", "z(b1)", "z(b2)"])
    # return list shows them in the order of the manual's Stored results (compat 1906 for the
    # Doornik-Hansen results; VERIFY for the others); the dictionary is filled in reverse
    order = ["p_dh", "df_dh", "chi2_dh", "rank_hz", "p_hz", "z_hz", "V_hz", "E_hz", "hz", "rank_mkurt", "p_mkurt",
             "z_mkurt", "chi2_mkurt", "mkurt", "rank_mskew", "p_mskew", "df_mskew", "chi2_mskew", "mskew",
             "U_dh", "Btest", "Utest"]
    s.r = {k_: res[k_] for k_ in reversed(order) if k_ in res}

"""Computations behind the small base-statistics commands (bitest, ksmirnov, nptrend, runtest,
tetrachoric, ktau, pcorr, icc, loneway, esize, ameans, lnskew0, boxcox).

Plain functions on arrays and numbers: no Session, no variables, no printing. The formulas
are the ones in the manual's Methods and formulas sections or the cited papers; layout-level
choices live in the command modules."""

from __future__ import annotations

import math

import numpy as np

from ..lazy import lazy

st = lazy("scipy.stats")
so = lazy("scipy.optimize")


# ---------------------------------------------------------------------------
# binomial test
# ---------------------------------------------------------------------------

def binomial_test(n: int, k: int, p0: float) -> tuple[float, float, float]:
    """(Pr(X >= k), Pr(X <= k), two-sided) for X ~ Binomial(n, p0). The two-sided p-value adds
    the probabilities of all outcomes not more likely than k (VERIFY against Stata)."""
    upper = float(st.binom.sf(k - 1, n, p0))
    lower = float(st.binom.cdf(k, n, p0))
    pk = st.binom.pmf(k, n, p0)
    pm = st.binom.pmf(np.arange(n + 1), n, p0)
    two = float(pm[pm <= pk * (1 + 1e-7)].sum())
    return upper, lower, min(two, 1.0)


# ---------------------------------------------------------------------------
# Kolmogorov-Smirnov
# ---------------------------------------------------------------------------

def ks_pvalue_combined(z: float) -> float:
    """Asymptotic two-sided Kolmogorov distribution: 2 sum (-1)^(j-1) exp(-2 j^2 z^2)."""
    if z <= 0:
        return 1.0
    j = np.arange(1, 101)
    return float(min(1.0, max(0.0, 2 * np.sum((-1.0) ** (j - 1) * np.exp(-2 * j * j * z * z)))))


def ks_two_sample(x: np.ndarray, y: np.ndarray) -> dict:
    """D+ (x smaller), D- (y smaller) and D with asymptotic p-values, plus Smirnov's exact
    combined p-value when feasible."""
    n1, n2 = len(x), len(y)
    allv = np.sort(np.concatenate([x, y]))
    F1 = np.searchsorted(np.sort(x), allv, side="right") / n1
    F2 = np.searchsorted(np.sort(y), allv, side="right") / n2
    d1 = float(max(0.0, (F1 - F2).max()))     # group 1 has the smaller values
    d2 = float(max(0.0, (F2 - F1).max()))
    d = max(d1, d2)
    ne = n1 * n2 / (n1 + n2)
    p1 = math.exp(-2 * ne * d1 * d1)
    p2 = math.exp(-2 * ne * d2 * d2)
    pc = ks_pvalue_combined(math.sqrt(ne) * d)
    try:
        pexact = float(st.ks_2samp(x, y, method="exact").pvalue)
    except Exception:
        pexact = math.nan
    return {"d1": d1, "d2": d2, "d": d, "p1": p1, "p2": p2, "p": pc, "pexact": pexact,
            "n1": n1, "n2": n2}


def ks_one_sample(x: np.ndarray, cdf) -> dict:
    """One-sample test of x against the theoretical cdf (a callable)."""
    n = len(x)
    xs = np.sort(x)
    F = cdf(xs)
    i = np.arange(1, n + 1)
    dplus = float(max(0.0, (i / n - F).max()))      # empirical above theoretical
    dminus = float(max(0.0, (F - (i - 1) / n).max()))
    d = max(dplus, dminus)
    p1 = math.exp(-2 * n * dplus ** 2)
    p2 = math.exp(-2 * n * dminus ** 2)
    pc = ks_pvalue_combined(math.sqrt(n) * d)
    try:
        pexact = float(st.kstwo.sf(d, n))
    except Exception:
        pexact = math.nan
    return {"dplus": dplus, "dminus": dminus, "d": d, "p1": p1, "p2": p2, "p": pc,
            "pexact": pexact, "n": n}


# ---------------------------------------------------------------------------
# nptrend (Cuzick)
# ---------------------------------------------------------------------------

def cuzick(y: np.ndarray, g: np.ndarray, score: dict) -> dict:
    """Cuzick's (1985) trend test; y values, g group codes, score[group] -> score. The
    variance is the permutation variance of the sum of rank-weighted scores."""
    N = len(y)
    r = st.rankdata(y)
    groups = sorted(score)
    n = np.array([(g == k).sum() for k in groups], dtype=float)
    R = np.array([r[g == k].sum() for k in groups])
    l = np.array([score[k] for k in groups], dtype=float)
    T = float((l * R).sum())
    ET = (N + 1) / 2 * float((l * n).sum())
    var_r = float(((r - (N + 1) / 2) ** 2).sum() / (N - 1))
    var = var_r * (N * float((l * l * n).sum()) - float((l * n).sum()) ** 2) / N
    z = (T - ET) / math.sqrt(var) if var > 0 else math.nan
    return {"z": z, "p": float(2 * st.norm.sf(abs(z))), "T": T, "ET": ET, "var": var,
            "groups": groups, "n": n, "R": R, "l": l}


# ---------------------------------------------------------------------------
# runs test
# ---------------------------------------------------------------------------

def runs_test(above: np.ndarray, continuity: bool = False) -> dict:
    """Wald-Wolfowitz runs test on a boolean sequence (True above the threshold)."""
    n = len(above)
    n1 = int(above.sum())
    n0 = n - n1
    runs = 1 + int((above[1:] != above[:-1]).sum()) if n else 0
    mu = 2.0 * n0 * n1 / n + 1
    var = 2.0 * n0 * n1 * (2.0 * n0 * n1 - n) / (n * n * (n - 1))
    diff = runs - mu
    if continuity and diff != 0:
        diff -= math.copysign(0.5, diff)
    z = diff / math.sqrt(var) if var > 0 else math.nan
    return {"n0": n0, "n1": n1, "runs": runs, "mu": mu, "var": var, "z": z,
            "p": float(2 * st.norm.sf(abs(z)))}


# ---------------------------------------------------------------------------
# Kendall's tau
# ---------------------------------------------------------------------------

def kendall(x: np.ndarray, y: np.ndarray) -> dict:
    """S, tau-a, tau-b and the tie-corrected standard error of S with the continuity-corrected
    normal p-value."""
    n = len(x)
    dx = np.sign(x[:, None] - x[None, :])
    dy = np.sign(y[:, None] - y[None, :])
    S = float(np.triu(dx * dy, 1).sum())
    n0 = n * (n - 1) / 2
    _, tc = np.unique(x, return_counts=True)
    _, uc = np.unique(y, return_counts=True)
    t_ = tc * (tc - 1) / 2
    u_ = uc * (uc - 1) / 2
    taua = S / n0
    denom = math.sqrt((n0 - t_.sum()) * (n0 - u_.sum()))
    taub = S / denom if denom > 0 else math.nan
    t = tc.astype(float)
    u = uc.astype(float)
    v0 = n * (n - 1) * (2 * n + 5)
    vt = (t * (t - 1) * (2 * t + 5)).sum()
    vu = (u * (u - 1) * (2 * u + 5)).sum()
    v1 = (t * (t - 1)).sum() * (u * (u - 1)).sum()
    v2 = (t * (t - 1) * (t - 2)).sum() * (u * (u - 1) * (u - 2)).sum()
    var = (v0 - vt - vu) / 18 + v1 / (2 * n * (n - 1))
    if n > 2:
        var += v2 / (9 * n * (n - 1) * (n - 2))
    se = math.sqrt(var)
    z = (abs(S) - 1) / se if se > 0 else math.nan
    p = float(min(1.0, 2 * st.norm.sf(max(z, 0.0)))) if not math.isnan(z) else math.nan
    return {"S": S, "taua": taua, "taub": taub, "se": se, "p": p, "n": n}


# ---------------------------------------------------------------------------
# tetrachoric correlation
# ---------------------------------------------------------------------------

_GL_X, _GL_W = np.polynomial.legendre.leggauss(24)


def bvn_cdf(h: float, k: float, r: float) -> float:
    """P(X <= h, Y <= k) for standard normals with correlation r:
    Phi(h)Phi(k) + integral_0^r phi2(h,k;t) dt (Gauss-Legendre)."""
    t = 0.5 * r * (_GL_X + 1)
    w = 0.5 * abs(r) * _GL_W * (1 if r >= 0 else -1)
    q = 1 - t * t
    f = np.exp(-(h * h - 2 * t * h * k + k * k) / (2 * q)) / (2 * math.pi * np.sqrt(q))
    return float(st.norm.cdf(h) * st.norm.cdf(k) + (w * f).sum())


def tetrachoric(table: np.ndarray) -> dict:
    """ML tetrachoric correlation from a 2x2 table [[n00, n01], [n10, n11]] (rows: first
    variable). Thresholds and rho are estimated jointly; the standard error comes from the
    numerical observed information. Returns rho, se, and the exact-test p-value (Fisher)."""
    t = np.asarray(table, dtype=float)
    N = t.sum()
    if (t.sum(0) == 0).any() or (t.sum(1) == 0).any():
        return {"rho": math.nan, "se": math.nan, "p": math.nan, "N": N}
    a0 = t[0].sum() / N
    b0 = t[:, 0].sum() / N
    h0 = float(st.norm.ppf(a0))
    k0 = float(st.norm.ppf(b0))

    def probs(th):
        h, k, r = th
        r = max(min(r, 0.999999), -0.999999)
        p00 = bvn_cdf(h, k, r)
        p0_ = float(st.norm.cdf(h))
        p_0 = float(st.norm.cdf(k))
        return np.array([[p00, p0_ - p00], [p_0 - p00, 1 - p0_ - p_0 + p00]])

    def nll(th):
        pr = np.clip(probs(th), 1e-300, None)
        return -float((t * np.log(pr)).sum())

    # zero cells: the ML estimate is +-1; report it without a standard error
    if (t == 0).any() and ((t[0, 0] * t[1, 1] == 0) or (t[0, 1] * t[1, 0] == 0)):
        rho = 1.0 if t[0, 0] * t[1, 1] > t[0, 1] * t[1, 0] else -1.0
        return {"rho": rho, "se": math.nan, "p": float(st.fisher_exact(t)[1]), "N": N}
    r0 = math.cos(math.pi / (1 + math.sqrt(max(t[0, 0] * t[1, 1], 1e-12) /
                                           max(t[0, 1] * t[1, 0], 1e-12))))
    # table rows/cols code 0 first: P(X=0, Y=0) is the lower tail, so the sign of rho is that of
    # the odds ratio
    res = so.minimize(nll, [h0, k0, r0], method="Nelder-Mead",
                      options={"xatol": 1e-10, "fatol": 1e-12, "maxiter": 4000})
    th = res.x
    H = _num_hessian(nll, th)
    try:
        V = np.linalg.inv(H)
        se = math.sqrt(V[2, 2]) if V[2, 2] > 0 else math.nan
    except np.linalg.LinAlgError:
        se = math.nan
    return {"rho": float(th[2]), "se": se, "p": float(st.fisher_exact(t)[1]), "N": N,
            "h": float(th[0]), "k": float(th[1]), "ll": -float(res.fun)}


def _num_hessian(f, x: np.ndarray, h: float = 1e-4) -> np.ndarray:
    n = len(x)
    H = np.zeros((n, n))
    f0 = f(x)
    for i in range(n):
        for j in range(i, n):
            ei = np.zeros(n)
            ej = np.zeros(n)
            ei[i] = h
            ej[j] = h
            if i == j:
                H[i, i] = (f(x + ei) - 2 * f0 + f(x - ei)) / (h * h)
            else:
                H[i, j] = H[j, i] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej)
                                     + f(x - ei - ej)) / (4 * h * h)
    return H


# ---------------------------------------------------------------------------
# partial correlations
# ---------------------------------------------------------------------------

def partial_corr(y: np.ndarray, X: np.ndarray, w: np.ndarray | None = None) -> dict:
    """Partial and semipartial correlations of y with each column of X, from the t statistics of
    the regression of y on X (with a constant)."""
    n, k = X.shape
    Z = np.column_stack([X, np.ones(n)])
    sw = np.ones(n) if w is None else np.sqrt(w)
    Zw, yw = Z * sw[:, None], y * sw
    XtX_inv = np.linalg.inv(Zw.T @ Zw)
    b = XtX_inv @ Zw.T @ yw
    e = yw - Zw @ b
    N = n if w is None else float(w.sum())
    df = N - k - 1
    rss = float(e @ e)
    tss = float(((yw - (sw * (np.sum(sw * yw) / np.sum(sw * sw)))) ** 2).sum())
    r2 = 1 - rss / tss
    s2 = rss / df
    t = b[:k] / np.sqrt(s2 * np.diag(XtX_inv)[:k])
    pc = t / np.sqrt(t * t + df)
    sc = t * np.sqrt((1 - r2) / df)
    p = 2 * st.t.sf(np.abs(t), df)
    return {"partial": pc, "semipartial": sc, "p": p, "df": df, "r2": r2}


# ---------------------------------------------------------------------------
# intraclass correlation
# ---------------------------------------------------------------------------

def anova_two_way(Y: np.ndarray) -> dict:
    """Mean squares of the n targets x k raters table (no replication)."""
    n, k = Y.shape
    gm = Y.mean()
    ssr = k * ((Y.mean(1) - gm) ** 2).sum()
    ssc = n * ((Y.mean(0) - gm) ** 2).sum()
    sst = ((Y - gm) ** 2).sum()
    sse = sst - ssr - ssc
    return {"n": n, "k": k, "BMS": ssr / (n - 1), "JMS": ssc / (k - 1),
            "EMS": sse / ((n - 1) * (k - 1)), "WMS": (sst - ssr) / (n * (k - 1))}


def icc(Y: np.ndarray, model: str = "random", kind: str = "consistency",
        level: float = 95.0) -> dict:
    """Intraclass correlations after McGraw and Wong (1996). model: oneway, random or mixed;
    kind: consistency or absolute. Returns the individual and average ICC with F-based
    confidence limits, and the F test of ICC = 0."""
    a = anova_two_way(Y)
    n, k = a["n"], a["k"]
    BMS, JMS, EMS, WMS = a["BMS"], a["JMS"], a["EMS"], a["WMS"]
    alpha = 1 - level / 100
    if model == "oneway":
        F = BMS / WMS
        df1, df2 = n - 1, n * (k - 1)
        ind = (BMS - WMS) / (BMS + (k - 1) * WMS)
        FL = F / st.f.ppf(1 - alpha / 2, df1, df2)
        FU = F * st.f.ppf(1 - alpha / 2, df2, df1)
        lo = (FL - 1) / (FL + k - 1)
        hi = (FU - 1) / (FU + k - 1)
    elif kind == "consistency":
        F = BMS / EMS
        df1, df2 = n - 1, (n - 1) * (k - 1)
        ind = (BMS - EMS) / (BMS + (k - 1) * EMS)
        FL = F / st.f.ppf(1 - alpha / 2, df1, df2)
        FU = F * st.f.ppf(1 - alpha / 2, df2, df1)
        lo = (FL - 1) / (FL + k - 1)
        hi = (FU - 1) / (FU + k - 1)
    else:
        F = BMS / EMS
        df1, df2 = n - 1, (n - 1) * (k - 1)
        ind = (BMS - EMS) / (BMS + (k - 1) * EMS + k * (JMS - EMS) / n)
        r = max(ind, 0.0)
        aa = k * r / (n * (1 - r))
        bb = 1 + k * r * (n - 1) / (n * (1 - r))
        v = (aa * JMS + bb * EMS) ** 2 / ((aa * JMS) ** 2 / (k - 1) + (bb * EMS) ** 2 / df2)
        Fs = st.f.ppf(1 - alpha / 2, n - 1, v)
        Fss = st.f.ppf(1 - alpha / 2, v, n - 1)
        lo = n * (BMS - Fs * EMS) / (Fs * (k * JMS + (k * n - k - n) * EMS) + n * BMS)
        hi = n * (Fss * BMS - EMS) / (k * JMS + (k * n - k - n) * EMS + n * Fss * BMS)
    avg = k * ind / (1 + (k - 1) * ind)
    sb = lambda r: k * r / (1 + (k - 1) * r)   # noqa: E731  Spearman-Brown
    return {"ind": ind, "ind_lb": lo, "ind_ub": hi, "avg": avg, "avg_lb": sb(lo), "avg_ub": sb(hi),
            "F": F, "df1": df1, "df2": df2, "p": float(st.f.sf(F, df1, df2)), "n": n, "k": k,
            **a}


def loneway(y: np.ndarray, g: np.ndarray, level: float = 95.0) -> dict:
    """One-way ANOVA and the intraclass correlation of [R] loneway."""
    N = len(y)
    labels, inv, ng = np.unique(g, return_inverse=True, return_counts=True)
    G = len(labels)
    gm = y.mean()
    means = np.bincount(inv, y) / ng
    ssb = float((ng * (means - gm) ** 2).sum())
    sst = float(((y - gm) ** 2).sum())
    ssw = sst - ssb
    dfb, dfw = G - 1, N - G
    msb, msw = ssb / dfb, ssw / dfw
    F = msb / msw
    n0 = (N - float((ng ** 2).sum()) / N) / (G - 1)
    rho = (F - 1) / (F + n0 - 1)
    se = math.sqrt(2 * (1 - rho) ** 2 * (1 + (n0 - 1) * rho) ** 2 / (n0 * (n0 - 1) * (G - 1))) \
        if n0 > 1 else math.nan
    alpha = 1 - level / 100
    FL = F / st.f.ppf(1 - alpha / 2, dfb, dfw)
    FU = F / st.f.ppf(alpha / 2, dfb, dfw)
    lb = (FL - 1) / (FL + n0 - 1)
    ub = (FU - 1) / (FU + n0 - 1)
    sd_b = math.sqrt(max((msb - msw) / n0, 0.0))
    sd_w = math.sqrt(msw)
    return {"N": N, "G": G, "ssb": ssb, "ssw": ssw, "sst": sst, "dfb": dfb, "dfw": dfw,
            "msb": msb, "msw": msw, "F": F, "p": float(st.f.sf(F, dfb, dfw)), "r2": ssb / sst,
            "n0": n0, "rho": rho, "se": se, "lb": lb, "ub": ub, "sd_b": sd_b, "sd_w": sd_w,
            "sd": math.sqrt(sd_b ** 2 + sd_w ** 2), "rel": (F - 1) / F, "nbar": N / G}


# ---------------------------------------------------------------------------
# effect sizes
# ---------------------------------------------------------------------------

def _nct_ci(t: float, df: float, level: float) -> tuple[float, float]:
    """Noncentrality parameters whose noncentral t puts t at the two tail probabilities."""
    a = (1 - level / 100) / 2

    def g(delta, target):
        return float(st.nct.cdf(t, df, delta)) - target

    span = max(10.0, abs(t) + 10)
    lo = so.brentq(g, -span - 20, span + 20, args=(1 - a,))
    hi = so.brentq(g, -span - 20, span + 20, args=(a,))
    return lo, hi


def effect_sizes(n1: int, m1: float, s1: float, n2: int, m2: float, s2: float,
                 level: float = 95.0, unequal: bool = False) -> dict:
    """Cohen's d, Hedges's g, Glass's Delta (both controls) and the point-biserial correlation
    of a two-group mean comparison, with noncentral-t confidence intervals."""
    df = n1 + n2 - 2
    sp = math.sqrt(((n1 - 1) * s1 * s1 + (n2 - 1) * s2 * s2) / df)
    diff = m1 - m2
    d = diff / sp
    J = math.exp(math.lgamma(df / 2) - math.lgamma((df - 1) / 2)) / math.sqrt(df / 2)
    cn = math.sqrt(n1 * n2 / (n1 + n2))
    t = d * cn
    lo, hi = _nct_ci(t, df, level)
    out = {"cohen": (d, lo / cn, hi / cn), "hedges": (J * d, J * lo / cn, J * hi / cn)}
    z = float(st.norm.ppf(1 - (1 - level / 100) / 2))
    for key, sd, nc in (("glass1", s1, n1), ("glass2", s2, n2)):
        D = diff / sd
        se = math.sqrt((n1 + n2) / (n1 * n2) + D * D / (2 * (nc - 1)))
        out[key] = (D, D - z * se, D + z * se)
    rr = lambda x: x / math.sqrt(x * x + df)   # noqa: E731
    out["pbcorr"] = (rr(t), rr(lo), rr(hi))
    out["df"] = df
    out["t"] = t
    return out


# ---------------------------------------------------------------------------
# means
# ---------------------------------------------------------------------------

def three_means(x: np.ndarray, level: float = 95.0) -> dict:
    """Arithmetic, geometric and harmonic means with t-based intervals computed on the
    original, log and reciprocal scales."""
    out = {}
    a = (1 - level / 100) / 2

    def ci(v):
        n = len(v)
        m = v.mean()
        se = v.std(ddof=1) / math.sqrt(n) if n > 1 else math.nan
        q = st.t.ppf(1 - a, n - 1) if n > 1 else math.nan
        return n, m, se, m - q * se, m + q * se

    n, m, se, lo, hi = ci(x)
    out["a"] = (n, m, se, lo, hi)
    pos = x[x > 0]
    if len(pos):
        n, m, se, lo, hi = ci(np.log(pos))
        out["g"] = (n, math.exp(m), math.exp(m) * se, math.exp(lo), math.exp(hi))
        n, m, se, lo, hi = ci(1 / pos)
        out["h"] = (n, 1 / m, se / m ** 2, 1 / hi if hi > 0 else math.nan, 1 / lo if lo > 0 else math.nan)
    return out


# ---------------------------------------------------------------------------
# skewness-zeroing shift (lnskew0)
# ---------------------------------------------------------------------------

def skewness(x: np.ndarray) -> float:
    d = x - x.mean()
    m2 = (d ** 2).mean()
    return float((d ** 3).mean() / m2 ** 1.5) if m2 > 0 else 0.0


def skew_z(x: np.ndarray) -> float:
    """D'Agostino's z for skewness (as in sktest)."""
    n = len(x)
    b1 = skewness(x)
    Y = b1 * math.sqrt((n + 1) * (n + 3) / (6 * (n - 2)))
    beta2 = 3 * (n * n + 27 * n - 70) * (n + 1) * (n + 3) / ((n - 2) * (n + 5) * (n + 7) * (n + 9))
    W2 = -1 + math.sqrt(2 * (beta2 - 1))
    delta = 1 / math.sqrt(math.log(math.sqrt(W2)))
    alpha = math.sqrt(2 / (W2 - 1))
    return delta * math.log(Y / alpha + math.sqrt((Y / alpha) ** 2 + 1))


def lnskew(x: np.ndarray, zero: float | None = None, level: float = 95.0):
    """k with skewness(ln(x - k)) = 0 (ln(k - x) for negatively skewed x) and the limits of k
    where the skewness z-test just fails to reject at `level`. Returns (k, lo, hi, sign) with
    sign = +1 for ln(x - k), -1 for ln(k - x)."""
    sign = 1 if skewness(x) >= 0 else -1
    zc = float(st.norm.ppf(1 - (1 - level / 100) / 2))
    xmin, xmax = float(x.min()), float(x.max())
    span = xmax - xmin or 1.0
    # parametrize k = xmin - span*e^u (sign 1) or xmax + span*e^u (sign -1)
    def kk(u):
        return xmin - span * math.exp(u) if sign == 1 else xmax + span * math.exp(u)

    def tr(k):
        return np.log(x - k) if sign == 1 else np.log(k - x)

    def fz(u, c):
        return skew_z(tr(kk(u))) - c

    def fs(u):
        return skewness(tr(kk(u)))

    lo_u, hi_u = -30.0, 30.0
    try:
        u0 = so.brentq(fs, lo_u, hi_u)
    except ValueError:
        return math.nan, math.nan, math.nan, sign
    ks = [kk(u0)]
    for c in (zc, -zc):
        try:
            ks.append(kk(so.brentq(fz, lo_u, hi_u, args=(c,))))
        except ValueError:
            ks.append(math.nan)
    k0, a, b = ks
    lo, hi = sorted([a, b], key=lambda v: (math.isnan(v), v))
    return k0, lo, hi, sign


# ---------------------------------------------------------------------------
# Box-Cox (lhsonly)
# ---------------------------------------------------------------------------

def _bc(y: np.ndarray, th: float) -> np.ndarray:
    return np.log(y) if abs(th) < 1e-9 else (y ** th - 1) / th


def boxcox_ll(th: float, y: np.ndarray, X: np.ndarray) -> float:
    """Log likelihood of (y^(theta) - 1)/theta = X b + e concentrated over b and sigma."""
    n = len(y)
    z = _bc(y, th)
    b, *_ = np.linalg.lstsq(X, z, rcond=None)
    rss = float(((z - X @ b) ** 2).sum())
    return -n / 2 * (math.log(2 * math.pi) + 1) - n / 2 * math.log(rss / n) \
        + (th - 1) * float(np.log(y).sum())


def boxcox_fit(y: np.ndarray, X: np.ndarray) -> dict:
    """ML estimate of theta (lhsonly model) with the standard error from the profile curvature,
    which equals the theta block of the full inverse information."""
    res = so.minimize_scalar(lambda t: -boxcox_ll(t, y, X), bounds=(-10, 10), method="bounded",
                             options={"xatol": 1e-10})
    th = float(res.x)
    h = 1e-3
    f0, fp, fm = (boxcox_ll(th, y, X), boxcox_ll(th + h, y, X), boxcox_ll(th - h, y, X))
    d2 = (fp - 2 * f0 + fm) / (h * h)
    se = math.sqrt(-1 / d2) if d2 < 0 else math.nan
    z = _bc(y, th)
    b, *_ = np.linalg.lstsq(X, z, rcond=None)
    rss = float(((z - X @ b) ** 2).sum())
    return {"theta": th, "se": se, "ll": f0, "b": b, "sigma": math.sqrt(rss / len(y)),
            "ll_at": lambda t: boxcox_ll(t, y, X), "z": z}


# ---------------------------------------------------------------------------
# Tukey smoothers ([R] smooth)
# ---------------------------------------------------------------------------

def _med_odd(x: np.ndarray, k: int) -> np.ndarray:
    """Running median of odd span k; near the ends the window shrinks symmetrically and the
    end points are copied."""
    n = len(x)
    h = k // 2
    y = x.copy()
    for i in range(1, n - 1):
        w = min(i, n - 1 - i, h)
        y[i] = np.median(x[i - w:i + w + 1])
    return y


def _endpoint(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Tukey's end-point rule on both ends: y0 = median(x0, y1, 3 y1 - 2 y2)."""
    z = y.copy()
    n = len(y)
    if n >= 3:
        z[0] = np.median([x[0], y[1], 3 * y[1] - 2 * y[2]])
        z[-1] = np.median([x[-1], y[-2], 3 * y[-2] - 2 * y[-3]])
    return z


def _split2(y: np.ndarray) -> np.ndarray:
    """Splitting of the 2-flats that are local peaks or valleys."""
    n = len(y)
    z = y.copy()
    for i in range(2, n - 3):
        if y[i] == y[i + 1] and ((y[i - 1] < y[i] > y[i + 2]) or (y[i - 1] > y[i] < y[i + 2])):
            z[i] = np.median([y[i], y[i - 1], 3 * y[i - 1] - 2 * y[i - 2]])
            z[i + 1] = np.median([y[i + 1], y[i + 2], 3 * y[i + 2] - 2 * y[i + 3]])
    return z


def tukey_smooth(x: np.ndarray, spec: str, twice: bool = False) -> np.ndarray:
    """Applies the smoother named by spec (digits: running medians of that span, an even span
    is followed by the 2-centering; R: repeat the previous median until it stops changing; S:
    split the 2-flats and resmooth with 3R; E: end-point rule; H: hanning). VERIFY the details
    against Stata's smooth."""
    spec = spec.upper()
    y = np.asarray(x, dtype=float).copy()
    last = None
    shifted = None
    i = 0
    while i < len(spec):
        c = spec[i]
        if c == "2" and shifted is not None:
            z = y.copy()
            for j in range(1, len(y) - 1):
                z[j] = (shifted[j - 1] + shifted[j]) / 2
            y = z
            shifted = None
        elif c.isdigit():
            k = int(c)
            last = k
            if k % 2:
                y = _med_odd(y, k) if k > 1 else y
            else:
                n = len(y)
                h = np.array([np.median(y[max(0, j - k // 2 + 1):min(n, j + k // 2 + 1)])
                              for j in range(n - 1)])
                shifted = h
            if k % 2 == 0 and not (i + 1 < len(spec) and spec[i + 1] == "2"):
                raise ValueError("an even span must be followed by 2")
        elif c == "R":
            if last is None:
                raise ValueError("R must follow a running median")
            for _ in range(100):
                z = _med_odd(y, last) if last % 2 else y
                if np.array_equal(z, y):
                    break
                y = z
        elif c == "S":
            z = _split2(y)
            if not np.array_equal(z, y):
                for _ in range(100):
                    z2 = _med_odd(z, 3)
                    if np.array_equal(z2, z):
                        break
                    z = z2
            y = z
        elif c == "E":
            y = _endpoint(y, np.asarray(x, dtype=float))
        elif c == "H":
            z = y.copy()
            z[1:-1] = 0.25 * y[:-2] + 0.5 * y[1:-1] + 0.25 * y[2:]
            y = z
        else:
            raise ValueError(f"unknown smoother {c}")
        i += 1
    if twice:
        r = np.asarray(x, dtype=float) - y
        y = y + tukey_smooth(r, spec, False)
    return y

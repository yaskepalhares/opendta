"""Measures of association for 2x2 tables and their stratified (Mantel-Haenszel)
versions, as described in [ST] epitab, Methods and formulas, and in the
epidemiology literature (Rothman, Greenland and Lash, Modern Epidemiology;
Breslow and Day, Statistical Methods in Cancer Research).

Notation for a cohort or case-control table:

                 Exposed   Unexposed
    Cases           a          b
    Noncases        c          d

For incidence-rate data, a and b are cases and t1, t0 the person-time of the
exposed and the unexposed. Every function returns plain floats; a missing
result is float("nan").
"""

from __future__ import annotations

import math

from ..lazy import lazy

st = lazy("scipy.stats")
opt = lazy("scipy.optimize")

NAN = float("nan")


def z_of(level: float) -> float:
    return float(st.norm.ppf(1 - (1 - level / 100) / 2))


def _div(x: float, y: float) -> float:
    return x / y if y else NAN


# ---------------------------------------------------------------------------
# chi-squared and exact tests for a 2x2 table
# ---------------------------------------------------------------------------

def pearson_2x2(a: float, b: float, c: float, d: float) -> tuple[float, float]:
    """Pearson chi2 (no continuity correction) and its p-value."""
    n = a + b + c + d
    den = (a + b) * (c + d) * (a + c) * (b + d)
    if not den:
        return NAN, NAN
    chi2 = n * (a * d - b * c) ** 2 / den
    return chi2, float(st.chi2.sf(chi2, 1))


def fisher_2x2(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    """Fisher's exact test: one-sided Pr(A >= a or A <= a, the tail of the observed
    direction) and the two-sided p (sum of the probabilities no larger than the observed)."""
    n1, m1, n = a + b, a + c, a + b + c + d
    dist = st.hypergeom(n, m1, n1)
    lo, hi = max(0, n1 - (n - m1)), min(n1, m1)
    probs = {k: float(dist.pmf(k)) for k in range(lo, hi + 1)}
    p_obs = probs[a]
    two = sum(p for p in probs.values() if p <= p_obs * (1 + 1e-7))
    expected = n1 * m1 / n if n else 0
    one = sum(p for k, p in probs.items() if k >= a) if a >= expected else \
        sum(p for k, p in probs.items() if k <= a)
    return min(one, 1.0), min(two, 1.0)


# ---------------------------------------------------------------------------
# odds ratio
# ---------------------------------------------------------------------------

def odds_ratio(a, b, c, d) -> float:
    return _div(a * d, b * c)


def or_ci_woolf(a, b, c, d, level) -> tuple[float, float]:
    if min(a, b, c, d) <= 0:
        return NAN, NAN
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lo = math.log(a * d / (b * c))
    z = z_of(level)
    return math.exp(lo - z * se), math.exp(lo + z * se)


def ci_test_based(est: float, chi2: float, level: float) -> tuple[float, float]:
    """Test-based (Miettinen) limits: est^(1 -/+ z/chi)."""
    if not (est > 0) or not (chi2 > 0):
        return NAN, NAN
    k = z_of(level) / math.sqrt(chi2)
    lo, hi = est ** (1 - k), est ** (1 + k)
    return (lo, hi) if lo <= hi else (hi, lo)


def or_ci_cornfield(a, b, c, d, level) -> tuple[float, float]:
    """Cornfield's approximate limits: the values A of the a cell (margins fixed) where
    (a - A)^2 / Var(A) = z^2, Var(A) = 1/(1/A + 1/B + 1/C + 1/D) (no continuity correction;
    VERIFY whether Stata applies one)."""
    n1, m1, n = a + b, a + c, a + b + c + d
    m0, n0 = n - m1, n - n1
    z = z_of(level)
    lo_a, hi_a = max(0.0, n1 - m0), min(float(n1), float(m1))

    def cells(A):
        return A, n1 - A, m1 - A, n0 - m1 + A

    def f(A, sign):
        A_, B_, C_, D_ = cells(A)
        v = 1 / (1 / A_ + 1 / B_ + 1 / C_ + 1 / D_)
        return (a - A) - sign * z * math.sqrt(v)

    def solve(sign, lo, hi):
        eps = 1e-10 * max(1.0, hi - lo)
        try:
            A = opt.brentq(lambda t: f(t, sign), lo + eps, hi - eps)
        except ValueError:
            return NAN
        A_, B_, C_, D_ = cells(A)
        return A_ * D_ / (B_ * C_)

    lower = solve(1.0, lo_a, float(a)) if a > lo_a else 0.0
    upper = solve(-1.0, float(a), hi_a) if a < hi_a else math.inf
    return lower, upper


def or_ci_exact(a: int, b: int, c: int, d: int, level: float) -> tuple[float, float]:
    """Exact limits from the noncentral hypergeometric distribution of a given the margins:
    Pr(A >= a; lower) = alpha/2 and Pr(A <= a; upper) = alpha/2."""
    n1, m1, n = a + b, a + c, a + b + c + d
    alpha = (1 - level / 100) / 2
    lo_a, hi_a = max(0, n1 - (n - m1)), min(n1, m1)

    def tail_ge(logpsi):
        return float(st.nchypergeom_fisher(n, m1, n1, math.exp(logpsi)).sf(a - 1))

    def tail_le(logpsi):
        return float(st.nchypergeom_fisher(n, m1, n1, math.exp(logpsi)).cdf(a))

    def root(fn):
        lo, hi = -50.0, 50.0
        return math.exp(opt.brentq(lambda t: fn(t) - alpha, lo, hi, xtol=1e-12))

    lower = 0.0 if a <= lo_a else root(tail_ge)
    upper = math.inf if a >= hi_a else root(tail_le)
    return lower, upper


def or_ci(a, b, c, d, level, method: str, chi2: float = NAN):
    if method == "woolf":
        return or_ci_woolf(a, b, c, d, level)
    if method == "tb":
        return ci_test_based(odds_ratio(a, b, c, d), chi2, level)
    if method == "cornfield":
        return or_ci_cornfield(a, b, c, d, level)
    return or_ci_exact(int(a), int(b), int(c), int(d), level)


# ---------------------------------------------------------------------------
# risk ratio and risk difference (cohort data)
# ---------------------------------------------------------------------------

def risk_ratio(a, b, c, d) -> float:
    return _div(_div(a, a + c), _div(b, b + d))


def rr_ci(a, b, c, d, level, tb=False, chi2=NAN) -> tuple[float, float]:
    rr = risk_ratio(a, b, c, d)
    if tb:
        return ci_test_based(rr, chi2, level)
    if a <= 0 or b <= 0:
        return NAN, NAN
    se = math.sqrt(1 / a - 1 / (a + c) + 1 / b - 1 / (b + d))
    z = z_of(level)
    return math.exp(math.log(rr) - z * se), math.exp(math.log(rr) + z * se)


def risk_difference(a, b, c, d) -> float:
    return _div(a, a + c) - _div(b, b + d)


def rd_ci(a, b, c, d, level, tb=False, chi2=NAN) -> tuple[float, float]:
    r1, r0 = _div(a, a + c), _div(b, b + d)
    rd = r1 - r0
    z = z_of(level)
    if tb:
        if not chi2 > 0:
            return NAN, NAN
        se = abs(rd) / math.sqrt(chi2)
    else:
        se = math.sqrt(r1 * (1 - r1) / (a + c) + r0 * (1 - r0) / (b + d))
    return rd - z * se, rd + z * se


def af_from_ratio(ratio: float, lo: float, hi: float) -> tuple[float, float, float]:
    """Attributable fraction among the exposed, (ratio - 1)/ratio, with limits from the
    ratio's limits."""
    def f(r):
        if r is None or not r > 0 or math.isinf(r):
            return 1.0 if r is not None and math.isinf(r) else NAN
        return (r - 1) / r
    return f(ratio), f(lo), f(hi)


# ---------------------------------------------------------------------------
# incidence rates
# ---------------------------------------------------------------------------

def irr_exact_ci(a: int, b: int, t1: float, t0: float, level: float) -> tuple[float, float]:
    """Limits of the rate ratio from the exact binomial limits of a/(a+b) given a+b."""
    n = a + b
    alpha = (1 - level / 100) / 2
    pl = 0.0 if a == 0 else float(st.beta.ppf(alpha, a, n - a + 1))
    pu = 1.0 if a == n else float(st.beta.ppf(1 - alpha, a + 1, n - a))
    k = t0 / t1

    def ratio(p):
        return math.inf if p >= 1 else p / (1 - p) * k
    return ratio(pl), ratio(pu)


def irr_midp(a: int, b: int, t1: float, t0: float) -> tuple[float, float]:
    """One-sided mid-p (in the direction of the observed count) and twice it."""
    n = a + b
    p = t1 / (t1 + t0)
    dist = st.binom(n, p)
    if a >= n * p:
        one = float(dist.sf(a)) + 0.5 * float(dist.pmf(a))
    else:
        one = float(dist.cdf(a - 1)) + 0.5 * float(dist.pmf(a)) if a > 0 else 0.5 * float(dist.pmf(0))
    return one, min(1.0, 2 * one)


def irr_exact_p(a: int, b: int, t1: float, t0: float) -> tuple[float, float]:
    n = a + b
    p = t1 / (t1 + t0)
    dist = st.binom(n, p)
    one = float(dist.sf(a - 1)) if a >= n * p else float(dist.cdf(a))
    return one, min(1.0, 2 * one)


# ---------------------------------------------------------------------------
# Mantel-Haenszel pooling
# ---------------------------------------------------------------------------

def mh_or(tables, level) -> tuple[float, float, float]:
    """MH odds ratio with the Robins-Breslow-Greenland limits. tables: [(a, b, c, d)]."""
    R = S = PR = PS_QR = QS = 0.0
    for a, b, c, d in tables:
        n = a + b + c + d
        if not n:
            continue
        r, s_ = a * d / n, b * c / n
        P, Q = (a + d) / n, (b + c) / n
        R += r
        S += s_
        PR += P * r
        PS_QR += P * s_ + Q * r
        QS += Q * s_
    if not R or not S:
        return _div(R, S), NAN, NAN
    est = R / S
    var = PR / (2 * R * R) + PS_QR / (2 * R * S) + QS / (2 * S * S)
    z = z_of(level)
    se = math.sqrt(var)
    return est, math.exp(math.log(est) - z * se), math.exp(math.log(est) + z * se)


def mh_chi2(tables) -> tuple[float, float]:
    """MH test that the common odds ratio is 1 (hypergeometric variance, no correction)."""
    O = E = V = 0.0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n < 2:
            continue
        n1, m1 = a + b, a + c
        O += a
        E += n1 * m1 / n
        V += n1 * (n - n1) * m1 * (n - m1) / (n * n * (n - 1))
    if not V:
        return NAN, NAN
    chi2 = (O - E) ** 2 / V
    return chi2, float(st.chi2.sf(chi2, 1))


def _expected_a(a, b, c, d, psi: float):
    """Expected a cell A (margins fixed) under the common odds ratio psi and its variance
    1/(1/A + 1/B + 1/C + 1/D); None when the table has no free cell or a zero fitted cell."""
    n1, m1 = a + b, a + c
    n = a + b + c + d
    n0, m0 = n - n1, n - m1
    lo, hi = max(0.0, n1 - m0), min(float(n1), float(m1))
    if hi <= lo:
        return None
    if psi == 1:
        A = n1 * m1 / n
    else:
        # (A (n0 - m1 + A)) / ((n1 - A)(m1 - A)) = psi: quadratic in A
        qa = 1 - psi
        qb = (n0 - m1) + psi * (n1 + m1)
        qc = -psi * n1 * m1
        disc = math.sqrt(qb * qb - 4 * qa * qc)
        roots = [(-qb + disc) / (2 * qa), (-qb - disc) / (2 * qa)]
        A = next((r for r in roots if lo - 1e-9 <= r <= hi + 1e-9), roots[0])
    B, C = n1 - A, m1 - A
    D = n0 - C
    if min(A, B, C, D) <= 0:
        return None
    return A, 1 / (1 / A + 1 / B + 1 / C + 1 / D)


def breslow_day(tables, psi: float) -> tuple[float, float, int]:
    """Breslow-Day test of homogeneity of the odds ratios around psi (without Tarone's
    correction; VERIFY which homogeneity test Stata's cc reports)."""
    chi2 = 0.0
    k = 0
    for t in tables:
        r = _expected_a(*t, psi)
        if r is not None:
            chi2 += (t[0] - r[0]) ** 2 / r[1]
            k += 1
    df = k - 1
    return chi2, float(st.chi2.sf(chi2, df)) if df > 0 else NAN, df


def mh_rr(tables, level) -> tuple[float, float, float]:
    """MH risk ratio with the Greenland-Robins variance of its log."""
    num = den = vnum = 0.0
    for a, b, c, d in tables:
        n1, n0 = a + c, b + d
        n = n1 + n0
        if not n:
            continue
        num += a * n0 / n
        den += b * n1 / n
        vnum += ((a + b) * n1 * n0 - a * b * n) / (n * n)
    if not num or not den:
        return _div(num, den), NAN, NAN
    est = num / den
    se = math.sqrt(vnum / (num * den))
    z = z_of(level)
    return est, math.exp(math.log(est) - z * se), math.exp(math.log(est) + z * se)


def mh_irr(strata, level) -> tuple[float, float, float]:
    """MH rate ratio, strata = [(a, b, t1, t0)], Greenland-Robins variance."""
    num = den = vnum = 0.0
    for a, b, t1, t0 in strata:
        t = t1 + t0
        if not t:
            continue
        num += a * t0 / t
        den += b * t1 / t
        vnum += (a + b) * t1 * t0 / (t * t)
    if not num or not den:
        return _div(num, den), NAN, NAN
    est = num / den
    se = math.sqrt(vnum / (num * den))
    z = z_of(level)
    return est, math.exp(math.log(est) - z * se), math.exp(math.log(est) + z * se)


def homogeneity_log(estimates, variances, pooled) -> tuple[float, float, int]:
    """Chi2 = sum (log est_i - log pooled)^2 / var(log est_i) over the strata with a finite
    estimate (VERIFY: Stata's homogeneity test for cs and ir)."""
    chi2 = 0.0
    k = 0
    for e, v in zip(estimates, variances):
        if e > 0 and v > 0 and not math.isinf(e) and pooled > 0:
            chi2 += (math.log(e) - math.log(pooled)) ** 2 / v
            k += 1
    df = k - 1
    return chi2, float(st.chi2.sf(chi2, df)) if df > 0 else NAN, df


def tarone(tables, psi: float) -> tuple[float, float, int]:
    """Tarone's correction of the Breslow-Day statistic around psi:
    BD - (sum a_i - sum A_i)^2 / sum Var(A_i), with A_i, Var(A_i) the expected a cell and its
    variance under the common odds ratio psi."""
    bd, _, df = breslow_day(tables, psi)
    obs = exp = var = 0.0
    for t in tables:
        r = _expected_a(*t, psi)
        if r is not None:
            obs += t[0]
            exp += r[0]
            var += r[1]
    chi2 = bd - (obs - exp) ** 2 / var if var else NAN
    return chi2, float(st.chi2.sf(chi2, df)) if df > 0 else NAN, df


# ---------------------------------------------------------------------------
# standardized and directly pooled estimates ([ST] epitab: estandard, istandard,
# standard(), pool)
# ---------------------------------------------------------------------------

def standardized(m1, v1, m0, v0, w, level, difference: bool = False) -> tuple[float, float, float]:
    """Standardized ratio sum w m1 / sum w m0 (or difference sum w (m1 - m0) / sum w) of the
    stratum measures m1 (exposed) and m0 (unexposed) with variances v1, v0, and its limits:
    var(ln ratio) = sum w^2 v1 / (sum w m1)^2 + sum w^2 v0 / (sum w m0)^2 and
    var(difference) = sum w^2 (v1 + v0) / (sum w)^2 (Rothman and Greenland; VERIFY against
    Stata's limits)."""
    z = z_of(level)
    W = sum(w)
    S1 = sum(wi * x for wi, x in zip(w, m1))
    S0 = sum(wi * x for wi, x in zip(w, m0))
    V1 = sum(wi * wi * x for wi, x in zip(w, v1))
    V0 = sum(wi * wi * x for wi, x in zip(w, v0))
    if difference:
        if not W:
            return NAN, NAN, NAN
        est = (S1 - S0) / W
        se = math.sqrt(V1 + V0) / W
        return est, est - z * se, est + z * se
    if not S1 or not S0:
        return _div(S1, S0), NAN, NAN
    est = S1 / S0
    se = math.sqrt(V1 / (S1 * S1) + V0 / (S0 * S0))
    return est, est * math.exp(-z * se), est * math.exp(z * se)


def pooled(estimates, variances, level, log: bool = True):
    """Directly pooled estimate: the inverse-variance weighted mean of the stratum estimates
    (of their logs for ratios), with the test of homogeneity
    sum (est_i - pooled)^2 / var_i. Strata without a finite estimate or variance are left out
    (VERIFY). Returns (est, lb, ub, chi2, p, df)."""
    num = den = 0.0
    keep = []
    for e, v in zip(estimates, variances):
        if log and not (e > 0 and not math.isinf(e)):
            continue
        if not (v > 0) or math.isinf(v) or math.isnan(e):
            continue
        x = math.log(e) if log else e
        keep.append((x, v))
        num += x / v
        den += 1 / v
    if not den:
        return NAN, NAN, NAN, NAN, NAN, 0
    m = num / den
    se = math.sqrt(1 / den)
    z = z_of(level)
    chi2 = sum((x - m) ** 2 / v for x, v in keep)
    df = len(keep) - 1
    p = float(st.chi2.sf(chi2, df)) if df > 0 else NAN
    if log:
        return math.exp(m), math.exp(m - z * se), math.exp(m + z * se), chi2, p, df
    return m, m - z * se, m + z * se, chi2, p, df

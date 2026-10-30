"""Formulas for power and sample-size analysis, part 2: twoway, repeated, cmh, mcc, trend, cox,
exponential and logrank ([PSS] power twoway ... power logrank, Methods and formulas).

Plain functions: they take numbers and arrays and return numbers, with no Session, no variables and
no printing. A failed search raises ValueError; commands/power3.py turns it into Stata's error.
Rules the manual does not state are marked VERIFY (case compat/do/1919_other_suites2.do).
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np

from ..lazy import lazy

st = lazy("scipy.stats")
opt = lazy("scipy.optimize")


def zq(p: float) -> float:
    return float(st.norm.ppf(p))


# ---------------------------------------------------------------------------
# searches
# ---------------------------------------------------------------------------

def smallest_int(power_of: Callable[[int], float], target: float, lo: int = 1, cap: int = 10 ** 9) -> int:
    """Smallest integer k >= lo with power_of(k) >= target (power increasing in k; nan counts as
    not reached). VERIFY: a saw-toothed power function can have a smaller k that fails."""
    k = max(int(lo), 1)
    if power_of(k) >= target:
        return k
    hi = max(2 * k, 2)
    while not power_of(hi) >= target:
        k, hi = hi, hi * 2
        if hi > cap:
            raise ValueError("no solution")
    while hi - k > 1:
        mid = (hi + k) // 2
        if power_of(mid) >= target:
            hi = mid
        else:
            k = mid
    return hi


def solve_up(f: Callable[[float], float], lo: float, hi: float) -> float:
    """Root of an increasing function f on [lo, hi]; hi is doubled until f(hi) >= 0."""
    if f(lo) >= 0:
        return lo
    for _ in range(80):
        if f(hi) >= 0:
            return float(opt.brentq(f, lo, hi, xtol=1e-13, rtol=1e-13, maxiter=500))
        lo, hi = hi, hi * 2
    raise ValueError("no solution")


def sides_z(alpha: float, sides: int) -> float:
    return zq(1 - alpha / sides)


# ---------------------------------------------------------------------------
# ANOVA: noncentral F
# ---------------------------------------------------------------------------

def anova_power(lam: float, df1: float, df2: float, alpha: float) -> float:
    if df2 <= 0:
        return math.nan
    crit = float(st.f.isf(alpha, df1, df2))
    return float(min(max(st.ncf.sf(crit, df1, df2, lam), 0.0), 1.0))


def lambda_for_power(df1: float, df2: float, alpha: float, power: float) -> float:
    """Noncentrality that gives the requested power."""
    return solve_up(lambda lam: anova_power(lam, df1, df2, alpha) - power, 1e-9, 10.0)


def standard_weights(w: list[float]) -> list[float]:
    """Smallest weights with the same ratios: divided by the gcd of integers, else by the
    minimum. VERIFY: the manual refers to [PSS] power oneway for the standardization."""
    if all(float(x).is_integer() for x in w):
        g = 0
        for x in w:
            g = math.gcd(g, int(x))
        return [float(int(x) // g) for x in w]
    m = min(w)
    return [x / m for x in w]


# ---------------------------------------------------------------------------
# twoway
# ---------------------------------------------------------------------------

def twoway_variances(means: np.ndarray) -> tuple[float, float, float]:
    """Row, column and row-by-column variances explained: sum a_j^2 / J, ..."""
    m = np.asarray(means, dtype=float)
    gm = m.mean()
    a = m.mean(axis=1) - gm
    b = m.mean(axis=0) - gm
    ab = m - m.mean(axis=1)[:, None] - m.mean(axis=0)[None, :] + gm
    return float((a ** 2).mean()), float((b ** 2).mean()), float((ab ** 2).mean())


def twoway_df1(J: int, K: int, which: str) -> int:
    return {"row": J - 1, "column": K - 1}.get(which, (J - 1) * (K - 1))


def twoway_contrast(J: int, K: int, which: str) -> np.ndarray:
    dj = np.hstack([np.eye(J - 1), -np.ones((J - 1, 1))])
    dk = np.hstack([np.eye(K - 1), -np.ones((K - 1, 1))])
    if which == "row":
        return np.kron(dj, np.ones((1, K)) / K)
    if which == "column":
        return np.kron(np.ones((1, J)) / J, dk)
    return np.kron(dj, dk)


def twoway_noncentrality(means: np.ndarray, frac: np.ndarray, n: float, which: str, var_e: float) -> float:
    """lambda = n (C mu)' {C W^-1 C'}^-1 (C mu) / var_e, cells in row-major order."""
    J, K = np.asarray(means).shape
    c = twoway_contrast(J, K, which)
    cm = c @ np.asarray(means, dtype=float).reshape(-1)
    w = np.asarray(frac, dtype=float).reshape(-1)
    a = c @ np.diag(1.0 / w) @ c.T
    return float(n * cm @ np.linalg.solve(a, cm) / var_e)


def twoway_power(lam: float, J: int, K: int, n: float, which: str, alpha: float) -> float:
    return anova_power(lam, twoway_df1(J, K, which), n - J * K, alpha)


# ---------------------------------------------------------------------------
# repeated
# ---------------------------------------------------------------------------

def cs_covariance(K: int, var: float, rho: float) -> np.ndarray:
    return var * ((1 - rho) * np.eye(K) + rho * np.ones((K, K)))


def _within_basis(K: int) -> np.ndarray:
    q, _ = np.linalg.qr(np.column_stack([np.ones(K) / math.sqrt(K), np.eye(K)[:, :K - 1]]))
    return q[:, 1:]


def repeated_summary(J: int, K: int, cov: np.ndarray, which: str, weights: np.ndarray,
                     means: np.ndarray | None = None, var_effect: float | None = None) -> dict:
    """delta^2, degrees of freedom, epsilon and the variance terms of [PSS] power repeated.
    The columns of U are the orthonormal ones divided by sqrt(K), as in Stata's Var_be, Var_we and
    Var_w (checked against the Stata 14 log of case 1919)."""
    if which in ("between", "bwithin"):
        if J < 2:
            raise ValueError("at least two groups are required for a between-subjects effect")
        c = np.hstack([np.eye(J - 1), -np.ones((J - 1, 1))])
    else:
        c = np.ones((1, J)) / J
    u = np.ones((K, 1)) / math.sqrt(K) if which == "between" else _within_basis(K)
    u = u / math.sqrt(K)
    dc, du = c.shape[0], u.shape[1]
    sstar = u.T @ np.asarray(cov, dtype=float) @ u
    tr_s = float(np.trace(sstar))
    eps = tr_s ** 2 / (du * float(np.trace(sstar @ sstar)))
    var_e = tr_s / (dc * du)
    if means is not None:
        theta = c @ np.asarray(means, dtype=float) @ u
        a = c @ np.diag(1.0 / np.asarray(weights, dtype=float)) @ c.T
        tr_h = float(np.trace(theta.T @ np.linalg.solve(a, theta)))
        var_x = tr_h / dc
    else:
        var_x = float(var_effect)
    return {"dc": dc, "du": du, "eps": eps, "var_x": var_x, "var_xe": var_e, "delta2": var_x / var_e,
            "spherical": abs(eps - 1) < 1e-8}


def repeated_power(delta2: float, dc: int, du: int, p: int, n: float, eps: float, alpha: float,
                   lower: bool = False) -> float:
    """F test with the Greenhouse-Geisser correction: epsilon scales both degrees of freedom and
    the noncentrality. Stata's multiplier m of Muller and Barton (1989) is not documented; this
    form reproduces the spherical cases of the Stata 14 log exactly and a non-spherical one to
    1e-4 (compat/known_differences.txt)."""
    df1, df2 = dc * du * eps, du * (n - p) * eps
    if df2 <= 0:
        return math.nan
    crit = float(st.f.isf(alpha, df1, df2))
    tail = st.ncf.cdf if lower else st.ncf.sf       # beta from the lower tail: no cancellation
    return float(min(max(tail(crit, df1, df2, n * delta2 * eps), 0.0), 1.0))


# ---------------------------------------------------------------------------
# cmh
# ---------------------------------------------------------------------------

def odds_probability(p1: np.ndarray, oratio: float) -> np.ndarray:
    """Experimental-group probability whose odds are oratio times the odds of p1."""
    p1 = np.asarray(p1, dtype=float)
    return oratio * p1 / (1 - p1 + oratio * p1)


def cmh_power(n1k, n2k, p1k, p2k, alpha: float, sides: int, continuity: bool = False,
              upper: bool | None = None, tail: bool = False) -> float:
    """Asymptotic power of the Cochran-Mantel-Haenszel test (Nam 1992)."""
    n1k, n2k, p1k, p2k = (np.asarray(x, dtype=float) for x in (n1k, n2k, p1k, p2k))
    if (n1k <= 0).any() or (n2k <= 0).any():
        return math.nan
    nk = n1k + n2k
    om = n1k * n2k / nk
    pk = (n1k * p1k + n2k * p2k) / nk
    ew = float((om * (p2k - p1k)).sum())
    v0 = float((om * pk * (1 - pk)).sum())
    v = float((om ** 2 * (p1k * (1 - p1k) / n1k + p2k * (1 - p2k) / n2k)).sum())
    z = sides_z(alpha, sides)
    cc = 0.5 if continuity else 0.0
    uu = (z * math.sqrt(v0) - ew + cc) / math.sqrt(v)
    ul = (-z * math.sqrt(v0) - ew - cc) / math.sqrt(v)
    pu, pl = float(st.norm.sf(uu)), float(st.norm.cdf(ul))
    if sides == 2 and not tail:         # tail: only the tail in the direction of the effect
        return pu + pl
    up = (ew > 0) if upper is None else upper
    return pu if up else pl


def cmh_sizes(weights: list[float], s: list[float], k: float, fractional: bool = False):
    """Group sizes per stratum for the multiplier k: stratum size weight x k. Equal groups
    (ratio 1/2 everywhere) split it in halves, which may be fractional (N = 201 gives 33.5 per
    group and stratum in Stata 14); otherwise the experimental group is rounded up and the
    control group takes the rest (Stata 14, case 1919)."""
    w = np.asarray(weights, dtype=float)
    sv = np.asarray(s, dtype=float)
    nk = w * k
    if np.allclose(sv, 0.5):
        return nk / 2, nk / 2
    n2 = nk * sv
    if not fractional:
        n2 = np.ceil(n2 - 1e-9)
    return nk - n2, n2


# ---------------------------------------------------------------------------
# mcc
# ---------------------------------------------------------------------------

def mcc_p1(p0: float, theta: float, rho: float) -> float:
    """Case exposure probability with odds ratio theta = p10 / p01 at correlation rho."""
    if rho == 0:
        return float(odds_probability(np.array(p0), theta))
    q0 = 1 - p0

    def g(p1):
        s = rho * math.sqrt(p1 * (1 - p1) * p0 * q0)
        return (p1 * q0 - s) - theta * ((1 - p1) * p0 - s)
    eps = 1e-9
    try:
        return float(opt.brentq(g, eps, 1 - eps, xtol=1e-14))
    except ValueError:
        raise ValueError("no valid exposure probability for these values")


def mcc_terms(p0: float, p1: float, rho: float, M: int, theta: float):
    """e_theta, nu_theta and the same at theta = 1 (VERIFY: E_1 and nu_1 are taken at theta = 1)."""
    q0, q1 = 1 - p0, 1 - p1
    s = rho * math.sqrt(p1 * q1 * p0 * q0)
    p11, p01 = p1 * p0 + s, q1 * p0 - s
    pp, pm = p11 / p1, p01 / q1
    qp, qm = 1 - pp, 1 - pm
    if min(pp, pm, qp, qm) < -1e-12:
        raise ValueError("invalid probabilities")
    pp, pm, qp, qm = (min(max(x, 0.0), 1.0) for x in (pp, pm, qp, qm))
    t = [p1 * math.comb(M, m - 1) * pp ** (m - 1) * qp ** (M - m + 1)
         + q1 * math.comb(M, m) * pm ** m * qm ** (M - m) for m in range(1, M + 1)]
    e = nu = e1 = nu1 = 0.0
    for m, tm in enumerate(t, 1):
        d = m * theta + M - m + 1
        e += m * tm * theta / d
        nu += m * tm * theta * (M - m + 1) / d ** 2
        e1 += m * tm / (M + 1)
        nu1 += m * tm * (M - m + 1) / (M + 1) ** 2
    return e, nu, e1, nu1


def mcc_power(n: float, p0: float, theta: float, rho: float, M: int, alpha: float, sides: int,
              upper: bool | None = None) -> float:
    p1 = mcc_p1(p0, theta, rho)
    e, nu, e1, nu1 = mcc_terms(p0, p1, rho, M, theta)
    z = sides_z(alpha, sides)
    s1, st_ = math.sqrt(n * nu1), math.sqrt(n * nu)
    lo = (n * (e1 - e) - z * s1) / st_
    up = (n * (e1 - e) + z * s1) / st_
    pl, pu = float(st.norm.cdf(lo)), float(st.norm.sf(up))
    if sides == 2:
        return pl + pu
    return pu if (theta > 1 if upper is None else upper) else pl


def mcc_size_start(p0: float, theta: float, rho: float, M: int, alpha: float, sides: int, power: float) -> float:
    """Closed form for the tail in the direction of the effect: Stata 14 reports it as r(init),
    the starting value of its solver (it uses z of 1 - alpha/sides, case 1919)."""
    p1 = mcc_p1(p0, theta, rho)
    e, nu, e1, nu1 = mcc_terms(p0, p1, rho, M, theta)
    za, zb = sides_z(alpha, sides), zq(power)
    return (zb * math.sqrt(nu) + za * math.sqrt(nu1)) ** 2 / (e1 - e) ** 2


# ---------------------------------------------------------------------------
# trend
# ---------------------------------------------------------------------------

def trend_power(nj, pj, xj, alpha: float, sides: int, continuity: bool = False, tail: bool = False) -> float:
    """Power of the Cochran-Armitage test for logit(p) = a + b x (Nam 1987). VERIFY: the signs of
    the one-sided bounds are the statistically coherent ones; the manual's rendering is garbled."""
    nj, pj, xj = (np.asarray(v, dtype=float) for v in (nj, pj, xj))
    n = nj.sum()
    pbar = float((nj * pj).sum() / n)
    xbar = float((nj * xj).sum() / n)
    v0 = pbar * (1 - pbar) * float((nj * (xj - xbar) ** 2).sum())
    v = float((nj * pj * (1 - pj) * (xj - xbar) ** 2).sum())
    e = float((nj * pj * (xj - xbar)).sum())
    c = (xj[1] - xj[0]) / 2 if continuity else 0.0
    z = sides_z(alpha, sides)
    uu = (z * math.sqrt(v0) + c - e) / math.sqrt(v)
    ul = (-z * math.sqrt(v0) - c - e) / math.sqrt(v)
    pu, pl = float(st.norm.sf(uu)), float(st.norm.cdf(ul))
    if sides == 2 and not tail:
        return pu + pl
    return pu if e > 0 else pl


def is_monotonic(p) -> bool:
    d = np.diff(np.asarray(p, dtype=float))
    return bool((d > 0).all() or (d < 0).all())


# ---------------------------------------------------------------------------
# cox
# ---------------------------------------------------------------------------

def cox_events(b1: float, sd: float, r2: float, alpha: float, sides: int, power: float) -> float:
    """Events for the Wald test of a covariate (Hsieh and Lavori 2000)."""
    return (sides_z(alpha, sides) + zq(power)) ** 2 / (sd ** 2 * b1 ** 2 * (1 - r2))


def cox_power(b1: float, sd: float, r2: float, n: float, pe: float, alpha: float, sides: int) -> float:
    return float(st.norm.cdf(abs(b1) * sd * math.sqrt(n * pe * (1 - r2)) - sides_z(alpha, sides)))


def cox_effect(sd: float, r2: float, n: float, pe: float, alpha: float, sides: int, power: float) -> float:
    """|b1| that reaches the power."""
    return (sides_z(alpha, sides) + zq(power)) / (sd * math.sqrt(n * pe * (1 - r2)))


# ---------------------------------------------------------------------------
# exponential
# ---------------------------------------------------------------------------

def exp_pe(lam: float, eta: float, T: float | None, r: float, gamma: float) -> float:
    """Probability of an event: follow-up to the end of the study T, accrual over [0, r] with
    truncated-exponential shape gamma (uniform for gamma = 0), loss hazard eta. T = None: followed
    until failure."""
    if T is None:
        return lam / (lam + eta)
    v = lam + eta
    if r < 1e-12:
        return lam / v * (1 - math.exp(-v * T))
    if abs(gamma) < 1e-6:
        return lam / v * (1 - (math.exp(-v * (T - r)) - math.exp(-v * T)) / (v * r))
    if abs(v - gamma) < 1e-9:
        v += 1e-7
    return lam / v * (1 + gamma * math.exp(-v * T) * (1 - math.exp((v - gamma) * r))
                      / ((v - gamma) * (1 - math.exp(-gamma * r))))


def exp_gamma(r: float, ta: float, pa: float) -> float:
    """Accrual shape with a proportion pa accrued by time ta (0 when uniform)."""
    if abs(pa - ta / r) < 1e-12:
        return 0.0

    def g(gm):
        return (1 - math.exp(-gm * ta)) / (1 - math.exp(-gm * r)) - pa
    if pa > ta / r:
        return float(opt.brentq(g, 1e-9, 1000.0 / r, xtol=1e-14))
    return float(opt.brentq(g, -1000.0 / r, -1e-9, xtol=1e-14))


def exp_xi(l1: float, l2: float, eta1: float, eta2: float, R: float, T, r: float, gamma: float,
           logtest: bool, unconditional: bool = False) -> tuple[float, float]:
    """(xi_o, xi_a) of Lachin and Foulkes (1986), hazard-difference or log-hazard-ratio test."""
    p1, p2 = 1 / (1 + R), R / (1 + R)
    lam = p1 * l1 + p2 * l2
    if unconditional:
        lam = l1          # VERIFY: lambda = lambda_1 under H0 for unconditional tests

    def zeta(lm, eta):
        pe = exp_pe(lm, eta, T, r, gamma)
        return 1 / pe if logtest else lm ** 2 / pe
    return (zeta(lam, eta1) / p1 + zeta(lam, eta2) / p2, zeta(l1, eta1) / p1 + zeta(l2, eta2) / p2)


def exp_size(delta: float, xo: float, xa: float, alpha: float, sides: int, power: float,
             unconditional: bool = False) -> float:
    za, zb = sides_z(alpha, sides), zq(power)
    if unconditional:
        return (za + zb) ** 2 * xa / delta ** 2
    return (za * math.sqrt(xo) + zb * math.sqrt(xa)) ** 2 / delta ** 2


def exp_power(delta: float, xo: float, xa: float, n: float, alpha: float, sides: int,
              unconditional: bool = False) -> float:
    za = sides_z(alpha, sides)
    if unconditional:
        return float(st.norm.cdf(abs(delta) * math.sqrt(n / xa) - za))
    return float(st.norm.cdf((abs(delta) * math.sqrt(n) - za * math.sqrt(xo)) / math.sqrt(xa)))


# ---------------------------------------------------------------------------
# logrank
# ---------------------------------------------------------------------------

def logrank_psi(hr: float, R: float, schoenfeld: bool) -> float:
    return (1 + R) / math.log(hr) if schoenfeld else (R * hr + 1) / (hr - 1)


def logrank_events(hr: float, R: float, alpha: float, sides: int, power: float, schoenfeld: bool) -> float:
    return (sides_z(alpha, sides) + zq(power)) ** 2 * logrank_psi(hr, R, schoenfeld) ** 2 / R


def logrank_power(hr: float, R: float, n: float, pe: float, alpha: float, sides: int,
                  schoenfeld: bool) -> float:
    return float(st.norm.cdf(math.sqrt(R * n * pe) / abs(logrank_psi(hr, R, schoenfeld)) - sides_z(alpha, sides)))


def logrank_pe(hr: float, R: float, s1: float | None, s2: float | None = None,
               simpson: tuple[float, float, float] | None = None) -> float:
    """Probability of an event. No survival given: no censoring. Default: average follow-up
    (Freedman); simpson = control survival at f, 0.5 r + f and T."""
    if simpson is not None:
        se = [(s + R * s ** hr) / (1 + R) for s in simpson]
        return 1 - (se[0] + 4 * se[1] + se[2]) / 6
    if s1 is None:
        return 1.0
    s2 = s1 ** hr if s2 is None else s2
    return 1 - (s1 + R * s2) / (1 + R)

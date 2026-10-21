"""Random-effects panel likelihoods by Gauss-Hermite quadrature, and the closed-form
likelihoods of the negative binomial panel models ([XT] xtprobit, xtcloglog, xtologit,
xtoprobit, xttobit, xtintreg, xtnbreg, Methods and formulas).

Everything takes arrays and numbers and returns numbers or arrays: no Session, no variables,
no printing. A model enters the quadrature through `f(u)`, the log likelihood of every
observation at the linear predictor `u` (any array shape, vectorized); the random effect
v_i ~ N(0, s^2) shifts the predictor to u + v_i.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Callable

import numpy as np

from ..lazy import lazy

special = lazy("scipy.special")
st = lazy("scipy.stats")

LOG_SQRT_2PI = 0.5 * math.log(2 * math.pi)


@lru_cache(maxsize=None)
def gauss_hermite(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Abscissas and weights of the n-point Gauss-Hermite rule for exp(-x^2)."""
    x, w = np.polynomial.hermite.hermgauss(n)
    return x, w


class Groups:
    """Observations sorted by group: `order` puts them in group order, `starts` are the first
    positions of every group in that order, `inv` is the group number of each original row."""

    def __init__(self, inv: np.ndarray):
        self.inv = np.asarray(inv)
        self.G = int(self.inv.max()) + 1
        self.order = np.argsort(self.inv, kind="stable")
        self.sorted_inv = self.inv[self.order]
        self.counts = np.bincount(self.inv, minlength=self.G)
        self.starts = np.concatenate([[0], np.cumsum(self.counts)[:-1]])

    def sums(self, a: np.ndarray) -> np.ndarray:
        """Sum of `a` (rows in the original order; first axis) within each group."""
        return np.add.reduceat(a[self.order], self.starts, axis=0)


def _adapt(f: Callable, eta: np.ndarray, grp: Groups, sigma: float) -> tuple[np.ndarray, np.ndarray]:
    """Mode and scale of the integrand of every group, in the random effect v: Newton steps
    on sum_t f(eta_t + v) - v^2/(2 sigma^2), with finite differences for the derivatives."""
    mu = np.zeros(grp.G)
    d = 1e-3
    tau = np.full(grp.G, sigma)
    for _ in range(25):
        u = eta + mu[grp.inv]
        f0, fp, fm = f(u), f(u + d), f(u - d)
        d1 = grp.sums((fp - fm) / (2 * d)) - mu / sigma ** 2
        d2 = grp.sums((fp - 2 * f0 + fm) / d ** 2) - 1 / sigma ** 2
        d2 = np.minimum(d2, -1 / sigma ** 2)
        step = -d1 / d2
        mu = mu + step
        tau = 1 / np.sqrt(-d2)
        if np.max(np.abs(step)) < 1e-10:
            break
    return mu, tau


def re_loglik_groups(f: Callable, eta: np.ndarray, grp: Groups, sigma: float, nq: int = 12,
                     adaptive: bool = True) -> np.ndarray:
    """Log likelihood of every group, integrating the random effect by Gauss-Hermite with nq
    points: adaptive (centered at the mode of each group, scaled by its curvature:
    mvaghermite) or at zero scaled by sigma (ghermite)."""
    x, w = gauss_hermite(nq)
    if adaptive:
        mu, tau = _adapt(f, eta, grp, sigma)
        nodes = mu[:, None] + math.sqrt(2) * tau[:, None] * x[None, :]          # G x Q
        logj = np.log(math.sqrt(2) * tau)[:, None]
        prior = -0.5 * (nodes / sigma) ** 2 - math.log(sigma) - LOG_SQRT_2PI
        terms = np.log(w)[None, :] + (x ** 2)[None, :] + prior + logj
    else:
        nodes = np.broadcast_to(math.sqrt(2) * sigma * x, (grp.G, nq))
        terms = np.broadcast_to(np.log(w) - 0.5 * math.log(math.pi), (grp.G, nq))
    S = np.empty((grp.G, nq))
    for q in range(nq):
        S[:, q] = grp.sums(f(eta + nodes[grp.inv, q]))
    return special.logsumexp(S + terms, axis=1)


def re_loglik(f, eta, grp, sigma, nq=12, adaptive=True) -> float:
    return float(re_loglik_groups(f, eta, grp, sigma, nq, adaptive).sum())


def tau_grid(f, eta, grp, pooled_ll: float, nq: int = 12, step: float = 0.1,
             limit: float = 0.95, sigma_of=None) -> list[tuple[float, float]]:
    """The starting-value search of xtlogit-like commands: log likelihood with the pooled
    coefficients and plain quadrature for tau = 0, 0.1, ... (sigma^2 = tau/(1-tau)), until it
    falls. `sigma_of(tau)` converts tau to sigma (default tau/(1-tau))."""
    grid: list[tuple[float, float]] = []
    t = 0.0
    while t < limit:
        if t == 0:
            val = pooled_ll
        else:
            sg = math.sqrt(t / (1 - t)) if sigma_of is None else sigma_of(t)
            val = re_loglik(f, eta, grp, sg, nq, adaptive=False)
        grid.append((t, val))
        if len(grid) > 1 and grid[-1][1] < grid[-2][1]:
            break
        t = round(t + step, 1)
    return grid


# ---------------------------------------------------------------------------
# log likelihood of one observation at the linear predictor u
# ---------------------------------------------------------------------------

def ll_probit(u: np.ndarray, sign: np.ndarray) -> np.ndarray:
    """sign = 2y - 1."""
    return special.log_ndtr(sign * u)


def ll_logit(u: np.ndarray, sign: np.ndarray) -> np.ndarray:
    return special.log_expit(sign * u)


def ll_cloglog(u: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Pr(y=1) = 1 - exp(-exp(u))."""
    e = np.exp(np.clip(u, -700, 700))
    with np.errstate(divide="ignore"):
        p1 = np.log(-np.expm1(-e))
    return np.where(y != 0, p1, -e)


def ll_ordered(u: np.ndarray, yi: np.ndarray, cuts: np.ndarray, logistic: bool) -> np.ndarray:
    """Ordered outcome with categories 0..J-1 and J-1 increasing cut points:
    Pr(y=j) = F(k_j - u) - F(k_{j-1} - u)."""
    J = len(cuts) + 1
    ext = np.concatenate([[-np.inf], cuts, [np.inf]])
    up = ext[yi + 1]
    lo = ext[yi]
    if u.ndim > yi.ndim:
        up, lo = up[:, None], lo[:, None]
    F = special.expit if logistic else st.norm.cdf
    with np.errstate(invalid="ignore"):
        p = F(up - u) - F(lo - u)
    return np.log(np.maximum(p, 1e-300))


def ll_interval(u: np.ndarray, y1: np.ndarray, y2: np.ndarray, sigma_e: float) -> np.ndarray:
    """Interval data (xtintreg): y1 and y2 are the lower and upper limits; NaN lower is left
    censoring, NaN upper is right censoring, equal limits are an exact observation."""
    if u.ndim > y1.ndim:
        y1, y2 = y1[:, None], y2[:, None]
    exact = (y1 == y2)
    lo_miss, hi_miss = np.isnan(y1), np.isnan(y2)
    z1 = (np.where(lo_miss, 0.0, y1) - u) / sigma_e
    z2 = (np.where(hi_miss, 0.0, y2) - u) / sigma_e
    pt = -0.5 * z1 ** 2 - LOG_SQRT_2PI - math.log(sigma_e)
    left = special.log_ndtr(z2)
    right = special.log_ndtr(-z1)
    both = np.log(np.maximum(st.norm.cdf(z2) - st.norm.cdf(z1), 1e-300))
    return np.where(exact, pt, np.where(lo_miss, left, np.where(hi_miss, right, both)))


def ll_tobit(u: np.ndarray, y: np.ndarray, sigma_e: float, ll: float | None, ul: float | None) -> np.ndarray:
    """Tobit with the limits ll and ul (None: no limit). Values at or beyond a limit are
    censored."""
    if u.ndim > y.ndim:
        y = y[:, None]
    z = (y - u) / sigma_e
    out = -0.5 * z ** 2 - LOG_SQRT_2PI - math.log(sigma_e)
    if ll is not None:
        out = np.where(y <= ll, special.log_ndtr((ll - u) / sigma_e), out)
    if ul is not None:
        out = np.where(y >= ul, special.log_ndtr(-(ul - u) / sigma_e), out)
    return out


# ---------------------------------------------------------------------------
# xtnbreg: Hausman, Hall and Griliches (1984)
# ---------------------------------------------------------------------------

def nbreg_re_group_ll(mu: np.ndarray, y: np.ndarray, grp: Groups, ln_r: float, ln_s: float) -> np.ndarray:
    """Random-effects negative binomial: the inverse of one plus the dispersion has a Beta(r, s)
    distribution within the group; mu = exp(xb) is the first parameter of the gamma-Poisson
    mixture of every observation."""
    r, s = math.exp(ln_r), math.exp(ln_s)
    lg = special.gammaln
    sy = grp.sums(y)
    smu = grp.sums(mu)
    term = grp.sums(lg(mu + y) - lg(mu) - lg(y + 1))
    return lg(r + s) + lg(r + smu) + lg(s + sy) - lg(r) - lg(s) - lg(r + s + smu + sy) + term


def nbreg_fe_group_ll(mu: np.ndarray, y: np.ndarray, grp: Groups) -> np.ndarray:
    """Conditional fixed-effects negative binomial: the likelihood given the group totals."""
    lg = special.gammaln
    sy = grp.sums(y)
    smu = grp.sums(mu)
    term = grp.sums(lg(mu + y) - lg(mu) - lg(y + 1))
    return lg(smu) + lg(sy + 1) - lg(smu + sy) + term

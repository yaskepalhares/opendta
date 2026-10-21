"""Exact conditional inference for logistic and Poisson regression ([R] exlogistic, [R] expoisson):
the distribution of the sufficient statistics, and from it the conditional maximum likelihood
estimate, the median unbiased estimate, the exact confidence interval, the exact p-value and the
exact score test. No Session and no output."""

from __future__ import annotations

import math

import numpy as np
from scipy import optimize


def _key(vec) -> tuple:
    return tuple(round(float(v), 9) for v in vec)


def _groups(X: np.ndarray):
    """Distinct rows of X and how many times each occurs (the statistics only depend on the sums)."""
    rows, counts = np.unique(np.round(X, 9), axis=0, return_counts=True)
    return rows, counts


def suff_dist_logistic(X: np.ndarray) -> dict:
    """{t: number of 0/1 vectors y with sum_i y_i x_i = t}. X has the constant first. Identical
    rows are handled together: s of the n_g observations with y = 1 can be chosen in C(n_g, s)
    ways. Counts are rescaled along the way (only ratios are used)."""
    rows, counts = _groups(X)
    k = X.shape[1]
    dist: dict[tuple, float] = {(0.0,) * k: 1.0}
    for xg, ng in zip(rows, counts):
        new: dict[tuple, float] = {}
        for t, c in dist.items():
            for sg in range(int(ng) + 1):
                t1 = _key(np.array(t) + sg * xg)
                new[t1] = new.get(t1, 0.0) + c * math.comb(int(ng), sg)
        mx = max(new.values())
        dist = {t: c / mx for t, c in new.items()}
    return dist


def suff_dist_poisson(X: np.ndarray, total: int) -> dict:
    """{t: sum over the count vectors y with sum y = total and sum_i y_i x_i = t of prod 1/y_i!}.
    X has the constant first (its statistic is sum y). The n_g observations with the same row
    contribute m = their total count with weight n_g^m / m!. The covariates must be integers: the
    distribution lives on an integer lattice kept in a dense array."""
    rows, counts = _groups(X)
    k = X.shape[1]
    R = np.round(rows[:, 1:]).astype(int)
    if not np.allclose(rows[:, 1:], R):
        raise ValueError("noninteger")
    lo = R.min(axis=0) * total
    lo = np.minimum(lo, 0)
    hi = np.maximum(R.max(axis=0) * total, 0)
    shape = (total + 1,) + tuple(int(h - l + 1) for l, h in zip(lo, hi))
    if np.prod(shape, dtype=float) > 6e7:
        raise MemoryError
    A = np.zeros(shape)
    A[(0,) + tuple(int(-l) for l in lo)] = 1.0
    for r, ng in zip(R, counts):
        new = np.zeros(shape)
        lw = np.array([m * math.log(ng) - math.lgamma(m + 1) for m in range(total + 1)])
        lw -= lw.max()                       # a common factor: the distribution is renormalized
        for m in range(total + 1):
            w = float(np.exp(lw[m]))
            shifts = (m,) + tuple(int(m * v) for v in r)
            dst, src = [], []
            for sh, size in zip(shifts, shape):
                if abs(sh) >= size:
                    break
                dst.append(slice(max(sh, 0), size + min(sh, 0)))
                src.append(slice(max(-sh, 0), size - max(sh, 0)))
            else:
                new[tuple(dst)] += w * A[tuple(src)]
        mx = new.max()
        A = new / mx
    out = {}
    last = A[total]
    for idx in np.argwhere(last > 0):
        t = (float(total),) + tuple(float(i + l) for i, l in zip(idx, lo))
        out[t] = float(last[tuple(idx)])
    return out


def conditional(dist: dict, obs: tuple, j: int):
    """Support and weights of T_j given the other statistics at their observed values."""
    sup, wts = [], []
    for t, c in dist.items():
        if all(abs(t[m] - obs[m]) < 1e-7 for m in range(len(obs)) if m != j):
            sup.append(t[j])
            wts.append(c)
    order = np.argsort(sup)
    sup, wts = np.array(sup)[order], np.array(wts)[order]
    return sup, wts / wts.sum()


def _tail(sup, p, beta, t_obs, upper: bool, mid: bool = False) -> float:
    z = np.log(p) + beta * sup
    z -= z.max()
    q = np.exp(z)
    q /= q.sum()
    if upper:
        pr = q[sup >= t_obs - 1e-9].sum()
        eq = q[np.abs(sup - t_obs) < 1e-9].sum()
    else:
        pr = q[sup <= t_obs + 1e-9].sum()
        eq = q[np.abs(sup - t_obs) < 1e-9].sum()
    return float(pr - eq / 2) if mid else float(pr)


def cmle(sup, p, t_obs):
    """Conditional maximum likelihood estimate of the log odds ratio/rate ratio, or None when the
    observed statistic is at the end of its support."""
    if t_obs <= sup.min() + 1e-9 or t_obs >= sup.max() - 1e-9:
        return None

    def f(b):
        z = np.log(p) + b * sup
        z -= z.max()
        q = np.exp(z)
        q /= q.sum()
        return float(q @ sup) - t_obs

    return optimize.brentq(f, -60, 60, xtol=1e-12)


def _solve_tail(sup, p, t_obs, upper: bool, target: float, mid: bool):
    f = lambda b: _tail(sup, p, b, t_obs, upper, mid) - target
    lo, hi = -60.0, 60.0
    if f(lo) * f(hi) > 0:
        return None
    return optimize.brentq(f, lo, hi, xtol=1e-12)


def mue(sup, p, t_obs, mid: bool = False):
    """Median unbiased estimate: the average of the values of beta at which the upper and the lower
    tail probabilities of the observed statistic are one half (the finite one when the other is not)."""
    lo = _solve_tail(sup, p, t_obs, True, 0.5, mid) if t_obs > sup.min() + 1e-9 else None
    hi = _solve_tail(sup, p, t_obs, False, 0.5, mid) if t_obs < sup.max() - 1e-9 else None
    vals = [v for v in (lo, hi) if v is not None]
    return float(np.mean(vals)) if vals else 0.0


def exact_ci(sup, p, t_obs, alpha: float, mid: bool = False):
    """(lower, upper) limits of the exact confidence interval for beta; +-inf at the ends of the support."""
    lo = -np.inf
    hi = np.inf
    if t_obs > sup.min() + 1e-9:
        r = _solve_tail(sup, p, t_obs, True, alpha / 2, mid)
        lo = -np.inf if r is None else r
    if t_obs < sup.max() - 1e-9:
        r = _solve_tail(sup, p, t_obs, False, alpha / 2, mid)
        hi = np.inf if r is None else r
    return float(lo), float(hi)


def exact_p(sup, p, t_obs, mid: bool = False) -> float:
    """2 * the smaller one-sided probability under beta = 0 (at most 1)."""
    up = _tail(sup, p, 0.0, t_obs, True, mid)
    dn = _tail(sup, p, 0.0, t_obs, False, mid)
    return float(min(1.0, 2 * min(up, dn)))


def score_test(dist: dict, obs: tuple, mean_scale: float, X: np.ndarray, kind: str):
    """Exact score test of all slopes zero: the observed statistic and Pr(statistic >= observed)
    under beta = 0. `mean_scale`: p-bar (logistic) or mu-bar (Poisson)."""
    n, k = X.shape
    sx = X.sum(axis=0)
    Xc = X[:, 1:] - X[:, 1:].mean(axis=0)
    var = mean_scale * (1 - mean_scale) if kind == "logistic" else mean_scale
    I = var * (Xc.T @ Xc)
    Ii = np.linalg.pinv(I)
    keys = [t for t in dist if abs(t[0] - obs[0]) < 1e-7]
    T = np.array([t[1:] for t in keys])
    c = np.array([dist[t] for t in keys])
    U = T - mean_scale * sx[1:]
    stat = np.einsum("ij,jk,ik->i", U, Ii, U)
    s_obs = float((np.array(obs[1:]) - mean_scale * sx[1:]) @ Ii @ (np.array(obs[1:]) - mean_scale * sx[1:]))
    p = float(c[stat >= s_obs - 1e-9 * max(1.0, s_obs)].sum() / c.sum())
    return s_obs, p

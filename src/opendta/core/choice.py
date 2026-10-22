"""Choice-model likelihoods on arrays: the conditional logit over groups of alternatives with any
number of chosen ones (clogit, asclogit), the exploded logit for rankings with ties (rologit) and
the GHK simulator of normal rectangle probabilities (mprobit, asmprobit, asroprobit). No Session
and no output; the commands are in estimation/choice.py."""

from __future__ import annotations

import numpy as np
from scipy import special, stats


def group_index(gid: np.ndarray) -> list[np.ndarray]:
    """Row positions of each group, in order of first appearance of the sorted ids."""
    order = np.argsort(gid, kind="stable")
    sg = gid[order]
    cut = np.flatnonzero(np.r_[True, sg[1:] != sg[:-1], True])
    return [order[cut[i]:cut[i + 1]] for i in range(len(cut) - 1)]


def clogit_group(eta: np.ndarray, y: np.ndarray):
    """(ll, d ll / d eta) of one group: Pr(the chosen set) = exp(sum eta chosen) / e_k(exp(eta)),
    with e_k the elementary symmetric polynomial of degree k = number of chosen alternatives."""
    k = int(round(y.sum()))
    n = len(eta)
    m = float(eta.max())
    a = np.exp(eta - m)
    e = np.zeros(k + 1)
    e[0] = 1.0
    for aj in a:
        for t in range(k, 0, -1):
            e[t] += aj * e[t - 1]
    ll = float(y @ eta) - (np.log(e[k]) + k * m)
    d = np.empty(n)
    for j in range(n):
        f_prev = 1.0
        for t in range(1, k):
            f_prev = e[t] - a[j] * f_prev
        d[j] = y[j] - a[j] * f_prev / e[k]
    return ll, d


def clogit(eta: np.ndarray, y: np.ndarray, groups: list[np.ndarray]):
    """Per-observation arrays (ll on the first row of each group, d ll / d eta on every row)."""
    ll = np.zeros(len(eta))
    d = np.zeros(len(eta))
    for idx in groups:
        l, dd = clogit_group(eta[idx], y[idx])
        ll[idx[0]] = l
        d[idx] = dd
    return ll, d


def clogit_null(y: np.ndarray, groups: list[np.ndarray]) -> float:
    """Log likelihood with all coefficients zero: minus the log of the number of ways to choose."""
    tot = 0.0
    for idx in groups:
        n, k = len(idx), int(round(y[idx].sum()))
        tot -= float(special.gammaln(n + 1) - special.gammaln(k + 1) - special.gammaln(n - k + 1))
    return tot


def rologit_group(eta: np.ndarray, rank: np.ndarray, ties: str = "exactm"):
    """(ll, d ll / d eta) of one group. Higher rank = preferred; alternatives with the same rank
    are tied. Unranked alternatives (missing) are dropped by the caller. `ties`: exactm (all
    orderings of a tie set), breslow (the risk set stays whole for the whole tie set), efron
    (as Breslow with the risk set reduced by the average of the tied), none (ties not allowed)."""
    d = np.zeros(len(eta))
    ll = 0.0
    levels = np.unique(rank)[::-1]
    a = np.exp(eta - eta.max())
    shift = float(eta.max())
    remaining = np.ones(len(eta), dtype=bool)
    for lev in levels:
        T = np.flatnonzero(rank == lev)
        m = len(T)
        R = float(a[remaining].sum())
        if m == 1 or ties == "none":
            if m > 1:
                raise ValueError("ties")
            j = T[0]
            ll += eta[j] - (np.log(R) + shift)
            d[j] += 1.0
            d[remaining] -= a[remaining] / R
        elif ties == "breslow":
            ll += float(eta[T].sum()) - m * (np.log(R) + shift)
            d[T] += 1.0
            d[remaining] -= m * a[remaining] / R
        elif ties == "efron":
            aT = float(a[T].sum())
            ll += float(eta[T].sum())
            for r in range(m):
                den = R - r / m * aT
                ll -= np.log(den) + shift
                d[remaining] -= a[remaining] / den
                d[T] += r / m * aT / den * (a[T] / aT)
            d[T] += 1.0
        else:                                # exactm
            ll_t, d_t = _exact_tie(a, T, remaining, shift, eta)
            ll += ll_t
            d += d_t
        remaining[T] = False
    return float(ll), d


def _exact_tie(a, T, remaining, shift, eta):
    """log of the sum over the m! orderings of T of prod a_j / (risk-set sum), and its gradient
    (central differences: tie sets are small). The ratios do not depend on the common shift."""
    R = float(a[remaining].sum())
    ll = float(np.log(_h_full(a[T], R)))
    d = np.zeros(len(a))
    eps = 1e-5
    for j in np.flatnonzero(remaining):
        up, dn = a.copy(), a.copy()
        up[j] *= np.exp(eps)
        dn[j] *= np.exp(-eps)
        d[j] = (np.log(_h_full(up[T], float(up[remaining].sum())))
                - np.log(_h_full(dn[T], float(dn[remaining].sum())))) / (2 * eps)
    return ll, d


def _h_full(aT: np.ndarray, R: float) -> float:
    m = len(aT)
    size = 1 << m
    sums = np.zeros(size)
    for S in range(1, size):
        low = (S & -S).bit_length() - 1
        sums[S] = sums[S & (S - 1)] + aT[low]
    h = np.zeros(size)
    h[0] = 1.0
    for S in range(1, size):
        tot = 0.0
        for j in range(m):
            if S >> j & 1:
                prev = S ^ (1 << j)
                tot += h[prev] * aT[j] / (R - sums[prev])
        h[S] = tot
    return float(h[size - 1])


def rologit(eta: np.ndarray, rank: np.ndarray, groups: list[np.ndarray], ties: str = "exactm"):
    ll = np.zeros(len(eta))
    d = np.zeros(len(eta))
    for idx in groups:
        l, dd = rologit_group(eta[idx], rank[idx], ties)
        ll[idx[0]] = l
        d[idx] = dd
    return ll, d


# ---------------------------------------------------------------------------
# GHK simulator
# ---------------------------------------------------------------------------

def halton(n: int, dim: int, skip: int = 0) -> np.ndarray:
    """First n points (after skipping `skip`) of the Halton sequence in `dim` dimensions."""
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    out = np.zeros((n, dim))
    for d in range(dim):
        p = primes[d]
        for i in range(n):
            f, r, k = 1.0, 0.0, i + 1 + skip
            while k > 0:
                f /= p
                r += f * (k % p)
                k //= p
            out[i, d] = r
    return out


def ghk_upper(mu: np.ndarray, L: np.ndarray, U: np.ndarray) -> np.ndarray:
    """Pr(z_1 < mu_1, ..., z_m < mu_m) for z = L e, e standard normal, by the GHK simulator with
    uniform draws U (draws x m) shared by all rows. mu is (rows x m); returns the rows' probabilities."""
    rows, m = mu.shape
    draws = U.shape[0]
    p = np.ones((rows, draws))
    e = np.zeros((rows, draws, m))
    for i in range(m):
        s = np.zeros((rows, draws))
        for j in range(i):
            s += L[i, j] * e[:, :, j]
        upper = (mu[:, i][:, None] - s) / L[i, i]
        cdf = stats.norm.cdf(upper)
        p *= cdf
        u = U[None, :, i] * cdf
        e[:, :, i] = stats.norm.ppf(np.clip(u, 1e-300, 1 - 1e-16))
    return p.mean(axis=1)

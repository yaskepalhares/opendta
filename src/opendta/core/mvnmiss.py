"""Multivariate normal data with missing values: patterns, log likelihood and the EM algorithm.

Plain functions on arrays (no Session). `patterns` splits the rows by missing-value pattern;
each pattern is summarized by its observed columns, the number of rows, and the mean and
(divide-by-n) covariance of those columns. The log likelihood of a normal model with mean mu and
covariance Sigma is the sum over patterns of the usual normal log likelihood of the pattern's
summary data on the matching elements of mu and Sigma. `em_mvn` maximizes it by the EM algorithm,
optionally with a block-diagonal covariance (variables in different blocks uncorrelated), which is
what the baseline model of a structural equation model needs."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass
class Pattern:
    idx: np.ndarray        # observed columns
    rows: np.ndarray       # row numbers
    n: int
    mean: np.ndarray       # mean of the observed columns
    cov: np.ndarray        # covariance of the observed columns, divided by n


def patterns(X: np.ndarray) -> list[Pattern]:
    """Rows of X (NaN = missing) grouped by the set of observed columns, in order of first
    appearance. Rows with no observed value are dropped."""
    obs = ~np.isnan(X)
    seen: dict[bytes, list[int]] = {}
    for i in range(X.shape[0]):
        key = obs[i].tobytes()
        seen.setdefault(key, []).append(i)
    out = []
    for key, rows in seen.items():
        idx = np.nonzero(np.frombuffer(key, dtype=bool))[0]
        if len(idx) == 0:
            continue
        r = np.asarray(rows)
        Z = X[np.ix_(r, idx)]
        m = Z.mean(axis=0)
        D = Z - m
        out.append(Pattern(idx, r, len(r), m, D.T @ D / len(r)))
    return out


def loglik(pats: list[Pattern], mu: np.ndarray, Sig: np.ndarray) -> float:
    """Log likelihood of the pattern summaries under N(mu, Sig)."""
    ll = 0.0
    for p in pats:
        S = Sig[np.ix_(p.idx, p.idx)]
        sign, logdet = np.linalg.slogdet(S)
        if sign <= 0:
            return -np.inf
        d = p.mean - mu[p.idx]
        Si = np.linalg.inv(S)
        k = len(p.idx)
        ll += -0.5 * p.n * (k * math.log(2 * math.pi) + logdet + np.trace(p.cov @ Si) + d @ Si @ d)
    return float(ll)


def em_mvn(X: np.ndarray, blocks: list[list[int]] | None = None, maxiter: int = 20000,
           tol: float = 1e-12, start: tuple[np.ndarray, np.ndarray] | None = None):
    """ML estimates of the mean and covariance by EM. `blocks` lists groups of columns that may
    be correlated with each other (default: all columns in one block). Returns
    (mean, cov, loglik, iterations)."""
    n, p = X.shape
    pats = patterns(X)
    rows_used = sum(q.n for q in pats)
    if start is None:
        mu = np.nanmean(X, axis=0)
        var = np.nanvar(X, axis=0)
        Sig = np.diag(np.where(var > 0, var, 1.0))
    else:
        mu, Sig = start[0].copy(), start[1].copy()
    mask = np.ones((p, p), dtype=bool)
    if blocks is not None:
        mask = np.zeros((p, p), dtype=bool)
        for b in blocks:
            mask[np.ix_(b, b)] = True
        for i in range(p):
            mask[i, i] = True
    prev = loglik(pats, mu, Sig)
    it = 0
    for it in range(1, maxiter + 1):
        T1 = np.zeros(p)
        T2 = np.zeros((p, p))
        for q in pats:
            o = q.idx
            m = np.setdiff1d(np.arange(p), o)
            Z = X[np.ix_(q.rows, o)]
            if len(m):
                Soo = Sig[np.ix_(o, o)]
                B = np.linalg.solve(Soo, Sig[np.ix_(o, m)]).T          # Smo Soo^-1
                Xm = mu[m] + (Z - mu[o]) @ B.T
                C = Sig[np.ix_(m, m)] - B @ Sig[np.ix_(o, m)]
            full = np.empty((q.n, p))
            full[:, o] = Z
            if len(m):
                full[:, m] = Xm
            T1 += full.sum(axis=0)
            T2 += full.T @ full
            if len(m):
                T2[np.ix_(m, m)] += q.n * C
        mu = T1 / rows_used
        Sig = T2 / rows_used - np.outer(mu, mu)
        Sig = np.where(mask, Sig, 0.0)
        Sig = (Sig + Sig.T) / 2
        cur = loglik(pats, mu, Sig)
        if abs(cur - prev) < tol * (1 + abs(cur)):
            prev = cur
            break
        prev = cur
    return mu, Sig, prev, it

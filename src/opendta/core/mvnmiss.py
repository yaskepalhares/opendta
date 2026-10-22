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


def da_mvn(Y: np.ndarray, X: np.ndarray, normal, chi2, nimp: int, burnin: int,
           between: int) -> list[np.ndarray]:
    """Data augmentation for the multivariate normal regression Y = X B + E (rows of E with
    covariance Sigma) with missing values in Y (NaN) and the noninformative prior
    |Sigma|^-(p+1)/2. P-step: Sigma from an inverse Wishart with n - k degrees of freedom and
    scale S (the residual sums of squares of the completed data), B from its normal conditional;
    I-step: the missing values of every pattern from their conditional normal distribution.
    `normal(m)` returns m standard normal draws and `chi2(df)` one chi-squared draw (the caller's
    random numbers). Returns `nimp` completed copies of Y, the first after `burnin` iterations
    and the others `between` iterations apart."""
    n, p = Y.shape
    k = X.shape[1]
    miss = np.isnan(Y)
    Yf = np.where(miss, np.nanmean(np.where(np.isfinite(Y), Y, np.nan), axis=0)[None, :], Y)
    XtXi = np.linalg.pinv(X.T @ X)
    Lx = np.linalg.cholesky((XtXi + XtXi.T) / 2 + 1e-14 * np.eye(k))
    keys, inv = np.unique(miss, axis=0, return_inverse=True)
    inv = np.asarray(inv).ravel()
    groups = [(keys[g], np.nonzero(inv == g)[0]) for g in range(len(keys)) if keys[g].any()]
    nu = max(n - k, p + 1)
    out: list[np.ndarray] = []
    total = burnin + (nimp - 1) * between
    for it in range(1, total + 1):
        Bh = XtXi @ X.T @ Yf
        R = Yf - X @ Bh
        S = R.T @ R
        Si = np.linalg.pinv((S + S.T) / 2)
        Ls = np.linalg.cholesky((Si + Si.T) / 2 + 1e-14 * np.eye(p))
        A = np.zeros((p, p))
        for i in range(p):
            A[i, i] = math.sqrt(chi2(nu - i))
            for j in range(i):
                A[i, j] = normal(1)[0]
        La = Ls @ A
        Sig = np.linalg.inv(La @ La.T)
        Sig = (Sig + Sig.T) / 2
        Lsig = np.linalg.cholesky(Sig + 1e-14 * np.eye(p))
        B = Bh + Lx @ normal(k * p).reshape(k, p) @ Lsig.T
        XB = X @ B
        for mk, rows in groups:
            o = np.nonzero(~mk)[0]
            m = np.nonzero(mk)[0]
            if len(o):
                Soo = Sig[np.ix_(o, o)]
                Bc = np.linalg.solve(Soo, Sig[np.ix_(o, m)]).T
                mean = XB[np.ix_(rows, m)] + (Yf[np.ix_(rows, o)] - XB[np.ix_(rows, o)]) @ Bc.T
                C = Sig[np.ix_(m, m)] - Bc @ Sig[np.ix_(o, m)]
            else:
                mean = XB[np.ix_(rows, m)]
                C = Sig[np.ix_(m, m)]
            Lc = np.linalg.cholesky((C + C.T) / 2 + 1e-14 * np.eye(len(m)))
            Z = normal(len(rows) * len(m)).reshape(len(rows), len(m))
            Yf[np.ix_(rows, m)] = mean + Z @ Lc.T
        if it >= burnin and (it - burnin) % max(between, 1) == 0:
            out.append(Yf.copy())
    if burnin == 0 or not out:
        out.insert(0, Yf.copy())
    while len(out) < nimp:
        out.append(out[-1].copy())
    return out[:nimp]

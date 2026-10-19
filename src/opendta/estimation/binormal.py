"""Binormal ROC model ([R] rocfit): the scores of the negatives are N(0, 1) and those of the
positives N(mu, sigma^2) on a latent scale cut at K - 1 increasing thresholds, fit by maximum
likelihood to the 2 x K table of counts. Used by rocfit and by roccomp, binormal."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..lazy import lazy
from . import ml

st = lazy("scipy.stats")


@dataclass
class BinormalFit:
    res: ml.MLResult
    b: np.ndarray         # cuts, mu, lnsigma
    V: np.ndarray         # their covariance
    p0: np.ndarray        # fitted category probabilities of the negatives
    p1: np.ndarray        # and of the positives


def _num_grad_hess(f, t):
    k = len(t)
    h = 1e-4 * np.maximum(1.0, np.abs(t))
    g = np.zeros(k)
    H = np.zeros((k, k))
    f0 = f(t)
    for i in range(k):
        e = np.zeros(k)
        e[i] = h[i]
        g[i] = (f(t + e) - f(t - e)) / (2 * h[i])
        H[i, i] = (f(t + e) - 2 * f0 + f(t - e)) / h[i] ** 2
        for j in range(i):
            e2 = np.zeros(k)
            e2[j] = h[j]
            H[i, j] = H[j, i] = (f(t + e + e2) - f(t + e - e2) - f(t - e + e2) + f(t - e - e2)) / (4 * h[i] * h[j])
    return f0, g, H


def fit(n0: np.ndarray, n1: np.ndarray, maxiter: int = 16000) -> BinormalFit:
    """Fits the model to the counts per category (in increasing order) of the negatives n0
    and the positives n1."""
    K = len(n0)
    Phi = st.norm.cdf

    def probs(theta):
        c = np.cumsum(np.concatenate([[theta[0]], np.exp(theta[1:K - 1])]))   # increasing cuts
        mu, sig = theta[K - 1], math.exp(theta[K])
        c0 = np.concatenate([[-np.inf], c, [np.inf]])
        p0 = np.diff(Phi(c0))
        p1 = np.diff(Phi((c0 - mu) / sig))
        return c, p0, p1

    def ll(theta):
        _, p0, p1 = probs(theta)
        if np.any(p0[n0 > 0] <= 0) or np.any(p1[n1 > 0] <= 0):
            return -np.inf
        return float((n0[n0 > 0] * np.log(p0[n0 > 0])).sum() + (n1[n1 > 0] * np.log(p1[n1 > 0])).sum())

    # start: cuts from the negatives' cumulative shares, mu from the mean shift
    cum0 = np.clip(np.cumsum(n0)[:-1] / n0.sum(), 0.02, 0.98)
    c_start = np.maximum.accumulate(st.norm.ppf(cum0) + np.arange(K - 1) * 1e-3)
    cum1 = np.clip(np.cumsum(n1)[:-1] / n1.sum(), 0.02, 0.98)
    mu0 = float(np.mean(c_start - st.norm.ppf(cum1)))
    theta0 = np.concatenate([[c_start[0]], np.log(np.maximum(np.diff(c_start), 1e-2)), [mu0, 0.0]])

    def fun(t, todo):
        return _num_grad_hess(ll, t)
    res = ml.maximize(fun, theta0, maxiter=maxiter)
    # natural parameters (cuts, mu, sigma) and their covariance by the delta method
    t = res.b
    c, p0, p1 = probs(t)
    J = np.zeros((K + 1, K + 1))
    for i in range(K - 1):
        J[i, 0] = 1.0
        for j in range(1, i + 1):
            J[i, j] = math.exp(t[j])
    J[K - 1, K - 1] = 1.0
    J[K, K] = 1.0              # lnsigma kept on the log scale, as reported
    Vt = np.linalg.inv(-res.H)
    b_nat = np.concatenate([c, [t[K - 1], t[K]]])
    return BinormalFit(res, b_nat, J @ Vt @ J.T, p0, p1)


def indices(mu: float, lnsigma: float, V: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """ROC area Phi(a / sqrt(1 + b^2)), delta(m), d(e) and d(a) of the binormal curve
    (a = mu/sigma, b = 1/sigma) and their standard errors from the 2 x 2 covariance V of
    (mu, lnsigma)."""
    Phi = st.norm.cdf

    def idx(v):
        m_, ls = v[0], v[1]
        sg = math.exp(ls)
        a_, b_ = m_ / sg, 1 / sg
        return np.array([float(Phi(a_ / math.sqrt(1 + b_ * b_))), m_, 2 * a_ / (1 + b_),
                         math.sqrt(2) * a_ / math.sqrt(1 + b_ * b_)])
    v0 = np.array([mu, lnsigma])
    G = np.zeros((4, 2))
    for j in range(2):
        e = np.zeros(2)
        e[j] = 1e-6
        G[:, j] = (idx(v0 + e) - idx(v0 - e)) / 2e-6
    return idx(v0), np.sqrt(np.diag(G @ V @ G.T))

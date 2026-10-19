"""Binormal ROC model ([R] rocfit): the scores of the negatives are N(0, 1) and those of the
positives N(mu, sigma^2) on a latent scale cut at K - 1 increasing thresholds, fit by maximum
likelihood to the 2 x K table of counts in the parameters intercept (mu/sigma), slope
(1/sigma) and cutpoints. Used by rocfit and by roccomp, binormal."""

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
    b: np.ndarray         # intercept, slope, cuts
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


def _probs(theta, K):
    """Category probabilities of the negatives, Phi(c_k) differences, and of the positives,
    Phi(b c_k - a) differences."""
    Phi = st.norm.cdf
    a, b, c = theta[0], theta[1], np.asarray(theta[2:])
    c0 = np.concatenate([[-np.inf], c, [np.inf]])
    return np.diff(Phi(c0)), np.diff(Phi(b * c0 - a))


def fit(n0: np.ndarray, n1: np.ndarray, maxiter: int = 16000) -> BinormalFit:
    """Fits the model to the counts per category (in increasing order) of the negatives n0
    and the positives n1, in Stata 14's parameters: intercept a = mu/sigma, slope
    b = 1/sigma and the cutpoints. Starting values: the least-squares line through the
    probit-transformed empirical ROC points (compat 1402: the iteration 0 log likelihood
    matches)."""
    K = len(n0)

    def ll(theta):
        p0, p1 = _probs(theta, K)
        if np.any(p0[n0 > 0] <= 0) or np.any(p1[n1 > 0] <= 0):
            return -np.inf
        return float((n0[n0 > 0] * np.log(p0[n0 > 0])).sum() + (n1[n1 > 0] * np.log(p1[n1 > 0])).sum())

    cum0 = np.clip(np.cumsum(n0)[:-1] / n0.sum(), 1e-4, 1 - 1e-4)
    cum1 = np.clip(np.cumsum(n1)[:-1] / n1.sum(), 1e-4, 1 - 1e-4)
    c = st.norm.ppf(cum0)
    z1 = st.norm.ppf(cum1)
    slope, icpt = np.polyfit(c, z1, 1)
    theta0 = np.concatenate([[-icpt, slope], c])

    def fun(t, todo):
        return _num_grad_hess(ll, t)
    res = ml.maximize(fun, theta0, maxiter=maxiter)
    p0, p1 = _probs(res.b, K)
    V = np.linalg.inv(-res.H)
    return BinormalFit(res, res.b.copy(), (V + V.T) / 2, p0, p1)


def indices(a: float, b: float, V: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """ROC area Phi(a / sqrt(1 + b^2)), delta(m) = a/b, d(e) = 2a/(1 + b) and
    d(a) = sqrt(2) a / sqrt(1 + b^2), with standard errors by the delta method from the
    2 x 2 covariance V of (a, b)."""
    Phi = st.norm.cdf

    def idx(v):
        a_, b_ = v
        return np.array([float(Phi(a_ / math.sqrt(1 + b_ * b_))), a_ / b_, 2 * a_ / (1 + b_),
                         math.sqrt(2) * a_ / math.sqrt(1 + b_ * b_)])
    v0 = np.array([a, b])
    G = np.zeros((4, 2))
    for j in range(2):
        e = np.zeros(2)
        e[j] = 1e-6
        G[:, j] = (idx(v0 + e) - idx(v0 - e)) / 2e-6
    return idx(v0), np.sqrt(np.diag(G @ V @ G.T))

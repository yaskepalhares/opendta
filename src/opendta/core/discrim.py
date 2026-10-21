"""Computations behind discriminant analysis ([MV] discrim lda, qda, knn and [MV] candisc):
group statistics, posterior probabilities, classification, leave-one-out classification, error
rates and canonical discriminant functions.

Plain functions on arrays and numbers: no Session, no variables, no printing."""

from __future__ import annotations

import numpy as np

from . import cluster as K


def group_stats(X: np.ndarray, codes: np.ndarray, g: int, w: np.ndarray | None = None) -> dict:
    """Counts (sum of weights), means and covariance matrices of each group, and the pooled
    within-group covariance matrix (divisor n - g)."""
    n, p = X.shape
    if w is None:
        w = np.ones(n)
    cnt = np.zeros(g)
    means = np.zeros((g, p))
    covs = []
    for k in range(g):
        sel = codes == k
        wk = w[sel]
        cnt[k] = wk.sum()
        means[k] = (wk[:, None] * X[sel]).sum(0) / cnt[k]
    W = np.zeros((p, p))
    for k in range(g):
        sel = codes == k
        Xc = X[sel] - means[k]
        SS = (Xc * w[sel][:, None]).T @ Xc
        W += SS
        covs.append(SS / (cnt[k] - 1.0) if cnt[k] > 1 else np.full((p, p), np.nan))
    return {"n": cnt, "means": means, "covs": covs, "pooled": W / (cnt.sum() - g)}


def priors_from(spec, cnt: np.ndarray) -> np.ndarray:
    """Prior probabilities: 'equal', 'proportional' or a vector."""
    if isinstance(spec, str):
        if spec == "equal":
            return np.full(len(cnt), 1.0 / len(cnt))
        if spec == "proportional":
            return cnt / cnt.sum()
        raise ValueError(spec)
    v = np.asarray(spec, dtype=float)
    return v / v.sum()


def log_posteriors(score: np.ndarray) -> np.ndarray:
    """Posterior probabilities from log(prior) + log-density terms (rows: observations)."""
    m = score.max(1, keepdims=True)
    e = np.exp(score - m)
    return e / e.sum(1, keepdims=True)


def mahalanobis_sq(X: np.ndarray, means: np.ndarray, S: np.ndarray) -> np.ndarray:
    """Squared Mahalanobis distance of each row of X to each group mean (columns)."""
    Si = np.linalg.inv(S)
    out = np.empty((X.shape[0], means.shape[0]))
    for k in range(means.shape[0]):
        d = X - means[k]
        out[:, k] = np.einsum("ij,jk,ik->i", d, Si, d)
    return out


def lda_scores(X: np.ndarray, st: dict, priors: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(log prior - Mahalanobis^2/2 for every group, squared Mahalanobis distances)."""
    d2 = mahalanobis_sq(X, st["means"], st["pooled"])
    return np.log(priors)[None, :] - 0.5 * d2, d2


def qda_scores(X: np.ndarray, st: dict, priors: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    g = st["means"].shape[0]
    d2 = np.empty((X.shape[0], g))
    sc = np.empty_like(d2)
    for k in range(g):
        S = st["covs"][k]
        d = X - st["means"][k]
        d2[:, k] = np.einsum("ij,jk,ik->i", d, np.linalg.inv(S), d)
        sc[:, k] = np.log(priors[k]) - 0.5 * np.linalg.slogdet(S)[1] - 0.5 * d2[:, k]
    return sc, d2


def fit_scores(kind: str, X: np.ndarray, codes: np.ndarray, g: int, w, priors_spec, Xnew: np.ndarray):
    """Fits on (X, codes) and scores Xnew. Returns (log scores, d2, stats, priors)."""
    st = group_stats(X, codes, g, w)
    pri = priors_from(priors_spec, st["n"])
    f = lda_scores if kind == "lda" else qda_scores
    sc, d2 = f(Xnew, st, pri)
    return sc, d2, st, pri


def loo_scores(kind: str, X: np.ndarray, codes: np.ndarray, g: int, w, priors_spec) -> np.ndarray:
    """Leave-one-out log scores: observation i is classified by a rule fitted without it (the
    priors are those of the full sample, as Stata does for the default)."""
    n = X.shape[0]
    out = np.empty((n, g))
    full = group_stats(X, codes, g, w)
    pri = priors_from(priors_spec, full["n"])
    for i in range(n):
        keep = np.ones(n, dtype=bool)
        keep[i] = False
        st = group_stats(X[keep], codes[keep], g, None if w is None else w[keep])
        f = lda_scores if kind == "lda" else qda_scores
        sc, _ = f(X[i:i + 1], st, pri)
        out[i] = sc[0]
    return out


def knn_posteriors(X: np.ndarray, codes: np.ndarray, g: int, k: int, priors: np.ndarray, Xnew: np.ndarray,
                   measure: str = "L2", power: float | None = None, loo: bool = False) -> np.ndarray:
    """Posterior probabilities of the k-nearest-neighbors rule: proportional to
    prior_g * (neighbors in g / observations in g). With `loo` the observation itself is not
    one of its neighbors."""
    D = K.pairwise(Xnew, X, measure, power)
    if loo:
        np.fill_diagonal(D, np.inf)
    n_g = np.bincount(codes, minlength=g).astype(float)
    P = np.zeros((Xnew.shape[0], g))
    for i in range(Xnew.shape[0]):
        idx = np.argsort(D[i], kind="stable")[:k]
        cnt = np.bincount(codes[idx], minlength=g).astype(float)
        P[i] = priors * cnt / n_g
    tot = P.sum(1, keepdims=True)
    tot[tot == 0] = 1.0
    return P / tot


def classify(P: np.ndarray) -> np.ndarray:
    """Group of the largest posterior probability (the first one in a tie)."""
    return np.argmax(P, axis=1)


def class_table(true: np.ndarray, pred: np.ndarray, g: int, w: np.ndarray | None = None) -> np.ndarray:
    T = np.zeros((g, g))
    ww = np.ones(len(true)) if w is None else w
    np.add.at(T, (true, pred), ww)
    return T


def error_rates(T: np.ndarray, priors: np.ndarray) -> tuple[np.ndarray, float]:
    """Error rate of each true group (proportion misclassified) and the prior-weighted total."""
    n_g = T.sum(1)
    err = 1.0 - np.diag(T) / n_g
    return err, float((err * priors).sum())


def smoothed_error_rates(P: np.ndarray, codes: np.ndarray, g: int, priors: np.ndarray) -> tuple[np.ndarray, float]:
    """Posterior-probability (smoothed) error-rate estimates: the average over the observations
    of a group of 1 - P(own group | x)."""
    err = np.array([1.0 - P[codes == k, k].mean() for k in range(g)])
    return err, float((err * priors).sum())


def candisc(st: dict) -> dict:
    """Canonical discriminant functions from the group statistics: eigenvalues of W^-1 B, raw
    coefficients scaled to unit pooled within-group variance, and the canonical correlations."""
    n = st["n"]
    g, p = st["means"].shape
    N = n.sum()
    grand = (n[:, None] * st["means"]).sum(0) / N
    B = sum(n[k] * np.outer(st["means"][k] - grand, st["means"][k] - grand) for k in range(g))
    W = st["pooled"] * (N - g)
    L = np.linalg.cholesky(W)
    Li = np.linalg.inv(L)
    A = Li @ B @ Li.T
    ev, V = np.linalg.eigh((A + A.T) / 2)
    order = np.argsort(-ev)
    ev, V = ev[order], V[:, order]
    m = min(g - 1, p)
    ev, V = np.maximum(ev[:m], 0.0), V[:, :m]
    coef = (Li.T @ V) * np.sqrt(N - g)            # unit pooled within-group variance
    return {"ev": ev, "coef": coef, "canon_corr": np.sqrt(ev / (1.0 + ev)), "B": B, "W": W, "grand": grand}

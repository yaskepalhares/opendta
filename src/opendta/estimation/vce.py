"""Variance of maximum-likelihood estimates: the model (OIM) variance and the robust
sandwich ([R] vce_option, [U] 20.22 Obtaining robust variance estimates).

One implementation for every ML command (logit, probit, poisson, the models of
models2.py, stcox and streg), so that weights and clusters are treated the same way
everywhere.
"""

from __future__ import annotations

import numpy as np


def bread(H: np.ndarray) -> np.ndarray:
    """(-H)^-1, symmetric; the pseudo-inverse when -H is singular (collinear or flat directions)."""
    if H.size == 0:
        return np.zeros((0, 0))
    try:
        V = np.linalg.inv(-H)
    except np.linalg.LinAlgError:
        V = np.linalg.pinv(-H)
    return (V + V.T) / 2


def sandwich(V: np.ndarray, scores: np.ndarray, w: np.ndarray, *, fweight: bool, n: float,
             clusters: np.ndarray | None = None) -> tuple[np.ndarray, int | None]:
    """Robust variance V M V with the small-sample factor; returns (variance, number of clusters).

    scores: derivative of each observation's log likelihood (n x k, unweighted); w: the weights
    used in the fit. With clusters, M sums the weighted scores within each cluster and the factor
    is G/(G-1). Otherwise M = sum w^2 s s' (frequency weights: sum w s s', each observation
    counting w times) and the factor is n/(n-1), n being the reported number of observations
    (the sum of the frequency weights)."""
    ws = scores * w[:, None]
    n_clust = None
    if clusters is not None:
        uq, inv = np.unique(clusters, return_inverse=True)
        n_clust = len(uq)
        G = np.zeros((n_clust, scores.shape[1]))
        np.add.at(G, inv, ws)
        meat = G.T @ G * n_clust / (n_clust - 1)
    else:
        # VERIFY: fweight as expanded data in every ML command (compat 0514)
        meat = ws.T @ (scores if fweight else ws) * n / (n - 1)
    out = V @ meat @ V
    return (out + out.T) / 2, n_clust

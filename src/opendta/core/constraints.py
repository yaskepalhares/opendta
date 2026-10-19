"""Linear constraints on coefficients, R b' = r ([P] makecns): the matrices of the transformation
that estimation under constraints uses, b = bc T' + a."""

from __future__ import annotations

import numpy as np


def constraint_transform(C: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """From C = (R, r), with one row per constraint, T (p x (p - rank R)) and a (1 x p) such that
    every b = bc T' + a satisfies R b' = r. T has the orthonormal eigenvectors of the projector
    I - R'(RR')^-1 R for the eigenvalue 1."""
    R, r = C[:, :-1], C[:, -1]
    p = R.shape[1]
    rank = int(np.linalg.matrix_rank(R)) if R.size else 0
    RRt_inv = np.linalg.pinv(R @ R.T) if R.size else np.zeros((0, 0))
    P = np.eye(p) - R.T @ RRt_inv @ R
    _, vectors = np.linalg.eigh((P + P.T) / 2)
    k = p - rank
    T = vectors[:, p - k:][:, ::-1] if k else np.zeros((p, 0))
    a = R.T @ RRt_inv @ r
    return T, a.reshape(1, p)

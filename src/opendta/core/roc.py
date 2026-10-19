"""Statistics of a nonparametric ROC curve ([R] roctab, Methods and formulas): Bamber's
standard error of the area, the exact binomial interval of the area, and the Lorenz curve
with the Gini and Pietra indices. pos holds the scores of the positives (abnormal, X), neg
those of the negatives (normal, Y)."""

from __future__ import annotations

import math

import numpy as np

from ..lazy import lazy

st = lazy("scipy.stats")


def _below_above(x: np.ndarray, other: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """For each x: how many of other are strictly below it and strictly above it."""
    o = np.sort(other)
    return np.searchsorted(o, x, "left").astype(float), (len(o) - np.searchsorted(o, x, "right")).astype(float)


def se_bamber(pos: np.ndarray, neg: np.ndarray, area: float) -> float:
    """Bamber (1975): var = [p(X != Y) + (m-1) B_XXY + (n-1) B_YYX - 4 (m+n-1)(A - 1/2)^2]
    / (4 (m-1)(n-1)), with
    B_YYX = p(Y_j, Y_k < X_i) + p(X_i < Y_j, Y_k) - 2 p(Y_j < X_i < Y_k) and
    B_XXY = p(X_i, X_l < Y_j) + p(Y_j < X_i, X_l) - 2 p(X_i < Y_j < X_l),
    the probabilities taken over distinct pairs j != k (l != i)."""
    m, n = len(pos), len(neg)
    if m < 2 or n < 2:
        return float("nan")
    lo, hi = _below_above(pos, neg)          # Y below / above each X
    b_yyx = ((lo * (lo - 1)).sum() + (hi * (hi - 1)).sum() - 2 * (lo * hi).sum()) / (m * n * (n - 1))
    lo2, hi2 = _below_above(neg, pos)        # X below / above each Y
    b_xxy = ((lo2 * (lo2 - 1)).sum() + (hi2 * (hi2 - 1)).sum() - 2 * (lo2 * hi2).sum()) / (n * m * (m - 1))
    p_ne = float((lo + hi).sum()) / (m * n)
    var = (p_ne + (m - 1) * b_xxy + (n - 1) * b_yyx - 4 * (m + n - 1) * (area - 0.5) ** 2) / (4 * (m - 1) * (n - 1))
    return math.sqrt(max(var, 0.0))


def binomial_ci(area: float, n: int, level: float) -> tuple[float, float]:
    """Exact (Clopper-Pearson) limits of a proportion equal to the area with n trials, as in
    [R] ci; the number of successes is area*n rounded (VERIFY: the manual's example matches
    with 109 observations and 97 successes)."""
    k = round(area * n)
    alpha = (1 - level / 100) / 2
    lo = 0.0 if k <= 0 else float(st.beta.ppf(alpha, k, n - k + 1))
    hi = 1.0 if k >= n else float(st.beta.ppf(1 - alpha, k + 1, n - k))
    return lo, hi


def lorenz(pos: np.ndarray, neg: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Lorenz curve: the cumulative share of the negatives (x) and of the positives (y) at or
    below each score, from (0, 0) to (1, 1); the Gini index (twice the area between the curve
    and the diagonal) and the Pietra index (twice the largest triangle inscribed between them:
    the largest vertical distance x - y)."""
    cuts = np.unique(np.concatenate([pos, neg]))
    x = np.concatenate([[0.0], np.searchsorted(np.sort(neg), cuts, "right") / len(neg)])
    y = np.concatenate([[0.0], np.searchsorted(np.sort(pos), cuts, "right") / len(pos)])
    under = float(((x[1:] - x[:-1]) * (y[1:] + y[:-1]) / 2).sum())
    return x, y, 2 * (0.5 - under), float((x - y).max())

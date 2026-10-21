"""Computations behind cluster analysis ([MV] cluster, [MV] matrix dissimilarity): dissimilarity
measures, agglomerative linkage by the Lance-Williams recurrence, k-means and k-medians, cutting
a tree, stopping rules and the leaf order of a dendrogram.

Plain functions on arrays and numbers: no Session, no variables, no printing. Ties between equal
dissimilarities are broken by the smallest pair (row-major), which the manuals do not specify."""

from __future__ import annotations

import numpy as np

MEASURES = ("L2", "L2squared", "L1", "Linf", "Canberra", "correlation", "angular")


def pairwise(X: np.ndarray, Y: np.ndarray | None, measure: str, power: float | None = None) -> np.ndarray:
    """Dissimilarities between the rows of X and the rows of Y (X itself when Y is None).
    The similarity measures (correlation, angular) become dissimilarities 1 - s (VERIFY)."""
    Y = X if Y is None else Y
    D = X[:, None, :] - Y[None, :, :]
    m = measure
    if m == "L2":
        return np.sqrt((D ** 2).sum(2))
    if m == "L2squared":
        return (D ** 2).sum(2)
    if m == "L1":
        return np.abs(D).sum(2)
    if m == "Linf":
        return np.abs(D).max(2)
    if m == "Lp":
        return (np.abs(D) ** power).sum(2) ** (1.0 / power)
    if m == "Canberra":
        den = np.abs(X)[:, None, :] + np.abs(Y)[None, :, :]
        with np.errstate(invalid="ignore", divide="ignore"):
            t = np.where(den > 0, np.abs(D) / den, 0.0)
        return t.sum(2)
    if m == "correlation":
        Xc = X - X.mean(1, keepdims=True)
        Yc = Y - Y.mean(1, keepdims=True)
        s = (Xc @ Yc.T) / np.outer(np.sqrt((Xc ** 2).sum(1)), np.sqrt((Yc ** 2).sum(1)))
        return 1.0 - s
    if m == "angular":
        s = (X @ Y.T) / np.outer(np.sqrt((X ** 2).sum(1)), np.sqrt((Y ** 2).sum(1)))
        return 1.0 - s
    raise ValueError(measure)


def linkage(D: np.ndarray, method: str) -> list[tuple[int, int, float]]:
    """Agglomerative clustering of a symmetric dissimilarity matrix. Returns the merges as
    (id of cluster a, id of cluster b, height) with a < b; leaves are 0..n-1 and the cluster
    formed at merge r has id n + r."""
    n = D.shape[0]
    D = D.astype(float).copy()
    np.fill_diagonal(D, np.inf)
    size = np.ones(n)
    ids = list(range(n))
    active = np.ones(n, dtype=bool)
    merges = []
    for r in range(n - 1):
        sub = np.where(active[:, None] & active[None, :], D, np.inf)
        flat = int(np.argmin(sub))
        i, j = divmod(flat, n)
        if i > j:
            i, j = j, i
        h = float(D[i, j])
        ni, nj = size[i], size[j]
        for k in np.flatnonzero(active):
            if k == i or k == j:
                continue
            dki, dkj, dij, nk = D[k, i], D[k, j], D[i, j], size[k]
            if method == "single":
                d = min(dki, dkj)
            elif method == "complete":
                d = max(dki, dkj)
            elif method == "average":
                d = (ni * dki + nj * dkj) / (ni + nj)
            elif method == "waverage":
                d = 0.5 * dki + 0.5 * dkj
            elif method == "centroid":
                d = (ni * dki + nj * dkj) / (ni + nj) - ni * nj * dij / (ni + nj) ** 2
            elif method == "median":
                d = 0.5 * dki + 0.5 * dkj - 0.25 * dij
            elif method == "ward":
                t = ni + nj + nk
                d = ((ni + nk) * dki + (nj + nk) * dkj - nk * dij) / t
            else:
                raise ValueError(method)
            D[k, i] = D[i, k] = d
        merges.append((min(ids[i], ids[j]), max(ids[i], ids[j]), h))
        size[i] = ni + nj
        ids[i] = n + r
        active[j] = False
        D[j, :] = np.inf
        D[:, j] = np.inf
    return merges


def cut_tree(merges: list[tuple[int, int, float]], n: int, groups: int) -> np.ndarray:
    """Group labels (0-based, numbered by first appearance) after merging down to `groups` clusters."""
    parent = list(range(2 * n - 1))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for r, (a, b, _) in enumerate(merges[: n - groups]):
        parent[find(a)] = n + r
        parent[find(b)] = n + r
    roots = [find(i) for i in range(n)]
    lab, out = {}, np.zeros(n, dtype=int)
    for i, rt in enumerate(roots):
        out[i] = lab.setdefault(rt, len(lab))
    return out


def cut_height(merges: list[tuple[int, int, float]], n: int, height: float) -> np.ndarray:
    """Group labels after performing every merge with height <= `height`."""
    done = sum(1 for m in merges if m[2] <= height)
    return cut_tree(merges, n, n - done)


def leaf_order(merges: list[tuple[int, int, float]], n: int) -> list[int]:
    """Left-to-right order of the leaves of the dendrogram."""
    kids = {n + r: (a, b) for r, (a, b, _) in enumerate(merges)}

    def walk(c):
        if c < n:
            return [c]
        a, b = kids[c]
        return walk(a) + walk(b)

    return walk(n + len(merges) - 1) if merges else [0]


def kmeans(X: np.ndarray, centers: np.ndarray, measure: str = "L2", median: bool = False,
           maxiter: int = 10000) -> tuple[np.ndarray, np.ndarray, int]:
    """Lloyd's algorithm starting from `centers`: assign every row to the nearest center, replace
    each center by the mean (k-means) or the coordinate-wise median (k-medians) of its group, and
    repeat until the assignment stops changing. Returns (labels, centers, iterations)."""
    k = centers.shape[0]
    labels = None
    it = 0
    for it in range(1, maxiter + 1):
        D = pairwise(X, centers, measure)
        new = np.argmin(D, axis=1)
        if labels is not None and np.array_equal(new, labels):
            it -= 1
            break
        labels = new
        for g in range(k):
            sel = labels == g
            if sel.any():
                centers[g] = np.median(X[sel], axis=0) if median else X[sel].mean(0)
    return labels, centers, it


def calinski_harabasz(X: np.ndarray, labels: np.ndarray) -> float:
    n = X.shape[0]
    g = len(np.unique(labels))
    mean = X.mean(0)
    B = W = 0.0
    for u in np.unique(labels):
        Xg = X[labels == u]
        c = Xg.mean(0)
        B += len(Xg) * float(((c - mean) ** 2).sum())
        W += float(((Xg - c) ** 2).sum())
    return (B / (g - 1.0)) / (W / (n - g)) if g > 1 and W > 0 else float("nan")


def duda_hart(X: np.ndarray, merges: list[tuple[int, int, float]], groups: int) -> tuple[float, float]:
    """(Je(2)/Je(1), pseudo T-squared) for the split that goes from `groups - 1` to `groups`
    clusters, that is, the undoing of the merge n - groups."""
    n = X.shape[0]
    r = n - groups
    a, b, _ = merges[r]
    kids = {n + i: (x, y) for i, (x, y, _) in enumerate(merges)}

    def members(c):
        if c < n:
            return [c]
        x, y = kids[c]
        return members(x) + members(y)

    ma, mb = members(a), members(b)
    A, B = X[ma], X[mb]
    U = np.vstack([A, B])
    je1 = float(((U - U.mean(0)) ** 2).sum())
    je2 = float(((A - A.mean(0)) ** 2).sum() + ((B - B.mean(0)) ** 2).sum())
    ratio = je2 / je1 if je1 > 0 else float("nan")
    t2 = (1.0 / ratio - 1.0) * (len(U) - 2.0) if ratio > 0 else float("nan")
    return ratio, t2

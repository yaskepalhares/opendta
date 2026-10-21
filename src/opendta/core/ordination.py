"""Computations behind correspondence analysis (ca, mca), multidimensional scaling (mds, mdsmat),
procrustes analysis and biplot.

Plain functions on arrays and numbers: no Session, no variables, no printing. Sign conventions of
the singular vectors are not documented; here the first nonzero element of each row-coordinate
(or object-coordinate) vector is positive (VERIFY)."""

from __future__ import annotations

import numpy as np


def _fix_signs(A: np.ndarray, *others: np.ndarray) -> None:
    """Flips, in place, every column of A (and of `others`) whose first nonzero element is negative."""
    for k in range(A.shape[1]):
        nz = np.flatnonzero(np.abs(A[:, k]) > 1e-12)
        if len(nz) and A[nz[0], k] < 0:
            A[:, k] *= -1.0
            for o in others:
                o[:, k] *= -1.0


# ---------------------------------------------------------------------------
# correspondence analysis
# ---------------------------------------------------------------------------

def ca_fit(N: np.ndarray) -> dict:
    """Simple correspondence analysis of a contingency table. Returns the masses r and c, the
    singular values, the principal inertias (squared singular values), the total inertia and
    the standard coordinates A (rows) and B (columns) of every nontrivial dimension."""
    n = float(N.sum())
    P = N / n
    r, c = P.sum(1), P.sum(0)
    S = (P - np.outer(r, c)) / np.sqrt(np.outer(r, c))
    U, sv, Vt = np.linalg.svd(S, full_matrices=False)
    k = int(np.sum(sv > 1e-10 * max(sv[0], 1e-300)))
    U, sv, V = U[:, :k], sv[:k], Vt.T[:, :k]
    A, B = U / np.sqrt(r)[:, None], V / np.sqrt(c)[:, None]
    _fix_signs(A, B)
    return {"n": n, "r": r, "c": c, "sv": sv, "inertia": sv ** 2, "total": float(np.sum(sv ** 2)), "A": A, "B": B}


def ca_coordinates(fit: dict, normalize: str = "symmetric") -> tuple[np.ndarray, np.ndarray]:
    """Row and column coordinates in the chosen normalization: symmetric, principal, row, column
    or standard."""
    sv, A, B = fit["sv"], fit["A"], fit["B"]
    if normalize == "symmetric":
        return A * np.sqrt(sv), B * np.sqrt(sv)
    if normalize == "principal":
        return A * sv, B * sv
    if normalize == "row":
        return A * sv, B
    if normalize == "column":
        return A, B * sv
    if normalize == "standard":
        return A, B
    raise ValueError(normalize)


def ca_point_stats(mass: np.ndarray, std: np.ndarray, sv: np.ndarray, dims: int) -> dict:
    """Statistics of the points of one side (rows or columns): share of the total inertia,
    quality (squared correlation with the retained dimensions), and for every retained dimension
    the squared correlation and the contribution to the dimension's inertia."""
    prin = std * sv
    d2 = (prin ** 2).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        sqcorr = np.where(d2[:, None] > 0, prin ** 2 / d2[:, None], 0.0)
    inert = mass * d2
    return {"mass": mass, "inertia": inert / inert.sum() if inert.sum() > 0 else inert,
            "quality": sqcorr[:, :dims].sum(1), "sqcorr": sqcorr[:, :dims], "contrib": (mass[:, None] * std ** 2)[:, :dims]}


# ---------------------------------------------------------------------------
# multiple correspondence analysis
# ---------------------------------------------------------------------------

def mca_fit(codes: list[np.ndarray], sizes: list[int]) -> dict:
    """MCA through the correspondence analysis of the indicator matrix. `codes[q]` holds the
    0-based category of every individual for variable q. Returns the CA of the indicator matrix,
    the Benzecri-adjusted principal inertias (VERIFY) and the adjusted total inertia."""
    n, Q = len(codes[0]), len(codes)
    Z = np.zeros((n, sum(sizes)))
    off = 0
    for q in range(Q):
        Z[np.arange(n), off + codes[q]] = 1.0
        off += sizes[q]
    fit = ca_fit(Z)
    lam = fit["inertia"]
    J = sum(sizes)
    burt_total = float(np.sum(lam ** 2))
    keep = lam > 1.0 / Q
    adj = np.where(keep, (Q / (Q - 1.0) * (lam - 1.0 / Q)) ** 2, 0.0)
    fit.update({"Q": Q, "J": J, "adjusted": adj, "adj_total": Q / (Q - 1.0) * (burt_total - (J - Q) / Q ** 2)})
    return fit


# ---------------------------------------------------------------------------
# multidimensional scaling
# ---------------------------------------------------------------------------

def mds_classical(D: np.ndarray, dims: int) -> dict:
    """Classical (Torgerson) scaling of a symmetric dissimilarity matrix: eigenvalues of the
    double-centered squared dissimilarities and the principal coordinates."""
    n = D.shape[0]
    J = np.eye(n) - 1.0 / n
    B = -0.5 * J @ (D ** 2) @ J
    w, V = np.linalg.eigh((B + B.T) / 2)
    order = np.argsort(-w)
    w, V = w[order], V[:, order]
    dims = min(dims, int(np.sum(w > 0)))
    Y = V[:, :dims] * np.sqrt(w[:dims])
    _fix_signs(Y)
    return {"ev": w, "Y": Y, "dims": dims}


def mardia(ev: np.ndarray, dims: int) -> tuple[float, float]:
    """Mardia fit measures of a classical solution with `dims` dimensions."""
    return float(ev[:dims].sum() / np.abs(ev).sum()), float((ev[:dims] ** 2).sum() / (ev ** 2).sum())


def stress1(D: np.ndarray, Y: np.ndarray) -> float:
    """Kruskal stress-1 of the configuration Y for the dissimilarities D."""
    d = np.sqrt(((Y[:, None, :] - Y[None, :, :]) ** 2).sum(2))
    iu = np.triu_indices(D.shape[0], 1)
    return float(np.sqrt(((D[iu] - d[iu]) ** 2).sum() / (D[iu] ** 2).sum()))


def mds_modern(D: np.ndarray, Y0: np.ndarray, maxiter: int = 1000, tol: float = 1e-10) -> tuple[np.ndarray, int]:
    """Metric modern scaling by SMACOF majorization from the start configuration Y0, minimizing
    the raw stress. Returns the configuration (centered) and the iterations used."""
    n = D.shape[0]
    Y = Y0 - Y0.mean(0)
    prev = np.inf
    it = 0
    for it in range(1, maxiter + 1):
        d = np.sqrt(((Y[:, None, :] - Y[None, :, :]) ** 2).sum(2))
        with np.errstate(invalid="ignore", divide="ignore"):
            Bm = np.where(d > 0, -D / d, 0.0)
        np.fill_diagonal(Bm, 0.0)
        np.fill_diagonal(Bm, -Bm.sum(1))
        Y = Bm @ Y / n
        s = float(np.sum(np.triu((D - d) ** 2, 1)))
        if abs(prev - s) < tol * max(prev, 1e-300):
            break
        prev = s
    return Y, it


# ---------------------------------------------------------------------------
# procrustes
# ---------------------------------------------------------------------------

def procrustes_orthogonal(Y: np.ndarray, X: np.ndarray, scale: bool = True) -> dict:
    """Orthogonal Procrustes: the rotation T, the dilation c and the translation that bring the
    source X closest (in least squares) to the target Y."""
    my, mx = Y.mean(0), X.mean(0)
    Yc, Xc = Y - my, X - mx
    U, s, Vt = np.linalg.svd(Xc.T @ Yc, full_matrices=False)
    T = U @ Vt
    c = float(s.sum() / (Xc ** 2).sum()) if scale else 1.0
    fitted = c * Xc @ T + my
    rss = float(((Y - fitted) ** 2).sum())
    ss = float((Yc ** 2).sum())
    return {"T": T, "c": c, "trans": my - c * mx @ T, "fitted": fitted, "rss": rss, "ss": ss}


# ---------------------------------------------------------------------------
# biplot
# ---------------------------------------------------------------------------

def biplot_coords(X: np.ndarray, std: bool = False, alpha: float = 0.5) -> dict:
    """Rank-2 (or higher) biplot: singular value decomposition of the centered (and, with std,
    standardized) data. Rows get U*s^alpha and columns V*s^(1-alpha)."""
    Z = X - X.mean(0)
    if std:
        Z = Z / Z.std(0, ddof=1)
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    rows, cols = U * s ** alpha, Vt.T * s ** (1.0 - alpha)
    _fix_signs(rows, cols)
    return {"rows": rows, "cols": cols, "sv": s}

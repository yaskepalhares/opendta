"""Computations behind principal components, factor analysis, rotation and scale reliability
([MV] pca, [MV] factor, [MV] rotate, [MV] alpha).

Plain functions on arrays and numbers: no Session, no variables, no printing. The formulas come
from the manuals' Methods and formulas and the cited papers; the layout and the choices the
manual does not settle (signs of eigenvectors, default numbers of factors) live in the command
modules and carry a VERIFY mark there."""

from __future__ import annotations

import math

import numpy as np

from ..lazy import lazy

so = lazy("scipy.optimize")

# ---------------------------------------------------------------------------
# correlation and covariance with weights
# ---------------------------------------------------------------------------


def weighted_moments(X: np.ndarray, w: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """(means, covariance with divisor sum(w)-1, correlation, sum of weights) of the rows of X."""
    n = X.shape[0]
    if w is None:
        w = np.ones(n)
    sw = float(w.sum())
    mean = (w[:, None] * X).sum(0) / sw
    Xc = X - mean
    S = (Xc * w[:, None]).T @ Xc / (sw - 1.0)
    sd = np.sqrt(np.diag(S))
    with np.errstate(invalid="ignore", divide="ignore"):
        R = S / np.outer(sd, sd)
    np.fill_diagonal(R, 1.0)
    return mean, S, R, sw


def sorted_eigh(A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Eigenvalues in descending order and the matching unit-length eigenvectors (columns)."""
    w, v = np.linalg.eigh((A + A.T) / 2)
    order = np.argsort(-w, kind="stable")
    return w[order], v[:, order]


def fix_signs(V: np.ndarray) -> np.ndarray:
    """Sign of each column: the sum of its elements is positive (VERIFY against Stata, whose
    eigenvectors have the sign that LAPACK returns)."""
    V = V.copy()
    for j in range(V.shape[1]):
        if V[:, j].sum() < 0:
            V[:, j] = -V[:, j]
    return V


# ---------------------------------------------------------------------------
# principal components
# ---------------------------------------------------------------------------

def pca_fit(C: np.ndarray) -> dict:
    """Eigen decomposition of the correlation or covariance matrix C."""
    ev, V = sorted_eigh(C)
    return {"ev": ev, "V": fix_signs(V), "trace": float(np.trace(C))}


def pca_unexplained(C: np.ndarray, ev: np.ndarray, V: np.ndarray, k: int) -> np.ndarray:
    """Variance of each variable not explained by the first k components."""
    explained = (V[:, :k] ** 2 * ev[:k]).sum(1)
    return np.diag(C) - explained


def n_retained(ev: np.ndarray, components: int | None, mineigen: float | None, default_min: float) -> int:
    if components is not None:
        return max(1, min(int(components), len(ev)))
    lo = default_min if mineigen is None else mineigen
    return int(np.sum(ev > lo)) if mineigen is None else int(np.sum(ev >= lo))


# ---------------------------------------------------------------------------
# factor analysis
# ---------------------------------------------------------------------------

def smc(R: np.ndarray) -> np.ndarray:
    """Squared multiple correlation of each variable with the others: 1 - 1/diag(R^-1)."""
    d = np.diag(np.linalg.inv(R))
    return 1.0 - 1.0 / d


def _loadings_from_reduced(Rr: np.ndarray, k: int | None, mineigen: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ev, V = sorted_eigh(Rr)
    if k is None:
        k = int(np.sum(ev > mineigen))
    k = max(0, min(k, len(ev)))
    L = V[:, :k] * np.sqrt(np.maximum(ev[:k], 0.0))
    return L, ev, V


def factor_pf(R: np.ndarray, k: int | None, mineigen: float = 0.0, psi0: np.ndarray | None = None) -> dict:
    """Principal factors: eigen decomposition of R with the SMC (or `psi0`) on the diagonal."""
    p = R.shape[0]
    h = smc(R) if psi0 is None else 1.0 - psi0
    Rr = R.copy()
    Rr[np.diag_indices(p)] = h
    L, ev, _ = _loadings_from_reduced(Rr, k, mineigen)
    L = flip_columns(L)
    return {"L": L, "ev": ev, "psi": 1.0 - (L ** 2).sum(1), "trace": float(np.trace(Rr)), "iter": 0}


def factor_ipf(R: np.ndarray, k: int | None, mineigen: float = 0.0, tol: float = 1e-6,
               maxiter: int = 1000, heywood: bool = False) -> dict:
    """Iterated principal factors: the communalities are updated from the loadings until they
    stop changing (the maximum change is below `tol`)."""
    p = R.shape[0]
    h = smc(R)
    L = None
    for it in range(1, maxiter + 1):
        Rr = R.copy()
        Rr[np.diag_indices(p)] = h
        L, ev, _ = _loadings_from_reduced(Rr, k, mineigen)
        hn = (L ** 2).sum(1)
        if not heywood:
            hn = np.minimum(hn, 1.0 - 0.005)
        done = np.max(np.abs(hn - h)) < tol
        h = hn
        if done:
            break
    Rr = R.copy()
    Rr[np.diag_indices(p)] = h
    L, ev, _ = _loadings_from_reduced(Rr, k, mineigen)
    L = flip_columns(L)
    return {"L": L, "ev": ev, "psi": 1.0 - (L ** 2).sum(1), "trace": float(np.trace(Rr)), "iter": it}


def factor_pcf(R: np.ndarray, k: int | None, mineigen: float = 1.0) -> dict:
    """Principal-component factors: the loadings are eigenvectors of R scaled by sqrt(eigenvalue)."""
    ev, V = sorted_eigh(R)
    if k is None:
        k = int(np.sum(ev >= mineigen))
    k = max(0, min(k, len(ev)))
    L = flip_columns(V[:, :k] * np.sqrt(np.maximum(ev[:k], 0.0)))
    return {"L": L, "ev": ev, "psi": 1.0 - (L ** 2).sum(1), "trace": float(np.trace(R)), "iter": 0}


def flip_columns(L: np.ndarray) -> np.ndarray:
    """Sign convention of the loadings: the sum of each column is positive (VERIFY)."""
    L = L.copy()
    for j in range(L.shape[1]):
        if L[:, j].sum() < 0:
            L[:, j] = -L[:, j]
    return L


def ml_discrepancy(R: np.ndarray, L: np.ndarray, psi: np.ndarray) -> float:
    S = L @ L.T + np.diag(psi)
    sign, ld = np.linalg.slogdet(S)
    sign0, ld0 = np.linalg.slogdet(R)
    return float(ld + np.trace(np.linalg.solve(S, R)) - ld0 - R.shape[0])


def _ml_loadings(R: np.ndarray, psi: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    d = 1.0 / np.sqrt(psi)
    Rs = R * np.outer(d, d)
    th, V = sorted_eigh(Rs)
    L = (V[:, :k] * np.sqrt(np.maximum(th[:k] - 1.0, 0.0))) / d[:, None]
    return L, th


def factor_ml(R: np.ndarray, k: int, heywood: bool = False, tol: float = 1e-6, maxiter: int = 500,
              floor: float = 0.005) -> dict:
    """Maximum likelihood factor analysis on the correlation matrix: the uniquenesses minimize the
    concentrated discrepancy (Joreskog, 1967); the loadings follow from them. Uniquenesses are
    bounded below by `floor` unless `heywood` allows the boundary to be crossed."""
    p = R.shape[0]
    psi0 = 1.0 - smc(R)
    psi0 = np.clip(psi0, floor, 1.0)

    def f(psi):
        d = 1.0 / np.sqrt(psi)
        th = sorted_eigh(R * np.outer(d, d))[0]
        t = th[k:]
        return float(np.sum(t - np.log(t) - 1.0))

    def g(psi):
        L, _ = _ml_loadings(R, psi, k)
        S = L @ L.T + np.diag(psi)
        Si = np.linalg.inv(S)
        return np.diag(Si @ (S - R) @ Si)

    lo = floor if not heywood else 1e-8
    res = so.minimize(f, psi0, jac=g, method="L-BFGS-B", bounds=[(lo, 1.0)] * p,
                      options={"maxiter": maxiter, "ftol": 1e-14, "gtol": 1e-9})
    psi = res.x
    L, th = _ml_loadings(R, psi, k)
    L = flip_columns(L)
    # factor eigenvalues of the scaled reduced matrix, as reported by Stata for ml
    return {"L": L, "psi": psi, "ev": th, "F": ml_discrepancy(R, L, psi), "iter": int(res.nit),
            "heywood": bool(np.any(psi <= floor * (1 + 1e-6)))}


def lr_independence(R: np.ndarray, n: float) -> tuple[float, int]:
    """LR test of independence vs saturated: -(n - 1 - (2p+5)/6) ln|R| on p(p-1)/2 df."""
    p = R.shape[0]
    sign, ld = np.linalg.slogdet(R)
    return float(-(n - 1.0 - (2 * p + 5) / 6.0) * ld), p * (p - 1) // 2


def lr_factors(F: float, n: float, p: int, k: int) -> tuple[float, int]:
    """LR test of k factors vs saturated (Bartlett's correction) and its degrees of freedom."""
    chi2 = (n - 1.0 - (2 * p + 5) / 6.0 - 2.0 * k / 3.0) * F
    df = ((p - k) ** 2 - (p + k)) // 2
    return float(chi2), int(df)


def kmo(R: np.ndarray) -> tuple[np.ndarray, float]:
    """Kaiser-Meyer-Olkin measure of sampling adequacy: per variable and overall, from the
    correlations and the partial correlations (the anti-image correlations)."""
    Ri = np.linalg.inv(R)
    d = np.sqrt(np.diag(Ri))
    Q = -Ri / np.outer(d, d)
    np.fill_diagonal(Q, 0.0)
    R0 = R.copy()
    np.fill_diagonal(R0, 0.0)
    r2 = (R0 ** 2).sum(1)
    q2 = (Q ** 2).sum(1)
    return r2 / (r2 + q2), float(r2.sum() / (r2.sum() + q2.sum()))


def anti_image(R: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(anti-image correlation matrix, anti-image covariance matrix)."""
    Ri = np.linalg.inv(R)
    d = np.diag(Ri)
    corr = -Ri / np.sqrt(np.outer(d, d))
    np.fill_diagonal(corr, 1.0)
    cov = Ri / np.outer(d, d)
    return corr, cov


# ---------------------------------------------------------------------------
# rotation
# ---------------------------------------------------------------------------

def orthomax(A: np.ndarray, gamma: float, normalize: bool = False, tol: float = 1e-8,
             maxiter: int = 1000) -> tuple[np.ndarray, np.ndarray, int]:
    """Orthogonal rotation of the loadings A (p x k) by the orthomax family: gamma = 1 varimax,
    0 quartimax, k/2 equamax, p(k-1)/(p+k-2) parsimax. Returns (rotated loadings, rotation
    matrix T with rotated = A T, iterations)."""
    p, k = A.shape
    if normalize:
        h = np.sqrt((A ** 2).sum(1))
        h[h == 0] = 1.0
        A = A / h[:, None]
    T = np.eye(k)
    d = 0.0
    it = 0
    for it in range(1, maxiter + 1):
        L = A @ T
        B = A.T @ (L ** 3 - (gamma / p) * L @ np.diag((L ** 2).sum(0)))
        U, S, Vt = np.linalg.svd(B)
        T = U @ Vt
        dn = float(S.sum())
        if dn < d * (1 + tol):
            break
        d = dn
    L = A @ T
    if normalize:
        L = L * h[:, None]
    return L, T, it


def _oblimin_crit(L: np.ndarray, gamma: float) -> tuple[float, np.ndarray]:
    p, k = L.shape
    L2 = L ** 2
    X = L2 @ (np.ones((k, k)) - np.eye(k))
    if gamma != 0:
        X = (np.eye(p) - gamma / p * np.ones((p, p))) @ X
    return float((L2 * X).sum() / 4.0), L * X


def oblimin(A: np.ndarray, gamma: float = 0.0, normalize: bool = False, eps: float = 1e-6,
            maxiter: int = 1000) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Oblique rotation by the oblimin family with gradient projection (Bernaards and Jennrich,
    2005); gamma = 0 is quartimin. Returns (pattern loadings, rotation matrix T, factor
    correlations Phi = T'T, iterations)."""
    p, k = A.shape
    if normalize:
        h = np.sqrt((A ** 2).sum(1))
        h[h == 0] = 1.0
        A = A / h[:, None]
    T = np.eye(k)
    L = A @ np.linalg.inv(T).T
    f, Gq = _oblimin_crit(L, gamma)
    G = -(L.T @ Gq @ np.linalg.inv(T)).T
    al = 1.0
    it = 0
    for it in range(1, maxiter + 1):
        Gp = G - T @ np.diag((T * G).sum(0))
        s = float(np.linalg.norm(Gp))
        if s < eps:
            break
        al *= 2.0
        for _ in range(11):
            X = T - al * Gp
            Tt = X / np.sqrt((X ** 2).sum(0))
            L = A @ np.linalg.inv(Tt).T
            ft, Gqt = _oblimin_crit(L, gamma)
            if ft < f - 0.5 * s * s * al:
                break
            al /= 2.0
        T, f = Tt, ft
        G = -(L.T @ Gqt @ np.linalg.inv(T)).T
    L = A @ np.linalg.inv(T).T
    if normalize:
        L = L * h[:, None]
    return L, T, T.T @ T, it


def promax(A: np.ndarray, power: float = 3.0, normalize: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Oblique promax rotation (Hendrickson and White, 1964): a varimax solution is raised to
    `power` to form a target, and the loadings are rotated towards it by least squares.
    Returns (pattern loadings, rotation matrix, factor correlations)."""
    Lv, Tv, _ = orthomax(A, 1.0, normalize=normalize)
    P = Lv * np.abs(Lv) ** (power - 1.0)
    U = np.linalg.solve(Lv.T @ Lv, Lv.T @ P)
    d = np.diag(np.linalg.inv(U.T @ U))
    U = U * np.sqrt(d)
    L = Lv @ U
    T = Tv @ U
    Phi = np.linalg.inv(U.T @ U)
    return L, T, Phi


# ---------------------------------------------------------------------------
# Cronbach's alpha
# ---------------------------------------------------------------------------

def pairwise_cov(X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Covariance matrix (divisor n-1 within each pair) and the matrix of pair sample sizes,
    using for each pair the rows where both variables are present (NaN marks missing)."""
    k = X.shape[1]
    C = np.full((k, k), np.nan)
    N = np.zeros((k, k))
    for i in range(k):
        for j in range(i, k):
            ok = ~(np.isnan(X[:, i]) | np.isnan(X[:, j]))
            n = int(ok.sum())
            N[i, j] = N[j, i] = n
            if n > 1:
                a, b = X[ok, i], X[ok, j]
                C[i, j] = C[j, i] = float(((a - a.mean()) * (b - b.mean())).sum() / (n - 1))
    return C, N


def alpha_from_cov(C: np.ndarray, standardized: bool = False) -> dict:
    """Scale reliability from a covariance matrix: average interitem covariance, average
    interitem correlation and alpha (covariance form, or standardized when requested)."""
    k = C.shape[0]
    off = ~np.eye(k, dtype=bool)
    cbar = float(C[off].mean())
    v = np.diag(C)
    d = np.sqrt(v)
    Rm = C / np.outer(d, d)
    rbar = float(Rm[off].mean())
    if standardized:
        a = k * rbar / (1.0 + (k - 1.0) * rbar)
    else:
        a = k * cbar / (float(v.mean()) + (k - 1.0) * cbar)
    return {"alpha": float(a), "cbar": cbar, "rbar": rbar, "k": k}

"""Rotation criteria and the gradient projection algorithm ([MV] rotatemat, Methods and formulas).

Each criterion c(L) of the rotated matrix L (p x k) returns its value and its gradient dc/dL. The
criteria are the ones the manual writes with the inner product <A, B> = trace(A'B) and L2 = L.L
(elementwise square). The optimization is the gradient projection algorithm of Jennrich (2001,
2002) as written in Bernaards and Jennrich (2005): orthogonal rotations L = A T with T'T = I;
oblique rotations L = A inv(T') with diag(T'T) = 1.

Plain functions on arrays: no Session, no output.
"""

from __future__ import annotations

import numpy as np

ORTHOGONAL_ONLY = ("vgpf", "quartimax", "equamax", "parsimax", "entropy", "tandem1", "tandem2")
GPF_CRITERIA = ORTHOGONAL_ONLY + ("oblimin", "cf", "bentler", "oblimax", "quartimin", "target", "partial")


def criterion(name: str, L: np.ndarray, arg=None) -> tuple[float, np.ndarray]:
    """Value and gradient of a rotation criterion at L. arg: gamma (oblimin), kappa (cf), the
    target H (target) or the pair (H, W) (partial)."""
    p, k = L.shape
    L2 = L * L
    N = np.ones((k, k)) - np.eye(k)
    if name in ("quartimax",):
        return -float((L2 * L2).sum()) / 4.0, -L ** 3
    if name in ("oblimin", "quartimin", "vgpf", "equamax"):
        gamma = {"quartimin": 0.0, "vgpf": 1.0, "equamax": p / 2.0}.get(name, arg)
        X = L2 @ N
        if gamma:
            X = (np.eye(p) - gamma / p * np.ones((p, p))) @ X
        return float((L2 * X).sum()) / 4.0, L * X
    if name in ("cf", "parsimax"):
        kappa = (k - 1.0) / (p + k - 2.0) if name == "parsimax" else float(arg)
        M = np.ones((p, p)) - np.eye(p)
        X1, X2 = L2 @ N, M @ L2
        f = (1.0 - kappa) * float((L2 * X1).sum()) / 4.0 + kappa * float((L2 * X2).sum()) / 4.0
        return f, (1.0 - kappa) * L * X1 + kappa * L * X2
    if name == "entropy":
        logL2 = np.log(np.where(L2 > 0, L2, 1.0))
        return -float((L2 * logL2).sum()) / 2.0, -(L * logL2 + L)
    if name in ("tandem1", "tandem2"):
        LL = L @ L.T
        LL2 = LL * LL
        G2 = 4.0 * (LL * (L2 @ L2.T)) @ L
        if name == "tandem1":
            Q = LL2 @ L2
            return -float((L2 * Q).sum()), -4.0 * L * Q - G2
        Q = (np.ones((p, p)) - LL2) @ L2
        return float((L2 * Q).sum()), 4.0 * L * Q - G2
    if name == "bentler":
        # -log det{(L2)'L2} + log det[diag{(L2)'L2}]: the manual's criterion with the sign that
        # makes simple structure the minimum; Stata 14 reports it at this scale
        M_ = L2.T @ L2
        D = np.diag(np.diag(M_))
        f = -(np.linalg.slogdet(M_)[1] - np.linalg.slogdet(D)[1])
        return float(f), -4.0 * L * (L2 @ (np.linalg.inv(M_) - np.linalg.inv(D)))
    if name == "oblimax":
        s4, s2 = float((L2 * L2).sum()), float(L2.sum())
        return -(np.log(s4) - 2.0 * np.log(s2)), -(4.0 * L ** 3 / s4 - 4.0 * L / s2)
    if name == "target":
        D = L - arg
        return 0.5 * float((D * D).sum()), D
    if name == "partial":
        H, W = arg
        D = W * (L - H)
        return float((D * D).sum()), 2.0 * W * D
    raise ValueError(name)


def kaiser(A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Kaiser normalization: rows scaled to unit length. Returns (normalized A, row norms)."""
    h = np.sqrt((A ** 2).sum(1))
    h[h == 0] = 1.0
    return A / h[:, None], h


def gpf(A: np.ndarray, name: str, arg=None, oblique: bool = False, T0: np.ndarray | None = None,
        maxiter: int = 1000, tol: float = 1e-6, gtol: float = 1e-6, ltol: float = 1e-6,
        maxstep: int = 20) -> dict:
    """Minimize criterion `name` over the orthogonal (or normal, when oblique) rotations of A,
    starting at T0 (identity by default). Convergence: the norm of the projected gradient is
    below gtol. VERIFY: how the manual's three criteria combine; this rule, with the iteration
    count starting at 1, reproduces the counts of the manual's rotatemat examples (35 and 78),
    and requiring tol and ltol as well does not, so tol and ltol are only reported back.

    Returns {"T", "L" (rotated A), "f", "iter", "converged"}; L = A T (orthogonal) or
    A inv(T') (oblique)."""
    k = A.shape[1]
    T = np.eye(k) if T0 is None else np.array(T0, dtype=float)

    def at(T_):
        return A @ np.linalg.inv(T_).T if oblique else A @ T_

    def grad(T_, L_, Gq_):
        if oblique:
            return -(L_.T @ Gq_ @ np.linalg.inv(T_)).T
        return A.T @ Gq_

    L = at(T)
    f, Gq = criterion(name, L, arg)
    G = grad(T, L, Gq)
    al = 1.0
    it = 0
    converged = False
    for it in range(1, maxiter + 1):
        if oblique:
            Gp = G - T @ np.diag((T * G).sum(0))
        else:
            M = T.T @ G
            Gp = G - T @ ((M + M.T) / 2.0)
        s = float(np.linalg.norm(Gp))
        if s < gtol:
            converged = True
            break
        if it == maxiter:
            break
        al *= 2.0
        for _ in range(maxstep + 1):
            X = T - al * Gp
            if oblique:
                Tt = X / np.sqrt((X ** 2).sum(0))
            else:
                U, _, Vt = np.linalg.svd(X)
                Tt = U @ Vt
            Lt = at(Tt)
            ft, Gqt = criterion(name, Lt, arg)
            if ft < f - 0.5 * s * s * al:
                break
            al /= 2.0
        T, f, L, Gq = Tt, ft, Lt, Gqt
        G = grad(T, L, Gq)
    return {"T": T, "L": L, "f": float(f), "iter": it, "converged": converged}


def random_start(k: int, oblique: bool, Z: np.ndarray) -> np.ndarray:
    """A random orthogonal (or, oblique, normal) k x k matrix built from the standard normal
    deviates Z, for the random and protect() options."""
    if oblique:
        return Z / np.sqrt((Z ** 2).sum(0))
    Q, R_ = np.linalg.qr(Z)
    return Q * np.sign(np.diag(R_))


def orient(L: np.ndarray, T: np.ndarray, sort_by: np.ndarray | None = None):
    """Columns with positive sums ([MV] rotate, technical note), in decreasing order of sort_by
    (default: the column sums of squares of L, the order of the manual's rotatemat examples,
    VERIFY). The rotation matrix follows the same changes; this keeps L = A T and L = A inv(T')."""
    sgn = np.where(L.sum(0) < 0, -1.0, 1.0)
    L, T = L * sgn, T * sgn
    key = (L ** 2).sum(0) if sort_by is None else np.asarray(sort_by)
    order = np.argsort(-key, kind="stable")
    return L[:, order], T[:, order]

"""GMM estimation of dynamic panel-data models: Arellano-Bond, Arellano-Bover/Blundell-Bond
([XT] xtabond, [XT] xtdpd, [XT] xtdpdsys, Methods and formulas).

Plain functions on arrays: no Session, no variables, no printing. Every panel contributes a block
of equations (rows): regressors X_i, dependent y_i, instruments Z_i (columns shared by all the
panels) and the initial weight matrix H_i of its errors."""

from __future__ import annotations

import numpy as np


def tridiagonal(times: np.ndarray) -> np.ndarray:
    """H of first-differenced errors: 2 on the diagonal, -1 between equations of consecutive periods."""
    n = len(times)
    H = 2.0 * np.eye(n)
    for a in range(n - 1):
        if times[a + 1] - times[a] == 1:
            H[a, a + 1] = H[a + 1, a] = -1.0
    return H


def _pinv(A: np.ndarray) -> np.ndarray:
    return np.linalg.pinv((A + A.T) / 2, rcond=1e-10)


def gmm_fit(X, y, Z, H, twostep: bool = False, robust: bool = False, nd=None) -> dict:
    """X, y, Z, H: lists with one block per panel. One-step estimator with W = (sum Z'HZ)^-1;
    the two-step estimator weights with the inverse of sum Z'e1 e1'Z. The covariance matrix of the
    two-step estimator is the asymptotic one, or with `robust` the Windmeijer (2005) correction."""
    G = len(X)
    k, L = X[0].shape[1], Z[0].shape[1]
    XZ = sum(x.T @ z for x, z in zip(X, Z))                    # k x L
    Zy = sum(z.T @ yy for z, yy in zip(Z, y))
    ZHZ = sum(z.T @ h @ z for z, h in zip(Z, H))
    W1 = _pinv(ZHZ)
    A1 = XZ @ W1 @ XZ.T
    b1 = np.linalg.solve(A1, XZ @ W1 @ Zy)
    e1 = [yy - x @ b1 for yy, x in zip(y, X)]
    n = sum(len(e) for e in e1)
    # sigma^2 of the idiosyncratic differenced error under Var = sigma^2 H   (VERIFY)
    sig2 = sum(float(e @ np.linalg.solve(h, e)) for e, h in zip(e1, H)) / (n - k)
    if nd is not None:
        # VERIFY: Stata's sigma^2 is the sum of squared differenced residuals over 2 (n_d - k), as found
        # in the 1905 log
        sig2 = sum(float(e[:m] @ e[:m]) for e, m in zip(e1, nd)) / (2 * (sum(nd) - k))
    Om1 = sum(np.outer(z.T @ e, z.T @ e) for z, e in zip(Z, e1))
    A1i = np.linalg.inv(A1)
    V1 = sig2 * A1i
    V1r = A1i @ (XZ @ W1 @ Om1 @ W1 @ XZ.T) @ A1i
    out = {"b1": b1, "e1": e1, "V1": V1, "V1r": V1r, "W1": W1, "sigma2": sig2, "XZ": XZ, "ZHZ": ZHZ,
           "n": n, "k": k, "L": L, "G": G}
    if not twostep:
        out.update(b=b1, e=e1, V=V1r if robust else V1, W=W1, A=A1)
        return out
    W2 = _pinv(Om1)
    A2 = XZ @ W2 @ XZ.T
    A2i = np.linalg.inv(A2)
    b2 = A2i @ (XZ @ W2 @ Zy)
    e2 = [yy - x @ b2 for yy, x in zip(y, X)]
    V2 = A2i
    if robust:
        Ze2 = sum(z.T @ e for z, e in zip(Z, e2))
        D = np.zeros((k, k))
        for j in range(k):
            dOm = sum(np.outer(z.T @ (x[:, j] * 1.0), z.T @ e) + np.outer(z.T @ e, z.T @ x[:, j])
                      for z, x, e in zip(Z, X, e1))
            D[:, j] = A2i @ XZ @ W2 @ dOm @ W2 @ Ze2
        Vc = V2 + D @ V2 + V2 @ D.T + D @ V1r @ D.T
        V = (Vc + Vc.T) / 2
    else:
        V = V2
    out.update(b=b2, e=e2, V=V, W=W2, A=A2, V2=V2)
    return out


def overid_test(fit: dict, Z, twostep: bool) -> tuple[float, int]:
    """Sargan statistic (one-step, homoskedastic errors) or Hansen J (two-step); degrees of
    freedom L - k."""
    e = fit["e"]
    Ze = sum(z.T @ ee for z, ee in zip(Z, e))
    if twostep:
        stat = float(Ze @ fit["W"] @ Ze)
    else:
        stat = float(Ze @ fit["W1"] @ Ze) / fit["sigma2"]
    return stat, int(np.linalg.matrix_rank(fit["ZHZ"]) - fit["k"])


def ar_test(order: int, fit: dict, X, Z, times) -> tuple[float, float]:
    """Arellano-Bond test of serial correlation of order `order` in the first-differenced
    residuals (equations of the differenced model only: `times[i]` holds the period of each row of
    panel i, or -1 for the level equations). Returns (z, d)."""
    e = fit["e"]
    d = 0.0
    t1 = 0.0
    t2 = np.zeros(fit["XZ"].shape[1])
    wX = np.zeros(fit["k"])
    for ei, xi, zi, ti in zip(e, X, Z, times):
        diff = ti >= 0
        pos = {int(t): j for j, t in enumerate(ti) if t >= 0}
        w = np.zeros(len(ei))
        for j in np.flatnonzero(diff):
            jj = pos.get(int(ti[j]) - order)
            if jj is not None:
                w[j] = ei[jj]
        ew = float(w @ ei)
        d += ew
        t1 += ew ** 2
        t2 += zi.T @ ei * ew
        wX += w @ xi
    W = fit["W"]
    A = fit["A"]
    mid = wX @ np.linalg.solve(A, fit["XZ"] @ W @ t2)
    var = t1 - 2 * mid + wX @ fit["V"] @ wX
    if var <= 0:
        return float("nan"), d
    return d / float(np.sqrt(var)), d

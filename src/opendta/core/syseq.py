"""Linear regression systems and variants, on arrays: restricted least squares (cnsreg), seemingly
unrelated regression and three-stage least squares (sureg, reg3), variance-weighted least squares
(vwls), errors-in-variables regression (eivreg) and robust regression (rreg). No Session and no
output: the command module in estimation/sysmodels.py reads the data and prints."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


# ---------------------------------------------------------------------------
# cnsreg: least squares under R b = r
# ---------------------------------------------------------------------------

def restricted_ls(X: np.ndarray, y: np.ndarray, w: np.ndarray, R: np.ndarray, r: np.ndarray):
    """(b, A, rank R) where A = (X'WX)^-1 projected on R b = r: V = s2 * A with s2 = RSS / (N - k + rank)."""
    XtWX = (X * w[:, None]).T @ X
    inv = np.linalg.pinv(XtWX)
    b0 = inv @ ((X * w[:, None]).T @ y)
    if R.shape[0] == 0:
        return b0, inv, 0
    M = R @ inv @ R.T
    adj = inv @ R.T @ np.linalg.pinv(M)
    b = b0 + adj @ (r - R @ b0)
    A = inv - adj @ R @ inv
    return b, (A + A.T) / 2, int(np.linalg.matrix_rank(R))


# ---------------------------------------------------------------------------
# sureg and reg3
# ---------------------------------------------------------------------------

@dataclass
class SystemFit:
    b: list                          # coefficients of each equation
    V: np.ndarray                    # covariance of the stacked coefficients
    sigma: np.ndarray                # residual covariance used in the last step
    resid: np.ndarray                # n x m residuals with the final coefficients
    iterations: int = 0
    converged: bool = True
    chi2_bp: float = 0.0             # Breusch-Pagan test of independence
    df_bp: int = 0
    corr: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    trace: list = field(default_factory=list)   # change in the coefficients at each iteration after the first


def _sigma(E: np.ndarray, ks: list, dfk: str) -> np.ndarray:
    """Residual covariance: divisor n, or sqrt((n - k_i)(n - k_j)) with dfk = "small" (every
    equation) or n - max(k_i, k_j) with dfk = "dfk"."""
    n, m = E.shape
    S = E.T @ E
    D = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            if dfk == "small":
                D[i, j] = np.sqrt((n - ks[i]) * (n - ks[j]))
            elif dfk == "dfk":
                D[i, j] = n - max(ks[i], ks[j])
            else:
                D[i, j] = n
    return S / D


def _gls(Xs: list, Zs: list, ys: list, sigma: np.ndarray):
    """Stacked GLS with instruments Zs for the regressors (Zs = Xs for SUR)."""
    m = len(Xs)
    ks = [X.shape[1] for X in Xs]
    offs = np.concatenate([[0], np.cumsum(ks)]).astype(int)
    Si = np.linalg.inv(sigma)
    A = np.zeros((offs[-1], offs[-1]))
    c = np.zeros(offs[-1])
    for i in range(m):
        for j in range(m):
            A[offs[i]:offs[i + 1], offs[j]:offs[j + 1]] = Si[i, j] * (Zs[i].T @ Xs[j])
            c[offs[i]:offs[i + 1]] += Si[i, j] * (Zs[i].T @ ys[j])
    V = np.linalg.inv(A)
    b = V @ c
    return [b[offs[i]:offs[i + 1]] for i in range(m)], V


def _resid(Xs: list, ys: list, bs: list) -> np.ndarray:
    return np.column_stack([y - X @ b for X, y, b in zip(Xs, ys, bs)])


def sure(Xs: list, ys: list, *, iterate: bool = False, dfk: str = "", maxit: int = 100,
         tol: float = 1e-6) -> SystemFit:
    """Feasible GLS for seemingly unrelated regressions. The first sigma comes from the equation-by-
    equation OLS residuals; with `iterate` it is updated until the coefficients stop changing."""
    ks = [X.shape[1] for X in Xs]
    bs = [np.linalg.lstsq(X, y, rcond=None)[0] for X, y in zip(Xs, ys)]
    E = _resid(Xs, ys, bs)
    sigma = _sigma(E, ks, dfk)
    it, conv, trace = 0, True, []
    while True:
        bs_new, V = _gls(Xs, Xs, ys, sigma)
        E = _resid(Xs, ys, bs_new)
        it += 1
        if not iterate:
            bs = bs_new
            break
        if it > 1:                       # the first GLS step is not an iteration (Stata's log starts at its change)
            delta = max(float(np.max(np.abs(a - c)) / max(1.0, np.max(np.abs(c)))) for a, c in zip(bs_new, bs))
            trace.append(delta)
            if delta < tol:
                bs = bs_new
                break
        bs = bs_new
        if it >= maxit:
            conv = False
            break
        sigma = _sigma(E, ks, dfk)
    return _finish_system(bs, V, E, sigma, it, conv, ks, dfk, trace)


def three_sls(Xs: list, Zs: list, ys: list, *, method: str = "3sls", iterate: bool = False,
              dfk: str = "", maxit: int = 100, tol: float = 1e-6) -> SystemFit:
    """Reg3 estimators with the regressors Xs and the first-stage projections Zs (instruments
    projected on every equation's regressors, i.e. X-hat). `method`: ols, 2sls, sure, 3sls."""
    ks = [X.shape[1] for X in Xs]
    m = len(Xs)
    if method in ("ols", "2sls"):
        src = Xs if method == "ols" else Zs
        bs = [np.linalg.solve(Z.T @ X, Z.T @ y) for Z, X, y in zip(src, Xs, ys)]
        E = _resid(Xs, ys, bs)
        sigma = _sigma(E, ks, dfk)
        V = np.zeros((sum(ks), sum(ks)))
        o = 0
        for i in range(m):
            Zi = src[i]
            Ai = np.linalg.inv(Zi.T @ Xs[i])
            V[o:o + ks[i], o:o + ks[i]] = sigma[i, i] * (Ai @ (Zi.T @ Zi) @ Ai.T)
            o += ks[i]
        return _finish_system(bs, V, E, sigma, 1, True, ks, dfk)
    first = "ols" if method == "sure" else "2sls"
    start = three_sls(Xs, Zs, ys, method=first, dfk=dfk)
    sigma, bs = start.sigma, start.b
    src = Xs if method == "sure" else Zs
    it, conv = 0, True
    while True:
        bs_new, V = _gls(Xs, src, ys, sigma)
        E = _resid(Xs, ys, bs_new)
        it += 1
        if not iterate:
            bs = bs_new
            break
        delta = max(float(np.max(np.abs(a - c)) / max(1.0, np.max(np.abs(c)))) for a, c in zip(bs_new, bs))
        bs = bs_new
        if delta < tol:
            break
        if it >= maxit:
            conv = False
            break
        sigma = _sigma(E, ks, dfk)
    return _finish_system(bs, V, E, sigma, it, conv, ks, dfk)


def _finish_system(bs, V, E, sigma, it, conv, ks, dfk, trace=()) -> SystemFit:
    n, m = E.shape
    S = E.T @ E / n
    d = np.sqrt(np.diag(S))
    corr = S / np.outer(d, d)
    bp = n * float(np.sum(np.tril(corr, -1) ** 2))
    return SystemFit(bs, V, sigma, E, it, conv, bp, m * (m - 1) // 2, corr, list(trace))


# ---------------------------------------------------------------------------
# vwls
# ---------------------------------------------------------------------------

def vwls(X: np.ndarray, y: np.ndarray, sd: np.ndarray):
    """Weighted least squares with weights 1/sd^2 where the variances are known: (b, V, chi2 of fit)."""
    w = 1.0 / sd ** 2
    V = np.linalg.inv((X * w[:, None]).T @ X)
    b = V @ ((X * w[:, None]).T @ y)
    chi2 = float(np.sum(w * (y - X @ b) ** 2))
    return b, V, chi2


# ---------------------------------------------------------------------------
# eivreg
# ---------------------------------------------------------------------------

def eivreg(X: np.ndarray, y: np.ndarray, reliab: np.ndarray, constant: bool = True):
    """Errors-in-variables regression (Fuller 1987, 2.2): the moments of the regressors are
    reduced by the measurement error implied by their reliabilities. X holds the slopes only (no
    constant column). Returns (slopes, intercept, V of [slopes, intercept], sigma_ee)."""
    n, p = X.shape
    if constant:
        xm, ym = X.mean(axis=0), float(y.mean())
        Xc, yc = X - xm, y - ym
    else:
        xm, ym = np.zeros(p), 0.0
        Xc, yc = X, y
    den = (n - 1) if constant else n
    Mxx = Xc.T @ Xc / den
    myy = float(yc @ yc / den)
    mxy = Xc.T @ yc / den
    Svv = np.diag((1.0 - reliab) * np.diag(Mxx))
    Mc = Mxx - Svv
    Mi = np.linalg.inv(Mc)
    beta = Mi @ mxy
    see = (den / (n - p - (1 if constant else 0))) * max(myy - float(mxy @ beta), 0.0) \
        if n - p - (1 if constant else 0) > 0 else 0.0
    u = see + float(beta @ Svv @ beta)
    # VERIFY: variance of Fuller (1987, 2.2.14) as the manual gives it
    Vb = (u * Mi + Mi @ np.outer(Svv @ beta, Svv @ beta) @ Mi) / den
    if not constant:
        return beta, 0.0, Vb, see
    a = ym - float(xm @ beta)
    var_a = 1.0 / n * u + float(xm @ Vb @ xm)
    cov = -Vb @ xm
    V = np.zeros((p + 1, p + 1))
    V[:p, :p] = Vb
    V[:p, p] = cov
    V[p, :p] = cov
    V[p, p] = var_a
    return beta, a, V, see


# ---------------------------------------------------------------------------
# rreg
# ---------------------------------------------------------------------------

def cooks_distance(X: np.ndarray, e: np.ndarray, s2: float) -> np.ndarray:
    XtXi = np.linalg.pinv(X.T @ X)
    h = np.einsum("ij,jk,ik->i", X, XtXi, X)
    return e ** 2 * h / (X.shape[1] * s2 * (1 - h) ** 2)


def _wls(X, y, w):
    Xw = X * w[:, None]
    return np.linalg.pinv(Xw.T @ X) @ (Xw.T @ y)


def rreg(X: np.ndarray, y: np.ndarray, *, tune: float = 7.0, tol: float = 0.01, maxit: int = 100,
         huber_c: float = 1.339, drop_cook: bool = True):
    """Robust regression by iteratively reweighted least squares ([R] rreg): observations with
    Cook's D > 1 are dropped, then Huber weights iterate until the largest change in the weights is
    below `tol`, then biweight weights. X includes the constant column. Returns (b, V, weights,
    keep, log, scale) where log = [("Huber"|"Biweight", iteration, max difference)]."""
    n, k = X.shape
    keep = np.ones(n, dtype=bool)
    b = np.linalg.pinv(X.T @ X) @ (X.T @ y)
    e = y - X @ b
    if drop_cook:
        s2 = float(e @ e) / (n - k)
        keep = cooks_distance(X, e, s2) <= 1.0
    Xk, yk = X[keep], y[keep]
    w = np.ones(keep.sum())
    log = []
    b = _wls(Xk, yk, w)

    def scale(res):
        return float(np.median(np.abs(res - np.median(res)))) / 0.6745 or 1e-300

    for kind in ("Huber", "Biweight"):
        c = huber_c if kind == "Huber" else tune
        for it in range(1, maxit + 1):
            res = yk - Xk @ b
            u = res / scale(res)
            if kind == "Huber":
                wn = np.where(np.abs(u) <= c, 1.0, c / np.maximum(np.abs(u), 1e-300))
            else:
                wn = np.where(np.abs(u) < c, (1 - (u / c) ** 2) ** 2, 0.0)
            diff = float(np.max(np.abs(wn - w)))
            w = wn
            b = _wls(Xk, yk, w)
            log.append((kind, it, diff))
            if diff < tol:
                break
    res = yk - Xk @ b
    s = scale(res)
    u = res / s
    c = tune
    inside = np.abs(u) < c
    psi = np.where(inside, u * (1 - (u / c) ** 2) ** 2, 0.0)
    dpsi = np.where(inside, (1 - (u / c) ** 2) * (1 - 5 * (u / c) ** 2), 0.0)
    nk, p = Xk.shape
    m = float(dpsi.mean())
    kappa = 1 + (p / nk) * float(np.var(dpsi)) / (m * m) if m else 1.0
    sig2 = float(np.sum(psi ** 2)) / (nk - p)
    V = kappa ** 2 * sig2 / (m * m if m else 1.0) * np.linalg.pinv(Xk.T @ Xk)
    full_w = np.zeros(n)
    full_w[keep] = w
    return b, (V + V.T) / 2, full_w, keep, log, s

"""Panel unit-root tests: Harris-Tzavalis, Breitung, Hadri and Fisher-type tests
([XT] xtunitroot, Methods and formulas).

Plain functions on arrays: no Session, no variables, no printing. A balanced panel is a matrix
Y of N rows (panels) by T columns (consecutive periods)."""

from __future__ import annotations

import math

import numpy as np

from ..lazy import lazy

st = lazy("scipy.stats")


def adf_tau(y: np.ndarray, lags: int = 0, const: bool = True, trend: bool = False) -> tuple[float, int]:
    """t statistic of rho in dy_t = rho y_{t-1} + sum_j g_j dy_{t-j} (+ constant, + trend), and
    the number of observations used."""
    dy = np.diff(y)
    n = len(dy) - lags
    if n <= lags + 3:
        raise ValueError("insufficient observations")
    cols = [y[lags:-1]] + [dy[lags - j:len(dy) - j] for j in range(1, lags + 1)]
    if const:
        cols.append(np.ones(n))
    if trend:
        cols.append(np.arange(1, n + 1, dtype=float))
    X = np.column_stack(cols)
    yy = dy[lags:]
    b = np.linalg.lstsq(X, yy, rcond=None)[0]
    e = yy - X @ b
    s2 = float(e @ e) / (n - X.shape[1])
    V = s2 * np.linalg.inv(X.T @ X)
    return float(b[0] / math.sqrt(V[0, 0])), n


def harris_tzavalis(Y: np.ndarray, const: bool = True) -> dict:
    """rho is the least-squares estimate of the autoregressive parameter of the pooled panels
    (with panel means when const); z = sqrt(N)(rho - mean)/sd, left tail."""
    N, T = Y.shape
    x, y = Y[:, :-1], Y[:, 1:]
    if const:
        x = x - x.mean(axis=1, keepdims=True)
        y = y - y.mean(axis=1, keepdims=True)
    rho = float((x * y).sum() / (x * x).sum())
    if const:
        mean = 1 - 3 / (T + 1)
        var = 3 * (17 * T ** 2 - 20 * T + 17) / (5 * (T - 1) * (T + 1) ** 3)
    else:
        mean = 1.0
        var = 2 / (T * (T - 1))
    z = math.sqrt(N) * (rho - mean) / math.sqrt(var)
    return {"rho": rho, "z": z, "p": float(st.norm.cdf(z))}


def breitung(Y: np.ndarray, trend: bool = False) -> dict:
    """Breitung's lambda (2000), without prewhitening. The statistic is asymptotically N(0, 1),
    left tail."""
    N, T = Y.shape
    dY = np.diff(Y, axis=1)                                   # dy_2 ... dy_T
    # forward orthogonalization of dy_t, t = 2..T-1
    n_t = T - 1
    fo = np.empty((N, n_t - 1))
    for j in range(n_t - 1):
        rest = dY[:, j + 1:].mean(axis=1)
        fo[:, j] = math.sqrt((n_t - 1 - j) / (n_t - j)) * (dY[:, j] - rest)
    tt = np.arange(1, T + 1)
    if trend:
        lagv = Y[:, :-1] - Y[:, [0]] - ((tt[:-1] - 1) / (T - 1)) * (Y[:, [-1]] - Y[:, [0]])
    else:
        lagv = Y[:, :-1] - Y[:, [0]]
    ylag = lagv[:, :n_t - 1]                                   # y*_{t-1} for t = 2..T-1 aligned to fo
    s2 = float((fo ** 2).sum()) / (N * (T - 2))
    lam = float((fo * ylag).sum() / math.sqrt(s2 * (ylag ** 2).sum()))
    return {"lambda": lam, "p": float(st.norm.cdf(lam))}


def lr_variance(e: np.ndarray, kernel: str, lags: int) -> float:
    """Long-run variance of the residuals e: sum of the autocovariances (divisor T) weighted by
    the Bartlett, Parzen or quadratic-spectral kernel with bandwidth lags + 1."""
    T = len(e)
    g = [float(e[j:] @ e[:T - j]) / T for j in range(T)]
    out = g[0]
    b = lags + 1
    for j in range(1, T):
        x = j / b
        if kernel == "bartlett":
            w = 1 - x if x < 1 else 0.0
        elif kernel == "parzen":
            w = 1 - 6 * x * x + 6 * x ** 3 if x <= 0.5 else (2 * (1 - x) ** 3 if x <= 1 else 0.0)
        else:                                                   # quadratic spectral
            d = 6 * math.pi * x / 5
            w = 25 / (12 * math.pi ** 2 * x * x) * (math.sin(d) / d - math.cos(d))
        if kernel != "quadraticspectral" and w == 0.0:
            break
        out += 2 * w * g[j]
    return out


def hadri(Y: np.ndarray, trend: bool = False, robust: bool = False, kernel: str = "", lags: int = 0) -> dict:
    """Hadri (2000) LM test of stationarity; z is asymptotically N(0, 1), right tail. kernel: the
    variance of the residuals is their long-run variance (per panel with robust, else the mean
    over the panels; VERIFY)."""
    N, T = Y.shape
    tt = np.arange(1, T + 1, dtype=float)
    Z = np.column_stack([np.ones(T), tt]) if trend else np.ones((T, 1))
    P = Z @ np.linalg.pinv(Z)
    E = Y - Y @ P.T
    S = np.cumsum(E, axis=1)
    if kernel:
        lrv = np.array([lr_variance(E[i], kernel, lags) for i in range(N)])
        if robust:
            lm = float(((S ** 2).sum(axis=1) / (T ** 2 * lrv)).mean())
        else:
            lm = float((S ** 2).sum() / (N * T ** 2 * lrv.mean()))
    elif robust:
        sig2 = (E ** 2).mean(axis=1)
        lm = float(((S ** 2).sum(axis=1) / (T ** 2 * sig2)).mean())
    else:
        # VERIFY: the variance is divided by N (T - k), k = number of deterministic terms (1905 log)
        sig2 = float((E ** 2).sum()) / (N * (T - Z.shape[1]))
        lm = float((S ** 2).sum() / (N * T ** 2 * sig2))
    xi, zeta2 = (1 / 15, 11 / 6300) if trend else (1 / 6, 1 / 45)
    z = math.sqrt(N) * (lm - xi) / math.sqrt(zeta2)
    return {"lm": lm, "z": z, "p": float(st.norm.sf(z))}


def fisher_combine(p: np.ndarray) -> dict:
    """Maddala-Wu and Choi combinations of the panel p-values. Returns each statistic and its
    p-value, with the degrees of freedom of the chi2 and t references."""
    p = np.clip(np.asarray(p, dtype=float), 1e-300, 1 - 1e-16)
    N = len(p)
    P = float(-2 * np.log(p).sum())
    Z = float(st.norm.ppf(p).sum() / math.sqrt(N))
    k = 3 * (5 * N + 4) / (math.pi ** 2 * N * (5 * N + 2))
    L = float(math.sqrt(k) * np.log(p / (1 - p)).sum())
    Pm = (P - 2 * N) / (2 * math.sqrt(N))
    return {"P": P, "df_P": 2 * N, "p_P": float(st.chi2.sf(P, 2 * N)),
            "Z": Z, "p_Z": float(st.norm.cdf(Z)),
            "L": L, "df_L": 5 * N + 4, "p_L": float(st.t.cdf(L, 5 * N + 4)),
            "Pm": Pm, "p_Pm": float(st.norm.sf(Pm))}

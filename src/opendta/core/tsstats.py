"""Time-series filters and smoothers: plain functions on arrays ([TS] tsfilter, [TS] tssmooth).

Nothing here reads variables, uses a Session or prints; opendta.commands.tsfilters does that.
Missing values are NaN in the arguments and in the results.

tsfilter: Hodrick-Prescott, Baxter-King, Christiano-Fitzgerald and Butterworth filters, each
returning the cyclical component (the trend is y minus the cycle) from the formulas of the
manual. tssmooth: exponential, double exponential and Holt-Winters recursions with the starting
values and the penalized criterion described in the manual.
"""

from __future__ import annotations

import math

import numpy as np

from ..lazy import lazy

sparse = lazy("scipy.sparse")
spla = lazy("scipy.sparse.linalg")
optimize = lazy("scipy.optimize")


# ---------------------------------------------------------------------------
# tsfilter
# ---------------------------------------------------------------------------

def hp_cycle(y: np.ndarray, lam: float) -> np.ndarray:
    """Hodrick-Prescott cycle c = (I - M^-1) y, M = I + lam K'K with K the second-difference
    matrix; for T < 3 the cycle is zero."""
    T = len(y)
    if T < 3:
        return np.zeros(T)
    K = sparse.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(T - 2, T), format="csc")
    M = (sparse.identity(T, format="csc") + lam * (K.T @ K)).tocsc()
    return y - spla.spsolve(M, y)


def hp_gain(lam: float, omega: np.ndarray) -> np.ndarray:
    a = 4 * lam * (1 - np.cos(omega)) ** 2
    return a / (1 + a)


def bk_weights(minperiod: float, maxperiod: float, q: int, stationary: bool = False) -> np.ndarray:
    """Weights b~_0..b~_q of the Baxter-King filter: the ideal band-pass weights cut at q and,
    unless `stationary`, shifted by their mean so that they sum to zero."""
    wl, wh = 2 * math.pi / maxperiod, 2 * math.pi / minperiod
    j = np.arange(0, q + 1, dtype=float)
    b = np.empty(q + 1)
    b[0] = (wh - wl) / math.pi
    b[1:] = (np.sin(j[1:] * wh) - np.sin(j[1:] * wl)) / (j[1:] * math.pi)
    if not stationary:
        total = b[0] + 2 * b[1:].sum()
        b = b - total / (2 * q + 1)
    return b


def bk_cycle(y: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Symmetric moving average with the weights w[0..q]; the first and last q values are NaN."""
    q = len(w) - 1
    T = len(y)
    out = np.full(T, np.nan)
    if T <= 2 * q:
        return out
    full = np.concatenate([w[:0:-1], w])
    out[q:T - q] = np.convolve(y, full, mode="valid")
    return out


def bk_gain(w: np.ndarray, omega: np.ndarray) -> np.ndarray:
    j = np.arange(1, len(w))
    return np.abs(w[0] + 2 * (w[1:][None, :] * np.cos(np.outer(omega, j))).sum(axis=1))


def cf_ideal(minperiod: float, maxperiod: float, n: int) -> np.ndarray:
    """Ideal band-pass weights b_0..b_n."""
    wl, wh = 2 * math.pi / maxperiod, 2 * math.pi / minperiod
    j = np.arange(0, n + 1, dtype=float)
    b = np.empty(n + 1)
    b[0] = (wh - wl) / math.pi
    b[1:] = (np.sin(j[1:] * wh) - np.sin(j[1:] * wl)) / (j[1:] * math.pi)
    return b


def cf_cycle(y: np.ndarray, minperiod: float, maxperiod: float, *, stationary: bool = False,
             drift: bool = False, smaorder: int | None = None) -> np.ndarray:
    """Christiano-Fitzgerald cycle. Without smaorder the filter is asymmetric and uses the
    whole sample (end weights that sum the weights to zero, unless `stationary`); with
    smaorder q it is the symmetric moving average of order q (first and last q values NaN)."""
    T = len(y)
    z = np.asarray(y, dtype=float)
    if drift and T > 1:
        z = z - np.arange(T) * (z[-1] - z[0]) / (T - 1)
    if smaorder is not None:
        q = int(smaorder)
        b = cf_ideal(minperiod, maxperiod, q)
        out = np.full(T, np.nan)
        if T <= 2 * q:
            return out
        bq = b[q] if stationary else -b[0] / 2 - b[1:q].sum()
        w = np.concatenate([[bq], b[1:q][::-1], [b[0]], b[1:q], [bq]])
        out[q:T - q] = np.convolve(z, w[::-1], mode="valid")
        return out
    b = cf_ideal(minperiod, maxperiod, max(T - 1, 1))
    cs = np.concatenate([[0.0], np.cumsum(b[1:])])      # cs[k] = b_1 + ... + b_k

    def end(k: int) -> float:                           # ~b_k
        return b[k] if stationary else -b[0] / 2 - cs[max(k - 1, 0)]
    out = np.empty(T)
    for t in range(1, T + 1):                           # t is the 1-based period
        i = t - 1
        acc = b[0] * z[i]
        nlead, nlag = T - t - 1, t - 2
        if nlead > 0:
            acc += b[1:nlead + 1] @ z[i + 1:i + nlead + 1]
        if nlag > 0:
            acc += b[1:nlag + 1] @ z[i - nlag:i][::-1]
        if t < T:
            acc += end(T - t) * z[T - 1]
        if t > 1:
            acc += end(t - 1) * z[0]
        out[i] = acc
    return out


def bw_cycle(y: np.ndarray, maxperiod: float, order: int) -> np.ndarray:
    """Butterworth cycle c = lam Q (|Q'Q| + lam Q'Q)^-1 Q'y with Q' the (T-m) x T matrix of the
    coefficients of (1 - L)^m and lam = tan(pi/maxperiod)^(-2m)."""
    T = len(y)
    m = int(order)
    if T <= m + 1:
        return np.zeros(T)
    coef = np.array([(-1) ** k * math.comb(m, k) for k in range(m + 1)], dtype=float)
    Qt = np.zeros((T - m, T))
    for i in range(T - m):
        Qt[i, i:i + m + 1] = coef
    lam = math.tan(math.pi / maxperiod) ** (-2 * m)
    H = Qt @ Qt.T
    L = np.abs(H)
    d = np.linalg.solve(L + lam * H, Qt @ y)
    return lam * (Qt.T @ d)


def bw_gain(maxperiod: float, order: int, omega: np.ndarray) -> np.ndarray:
    wc = 2 * math.pi / maxperiod
    with np.errstate(divide="ignore", invalid="ignore"):
        return 1.0 / (1.0 + (math.tan(wc / 2) / np.tan(omega / 2)) ** (2 * order))


# ---------------------------------------------------------------------------
# tssmooth
# ---------------------------------------------------------------------------

def _first_half(n_obs: int) -> int:
    return max(n_obs // 2, 1)


def linear_start(x: np.ndarray, n0: int) -> tuple[float, float]:
    """Intercept and slope of the regression of x on t = 1..n0 (the first n0 periods, ignoring
    missing values): the level at time 0 and the trend."""
    t = np.arange(1, n0 + 1, dtype=float)
    v = x[:n0]
    ok = ~np.isnan(v)
    if ok.sum() < 2:
        return (float(v[ok][0]) if ok.any() else math.nan), 0.0
    b1, b0 = np.polyfit(t[ok], v[ok], 1)
    return float(b0), float(b1)


def exp_run(x: np.ndarray, alpha: float, s0: float):
    """Single exponential recursion. Returns pred (pred[t] = S_{t-1}, the forecast of x_t), S and
    the sum of squared forecast errors; a missing x_t leaves S unchanged."""
    n = len(x)
    pred = np.empty(n)
    S = np.empty(n)
    s = s0
    rss = 0.0
    for t in range(n):
        pred[t] = s
        if not math.isnan(x[t]):
            rss += (x[t] - s) ** 2
            s = alpha * x[t] + (1 - alpha) * s
        S[t] = s
    return pred, S, rss


def exp_smooth(x: np.ndarray, *, alpha: float | None = None, s0: float | None = None,
               samp0: int | None = None, forecast: int = 0) -> dict:
    n = len(x)
    n_obs = int((~np.isnan(x)).sum())
    n0 = samp0 if samp0 else _first_half(n_obs)
    if s0 is None:
        v = x[:n0]
        s0 = float(np.nanmean(v)) if (~np.isnan(v)).any() else math.nan
    if alpha is None:
        res = optimize.minimize_scalar(lambda a: exp_run(x, a, s0)[2], bounds=(1e-6, 1 - 1e-6),
                                       method="bounded", options={"xatol": 1e-10})
        alpha = float(res.x)
    pred, S, rss = exp_run(x, alpha, s0)
    out = np.concatenate([pred, np.full(forecast, S[-1] if n else math.nan)])
    return {"fit": out, "alpha": alpha, "rss": rss, "rmse": math.sqrt(rss / n_obs) if n_obs else math.nan,
            "N": n_obs, "N_pre": n0, "s1_0": s0, "S": S}


def dexp_run(x: np.ndarray, alpha: float, b0: float, b1: float):
    c = (1 - alpha) / alpha
    s1 = b0 - c * b1
    s2 = b0 - 2 * c * b1
    k = alpha / (1 - alpha)
    n = len(x)
    pred = np.empty(n)
    rss = 0.0
    for t in range(n):
        pred[t] = (2 + k) * s1 - (1 + k) * s2
        if not math.isnan(x[t]):
            rss += (x[t] - pred[t]) ** 2
            s1 = alpha * x[t] + (1 - alpha) * s1
            s2 = alpha * s1 + (1 - alpha) * s2
        else:
            s1 = alpha * pred[t] + (1 - alpha) * s1
            s2 = alpha * s1 + (1 - alpha) * s2
    return pred, rss, s1, s2


def dexp_smooth(x: np.ndarray, *, alpha: float | None = None, s0: tuple | None = None,
                samp0: int | None = None, forecast: int = 0) -> dict:
    n = len(x)
    n_obs = int((~np.isnan(x)).sum())
    n0 = samp0 if samp0 else _first_half(n_obs)
    if s0 is None:
        b0, b1 = linear_start(x, n0)
    else:
        # s0(S1 S2): translate the two starting smoothed values into the regression form
        b0, b1 = None, None
    if alpha is None:
        res = optimize.minimize_scalar(lambda a: dexp_run(x, a, b0, b1)[1], bounds=(1e-6, 1 - 1e-6),
                                       method="bounded", options={"xatol": 1e-10})
        alpha = float(res.x)
    pred, rss, s1, s2 = dexp_run(x, alpha, b0, b1)
    k = alpha / (1 - alpha)
    a_T, b_T = 2 * s1 - s2, k * (s1 - s2)
    fc = a_T + b_T * np.arange(1, forecast + 1)
    return {"fit": np.concatenate([pred, fc]), "alpha": alpha, "rss": rss,
            "rmse": math.sqrt(rss / n_obs) if n_obs else math.nan, "N": n_obs, "N_pre": n0,
            "s1_0": b0 - (1 - alpha) / alpha * b1, "s2_0": b0 - 2 * (1 - alpha) / alpha * b1,
            "constant": a_T, "linear": b_T}


def _logit_penalty(p: np.ndarray) -> float:
    """Penalty of the Holt-Winters criterion: zero unless a parameter is within e^-12 of 0 or 1
    (|ln(x(1-x))| > 12)."""
    pen = 0.0
    for v in p:
        v = min(max(v, 1e-300), 1 - 1e-16)
        f = abs(math.log(v * (1 - v)))
        if f > 12:
            pen += (f - 12) ** 2
    return pen


def hw_run(x: np.ndarray, alpha: float, beta: float, a0: float, b0: float):
    n = len(x)
    pred = np.empty(n)
    a, b = a0, b0
    rss = 0.0
    for t in range(n):
        pred[t] = a + b
        xt = x[t] if not math.isnan(x[t]) else pred[t]
        if not math.isnan(x[t]):
            rss += (x[t] - pred[t]) ** 2
        an = alpha * xt + (1 - alpha) * (a + b)
        b = beta * (an - a) + (1 - beta) * b
        a = an
    return pred, rss, a, b


def hw_start(x: np.ndarray, n0: int, diff: bool) -> tuple[float, float]:
    if diff:
        d = np.diff(x[:max(n0, 2)])
        d = d[~np.isnan(d)]
        b0 = float(d.mean()) if d.size else 0.0
        return float(x[0]) - b0, b0
    return linear_start(x, n0)


def hw_smooth(x: np.ndarray, *, parms=None, s0=None, samp0=None, forecast=0, diff=False,
              start=(0.5, 0.5), maximize=None) -> dict:
    """Holt-Winters nonseasonal. `maximize(fun, b0)` runs the optimizer on the logit scale and
    returns an MLResult (the command supplies it, to print the log)."""
    n_obs = int((~np.isnan(x)).sum())
    n0 = samp0 if samp0 else _first_half(n_obs)
    a0, b0 = s0 if s0 is not None else hw_start(x, n0, diff)
    res = None
    if parms is None:
        def fun(z, todo):
            al, be = 1 / (1 + np.exp(-z[0])), 1 / (1 + np.exp(-z[1]))
            return -(hw_run(x, al, be, a0, b0)[1] + _logit_penalty((al, be))), None, None
        res = maximize(fun, np.log(np.array(start) / (1 - np.array(start))))
        alpha, beta = [float(1 / (1 + math.exp(-v))) for v in res.b]
    else:
        alpha, beta = parms
    pred, rss, a, b = hw_run(x, alpha, beta, a0, b0)
    fc = a + b * np.arange(1, forecast + 1)
    prss = rss + _logit_penalty((alpha, beta))
    return {"fit": np.concatenate([pred, fc]), "alpha": alpha, "beta": beta, "rss": rss,
            "prss": prss, "rmse": math.sqrt(rss / n_obs) if n_obs else math.nan, "N": n_obs,
            "N_pre": n0, "s1_0": a0, "s2_0": b0, "constant": a, "linear": b, "opt": res}


def shw_start(x: np.ndarray, L: int, n0: int, additive: bool, normalize: bool, alt: bool):
    """Starting level, trend and seasonal terms. Default: from the first m complete years of
    the first half of the data (m years averaged), as in the manual; `alt`: regressions."""
    n = len(x)
    m = max(n0 // L, 1)
    if alt:
        t = np.arange(1, n0 + 1, dtype=float)
        ok = ~np.isnan(x[:n0])
        X = np.column_stack([np.ones(ok.sum()), t[ok]])
        coef = np.linalg.lstsq(X, x[:n0][ok], rcond=None)[0]
        a0, b0 = float(coef[0]), float(coef[1])
        season = (np.arange(n0) % L)
        dm = np.zeros(L)
        cnt = np.zeros(L)
        trend = a0 + b0 * t
        for i in range(n0):
            if not math.isnan(x[i]):
                dm[season[i]] += (x[i] / trend[i] if not additive else x[i] - trend[i])
                cnt[season[i]] += 1
        s = dm / np.maximum(cnt, 1)
    else:
        yrs = [x[j * L:(j + 1) * L] for j in range(m)]
        means = [float(np.nanmean(y)) for y in yrs]
        b0 = (means[m - 1] - means[0]) / ((m - 1) * L) if m > 1 else 0.0
        a0 = means[0] - (L / 2) * b0
        dev = np.zeros((m, L))
        for j in range(m):
            for i in range(L):
                tr = means[j] + (i + 1 - (L + 1) / 2) * b0
                dev[j, i] = (yrs[j][i] / tr) if not additive else yrs[j][i] - tr
        s = np.nanmean(dev, axis=0)
    if additive:
        s = s - s.mean()
    else:
        s = s * L / s.sum()
    return a0, b0, s


def shw_run(x: np.ndarray, alpha: float, beta: float, gamma: float, a0: float, b0: float,
            s0: np.ndarray, additive: bool):
    n, L = len(x), len(s0)
    seas = list(s0)                      # seas[i] is the seasonal term of period t - L + i
    pred = np.empty(n)
    a, b = a0, b0
    rss = 0.0
    for t in range(n):
        st_ = seas[t % L]
        pred[t] = (a + b) + st_ if additive else (a + b) * st_
        obs = not math.isnan(x[t])
        xt = x[t] if obs else pred[t]
        if obs:
            rss += (x[t] - pred[t]) ** 2
        if additive:
            an = alpha * (xt - st_) + (1 - alpha) * (a + b)
        else:
            an = alpha * xt / st_ + (1 - alpha) * (a + b)
        b = beta * (an - a) + (1 - beta) * b
        if additive:
            seas[t % L] = gamma * (xt - an) + (1 - gamma) * st_
        else:
            seas[t % L] = gamma * xt / an + (1 - gamma) * st_
        a = an
    return pred, rss, a, b, seas


def shw_smooth(x: np.ndarray, period: int, *, additive=False, parms=None, s0=None, sn0=None,
               samp0=None, forecast=0, normalize=False, altstarts=False, start=(0.5, 0.5, 0.5),
               maximize=None) -> dict:
    n_obs = int((~np.isnan(x)).sum())
    n0 = samp0 if samp0 else _first_half(n_obs)
    a0, b0, sea0 = shw_start(x, period, n0, additive, normalize, altstarts)
    if s0 is not None:
        a0, b0 = s0
    if sn0 is not None:
        sea0 = np.asarray(sn0, dtype=float)
    res = None
    if parms is None:
        def fun(z, todo):
            p = 1 / (1 + np.exp(-z))
            return -(shw_run(x, p[0], p[1], p[2], a0, b0, sea0, additive)[1] + _logit_penalty(p)), None, None
        res = maximize(fun, np.log(np.array(start) / (1 - np.array(start))))
        alpha, beta, gamma = [float(1 / (1 + math.exp(-v))) for v in res.b]
    else:
        alpha, beta, gamma = parms
    pred, rss, a, b, seas = shw_run(x, alpha, beta, gamma, a0, b0, sea0, additive)
    n = len(x)
    fc = np.empty(forecast)
    for tau in range(1, forecast + 1):
        sv = seas[(n + tau - 1) % period]
        fc[tau - 1] = (a + b * tau) + sv if additive else (a + b * tau) * sv
    prss = rss + _logit_penalty((alpha, beta, gamma))
    return {"fit": np.concatenate([pred, fc]), "alpha": alpha, "beta": beta, "gamma": gamma,
            "rss": rss, "prss": prss, "rmse": math.sqrt(rss / n_obs) if n_obs else math.nan,
            "N": n_obs, "N_pre": n0, "s1_0": a0, "s2_0": b0, "constant": a, "linear": b,
            "period": period, "opt": res, "seasonal": sea0}

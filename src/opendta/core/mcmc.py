"""Markov chain Monte Carlo: adaptive random-walk Metropolis-Hastings and the summaries of a chain.

Plain functions on arrays (no Session). `adaptive_mh` samples a log posterior in blocks of
parameters: the proposal of a block is normal with covariance scale^2 * Sigma_b, Sigma_b starts
from the curvature at the posterior mode (or an identity matrix) and, during the burn-in, is
re-estimated from the chain and the scale follows the acceptance rate of the block, as in the
adaptive algorithm of [BAYES] bayesmh (every 100 iterations; the exact rules are not documented
and the numbers will not match Stata's: VERIFY). The summaries are the effective sample size by
the initial positive sequence of autocorrelations, the Monte Carlo standard error
sd / sqrt(ESS), equal-tailed and highest-posterior-density intervals, and the
Laplace-Metropolis estimator of the log marginal likelihood."""

from __future__ import annotations

import math

import numpy as np


def autocorr(x: np.ndarray, maxlag: int | None = None) -> np.ndarray:
    """Autocorrelations at lags 0..maxlag by FFT."""
    n = len(x)
    maxlag = n - 1 if maxlag is None else min(maxlag, n - 1)
    xc = x - x.mean()
    v = float(xc @ xc)
    if v <= 0:
        out = np.zeros(maxlag + 1)
        out[0] = 1.0
        return out
    m = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(xc, m)
    ac = np.fft.irfft(f * np.conj(f), m)[:maxlag + 1]
    return ac / ac[0]


def ess(x: np.ndarray) -> float:
    """Effective sample size: n / (1 + 2 sum rho_k) with the sum cut by Geyer's initial positive
    sequence (pairs of consecutive autocorrelations that stay positive)."""
    n = len(x)
    if n < 4 or float(np.var(x)) == 0.0:
        return float(n)
    rho = autocorr(x, min(n - 1, max(100, n // 3)))
    s = 0.0
    k = 1
    while k + 1 < len(rho):
        pair = rho[k] + rho[k + 1]
        if pair < 0:
            break
        s += pair
        k += 2
    tau = 1 + 2 * s
    return float(min(n / max(tau, 1e-12), n * 10.0)) if tau > 0 else float(n)


def hpd(x: np.ndarray, level: float = 95.0) -> tuple[float, float]:
    """Shortest interval that holds level% of the draws."""
    xs = np.sort(x)
    n = len(xs)
    m = max(int(math.ceil(n * level / 100.0)), 1)
    if m >= n:
        return float(xs[0]), float(xs[-1])
    w = xs[m:] - xs[:n - m]
    i = int(np.argmin(w))
    return float(xs[i]), float(xs[i + m])


def equal_tailed(x: np.ndarray, level: float = 95.0) -> tuple[float, float]:
    a = (100.0 - level) / 2.0
    lo, hi = np.percentile(x, [a, 100.0 - a])
    return float(lo), float(hi)


def summarize(draws: np.ndarray, level: float = 95.0, use_hpd: bool = False):
    """Per parameter: mean, sd, mcse, median, interval, ess."""
    rows = []
    for j in range(draws.shape[1]):
        x = draws[:, j]
        e = ess(x)
        sd = float(np.std(x, ddof=1)) if len(x) > 1 else 0.0
        lo, hi = hpd(x, level) if use_hpd else equal_tailed(x, level)
        rows.append(dict(mean=float(x.mean()), sd=sd, mcse=sd / math.sqrt(e) if e > 0 else float("nan"),
                         median=float(np.median(x)), lo=lo, hi=hi, ess=e))
    return rows


def _num_hessian(f, x: np.ndarray, h0: float = 1e-4) -> np.ndarray:
    p = len(x)
    H = np.zeros((p, p))
    f0 = f(x)
    hs = h0 * (1 + np.abs(x))
    for i in range(p):
        for j in range(i, p):
            if i == j:
                xp, xm = x.copy(), x.copy()
                xp[i] += hs[i]
                xm[i] -= hs[i]
                H[i, i] = (f(xp) - 2 * f0 + f(xm)) / hs[i] ** 2
            else:
                vals = []
                for si, sj in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                    y = x.copy()
                    y[i] += si * hs[i]
                    y[j] += sj * hs[j]
                    vals.append(f(y))
                H[i, j] = H[j, i] = (vals[0] - vals[1] - vals[2] + vals[3]) / (4 * hs[i] * hs[j])
    return H


def find_mode(logpost, x0: np.ndarray):
    """Posterior mode by BFGS on -logpost; returns (x, ok)."""
    from scipy import optimize

    def f(x):
        v = logpost(x)
        return -v if np.isfinite(v) else 1e300
    try:
        res = optimize.minimize(f, x0, method="Nelder-Mead" if len(x0) == 1 else "BFGS",
                                options={"maxiter": 2000} if len(x0) == 1 else {"maxiter": 500, "gtol": 1e-6})
        if np.isfinite(logpost(res.x)) and logpost(res.x) >= logpost(x0) - 1e-9:
            return res.x, True
    except Exception:
        pass
    return x0.copy(), False


def proposal_cov(logpost, x: np.ndarray, blocks: list[list[int]]) -> list[np.ndarray]:
    """Initial proposal covariance of every block: inverse curvature at x when it is usable."""
    covs = []
    for b in blocks:
        def fb(z, b=b):
            y = x.copy()
            y[b] = z
            return logpost(y)
        try:
            H = _num_hessian(fb, x[b].copy())
            C = np.linalg.inv(-H)
            C = (C + C.T) / 2
            if np.all(np.isfinite(C)) and np.all(np.linalg.eigvalsh(C) > 1e-14):
                covs.append(C)
                continue
        except np.linalg.LinAlgError:
            pass
        covs.append(np.eye(len(b)) * 0.01 * (1 + float(np.mean(np.abs(x[b]))) ** 2))
    return covs


def adaptive_mh(logpost, x0: np.ndarray, blocks: list[list[int]], nburn: int, nsamp: int,
                thin: int, rng: np.random.Generator, covs: list[np.ndarray] | None = None,
                adapt_every: int = 100, adapt_max: int = 25, scale: float | None = None):
    """Samples logpost. Returns dict(draws, lp, accept (per block), acc_rate, x)."""
    d = len(x0)
    x = x0.copy()
    lp = logpost(x)
    if not np.isfinite(lp):
        raise ValueError("initial values not feasible")
    nb = len(blocks)
    covs = [c.copy() for c in covs] if covs else [np.eye(len(b)) * 0.01 for b in blocks]
    sc = [(scale or 2.38) / math.sqrt(len(b)) for b in blocks]
    chol = [np.linalg.cholesky(c + 1e-14 * np.eye(len(c))) for c in covs]
    target = [0.44 if len(b) == 1 else 0.234 for b in blocks]
    total = nburn + nsamp * thin
    draws = np.empty((nsamp, d))
    lps = np.empty(nsamp)
    acc = np.zeros(nb)
    win = np.zeros(nb)
    nprop = np.zeros(nb)
    hist = np.empty((max(nburn, 1), d))
    adapts = 0
    k = 0
    for it in range(total):
        for bi, b in enumerate(blocks):
            z = rng.standard_normal(len(b))
            y = x.copy()
            y[b] = x[b] + sc[bi] * (chol[bi] @ z)
            lpy = logpost(y)
            nprop[bi] += 1
            if np.isfinite(lpy) and math.log(rng.random()) < lpy - lp:
                x, lp = y, lpy
                acc[bi] += 1
                win[bi] += 1
        if it < nburn:
            hist[it] = x
            if (it + 1) % adapt_every == 0 and adapts < adapt_max:
                adapts += 1
                for bi, b in enumerate(blocks):
                    ar = win[bi] / adapt_every
                    sc[bi] *= math.exp(1.5 * (ar - target[bi]) / max(target[bi], 1e-3) * 0.5)
                    sc[bi] = min(max(sc[bi], 1e-6), 20.0)
                    seg = hist[max(0, it + 1 - 5 * adapt_every):it + 1][:, b]
                    if len(seg) > 2 * len(b) + 2:
                        C = np.atleast_2d(np.cov(seg.T))
                        if np.all(np.isfinite(C)) and np.all(np.linalg.eigvalsh(C) > 1e-14):
                            covs[bi] = 0.5 * covs[bi] + 0.5 * C if adapts < 3 else C
                            chol[bi] = np.linalg.cholesky(covs[bi] + 1e-14 * np.eye(len(b)))
                win[:] = 0
            if it == nburn - 1:
                acc[:] = 0
                nprop[:] = 0
        elif (it - nburn + 1) % thin == 0:
            draws[k] = x
            lps[k] = lp
            k += 1
    rates = acc / np.maximum(nprop, 1)
    return dict(draws=draws, lp=lps, block_rates=rates, acc_rate=float(rates.mean()), x=x)


def laplace_metropolis(draws: np.ndarray, logpost_vals: np.ndarray) -> float:
    """Laplace-Metropolis estimate of the log marginal likelihood: the log posterior
    (likelihood times prior) at the draw with the highest value, plus (d/2) log(2 pi) and
    half the log determinant of the sample covariance of the draws."""
    d = draws.shape[1]
    C = np.atleast_2d(np.cov(draws.T))
    sign, ld = np.linalg.slogdet(C)
    if sign <= 0:
        return float("nan")
    return float(np.max(logpost_vals) + 0.5 * d * math.log(2 * math.pi) + 0.5 * ld)

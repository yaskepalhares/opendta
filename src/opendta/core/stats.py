"""Estatísticas básicas com pesos, nas definições do manual [R] (Methods and
formulas de summarize, centile, pctile, correlate, ttest).

Pesos (w):
    fweight  frequências inteiras: cada obs. conta w vezes
    aweight  pesos analíticos: normalizados para somar N (nº de obs.)
    iweight  pesos "importance": usados como dados, sem normalizar
    pweight  pesos amostrais (tratados como aweight nas médias)
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from . import missing as M


def valid(x: np.ndarray) -> np.ndarray:
    return x < M.SYSMISS


@dataclass
class Moments:
    N: float          # nº de observações (fweight: soma dos pesos)
    sum_w: float
    mean: float
    var: float
    sd: float
    min: float
    max: float
    sum: float
    skewness: float = M.SYSMISS
    kurtosis: float = M.SYSMISS


def moments(x: np.ndarray, w: np.ndarray | None = None, wtype: str = "") -> Moments:
    """x e w já filtrados (sem missing). Variância com divisor N-1."""
    n_obs = len(x)
    if n_obs == 0:
        return Moments(0, 0, M.SYSMISS, M.SYSMISS, M.SYSMISS, M.SYSMISS, M.SYSMISS, 0.0)
    if w is None or wtype == "":
        w = np.ones(n_obs)
        wtype = "fweight"
    w = w.astype(np.float64)
    sum_w = float(w.sum())
    if wtype == "fweight":
        N = sum_w
        ww = w
    elif wtype == "iweight":
        N = sum_w
        ww = w
    else:                                  # aweight / pweight: normaliza para somar n
        N = float(n_obs)
        ww = w * n_obs / sum_w if sum_w else w
    mean = float((ww * x).sum() / ww.sum()) if ww.sum() else M.SYSMISS
    dev = x - mean
    m2 = float((ww * dev ** 2).sum())
    var = m2 / (N - 1) if N > 1 else M.SYSMISS
    sd = math.sqrt(var) if var < M.SYSMISS and var >= 0 else M.SYSMISS
    total = float((ww * x).sum()) if wtype != "aweight" else float((w * x).sum())
    skew = kurt = M.SYSMISS
    W = float(ww.sum())
    if W > 0 and m2 > 0:
        mm2 = m2 / W
        mm3 = float((ww * dev ** 3).sum()) / W
        mm4 = float((ww * dev ** 4).sum()) / W
        skew = mm3 / mm2 ** 1.5
        kurt = mm4 / mm2 ** 2
    return Moments(N, sum_w, mean, var, sd, float(x.min()), float(x.max()), total, skew, kurt)


def percentile(xs: np.ndarray, p: float, w: np.ndarray | None = None) -> float:
    """Percentil p (0–100) na definição de summarize/pctile: com P = W·p/100,
    é x_(i) para o primeiro i com W_i > P, ou a média de x_(i) e x_(i+1)
    quando W_i = P. xs ordenado; w alinhado com xs."""
    n = len(xs)
    if n == 0:
        return M.SYSMISS
    if w is None:
        w = np.ones(n)
    cw = np.cumsum(w)
    W = cw[-1]
    P = W * p / 100.0
    i = int(np.searchsorted(cw, P, side="left"))
    if i >= n:
        return float(xs[-1])
    if abs(cw[i] - P) <= 1e-12 * max(1.0, W) and i + 1 < n:
        return float((xs[i] + xs[i + 1]) / 2)
    if cw[i] <= P and i + 1 < n:
        i += 1
    return float(xs[i])


def centile_value(xs: np.ndarray, p: float) -> float:
    """Definição do comando centile: interpolação em R = (n+1)p/100."""
    n = len(xs)
    if n == 0:
        return M.SYSMISS
    R = (n + 1) * p / 100.0
    r = math.floor(R)
    f = R - r
    if r < 1:
        return float(xs[0])
    if r >= n:
        return float(xs[-1])
    return float(xs[r - 1] + f * (xs[r] - xs[r - 1]))


def t_ci(mean: float, se: float, df: float, level: float) -> tuple[float, float]:
    from scipy import stats as st
    if se >= M.SYSMISS or df <= 0:
        return M.SYSMISS, M.SYSMISS
    t = st.t.ppf(1 - (1 - level / 100) / 2, df)
    return mean - t * se, mean + t * se


def z_crit(level: float) -> float:
    from scipy import stats as st
    return float(st.norm.ppf(1 - (1 - level / 100) / 2))

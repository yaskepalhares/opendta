"""ROC regression ([R] rocreg; Pepe 2003, Janes and Pepe 2009, Alonso and Pepe 2002).

Each case (refvar != 0) is placed in the distribution of the controls' classifier values,
possibly adjusted for control covariates: its percentile value pv is the share of the controls
below it. The ROC curve is ROC(t) = Pr(1 - pv <= t) over the cases, and its summaries are
means and quantiles of pv. The parametric (GLM) model fits
Phi^-1 ROC(t) = i_cons + s_cons Phi^-1(t) + z'b by probit regression of the indicators
1(1 - pv_i <= t_k) on a grid of false-positive rates t_k.

Plain functions of arrays: the command is commands/rocreg.py.
"""

from __future__ import annotations

import numpy as np

from ..lazy import lazy

st = lazy("scipy.stats")


def pv_empirical(case_y: np.ndarray, ctrl_y: np.ndarray, tiecorrected: bool = False) -> np.ndarray:
    """Share of the controls whose value is below each case's; with tiecorrected, ties count
    one half (VERIFY: Stata's treatment of ties without tiecorrected)."""
    c = np.sort(ctrl_y)
    below = np.searchsorted(c, case_y, "left")
    if tiecorrected:
        below = below + 0.5 * (np.searchsorted(c, case_y, "right") - below)
    return below / len(c)


def pv_normal(case_y: np.ndarray, ctrl_y: np.ndarray) -> np.ndarray:
    """Percentile values from a normal distribution of the controls' values."""
    sd = float(np.std(ctrl_y, ddof=1)) if len(ctrl_y) > 1 else 0.0
    return st.norm.cdf((case_y - ctrl_y.mean()) / sd) if sd > 0 else np.full(len(case_y), np.nan)


def pv_linear(case_y, case_X, ctrl_y, ctrl_X, pvc: str = "empirical", tiecorrected: bool = False):
    """Control covariates through a linear model: the controls' values are regressed on the
    covariates (with a constant) and the cases are placed by their residuals, in the empirical
    distribution of the controls' residuals or in a normal one (residual standard deviation
    with n - k degrees of freedom; VERIFY)."""
    Xc = np.column_stack([ctrl_X, np.ones(len(ctrl_y))])
    g, *_ = np.linalg.lstsq(Xc, ctrl_y, rcond=None)
    rc = ctrl_y - Xc @ g
    rd = case_y - np.column_stack([case_X, np.ones(len(case_y))]) @ g
    if pvc == "normal":
        sd = float(np.sqrt((rc ** 2).sum() / max(len(rc) - Xc.shape[1], 1)))
        return st.norm.cdf(rd / sd)
    return pv_empirical(rd, rc, tiecorrected)


def roc_at(pv: np.ndarray, t: float) -> float:
    """ROC(t): the share of the cases with 1 - pv <= t."""
    return float(np.mean(1 - pv <= t + 1e-12))


def invroc_at(pv: np.ndarray, r: float) -> float:
    """The false-positive rate at which the ROC curve reaches the true-positive rate r: the
    smallest t with ROC(t) >= r (VERIFY: Stata's inverse of the empirical ROC)."""
    f = np.sort(1 - pv)
    k = int(np.ceil(r * len(f) - 1e-9))
    return float(f[min(max(k, 1), len(f)) - 1])


def pauc_at(pv: np.ndarray, t: float) -> float:
    """Partial area: the integral of ROC(u) for u from 0 to t, E[(t - (1 - pv))+]."""
    return float(np.mean(np.clip(t - (1 - pv), 0, None)))


def summaries(pv: np.ndarray, spec: list[tuple[str, float | None]]) -> np.ndarray:
    """Values of the requested summaries: ("auc", None), ("roc", t), ("invroc", r),
    ("pauc", t)."""
    out = []
    for kind, x in spec:
        if kind == "auc":
            out.append(float(np.mean(pv)))
        elif kind == "roc":
            out.append(roc_at(pv, x))
        elif kind == "invroc":
            out.append(invroc_at(pv, x))
        else:
            out.append(pauc_at(pv, x))
    return np.array(out)


def fpr_points(n: int) -> np.ndarray:
    """Equally spaced false-positive rates k/(n+1), k = 1..n (VERIFY: Stata's grid)."""
    return np.arange(1, n + 1) / (n + 1)


def probit_glm(pv: np.ndarray, Z: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Coefficients (b for each column of Z, i_cons, s_cons) of the probit regression of
    1(1 - pv_i <= t_k) on z_i, a constant and Phi^-1(t_k), over every case and point."""
    from . import ml
    from .models import _probit
    n, k = len(pv), Z.shape[1]
    q = st.norm.ppf(t)
    y = (1 - pv[:, None] <= t[None, :] + 1e-12).astype(float).ravel()
    X = np.column_stack([np.repeat(Z, len(t), axis=0), np.ones(n * len(t)), np.tile(q, n)])
    fun, _ = _probit(y, X, np.ones(len(y)))
    return ml.maximize(fun, np.zeros(k + 2)).b

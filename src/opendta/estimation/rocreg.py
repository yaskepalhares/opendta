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

import math

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


# ---------------------------------------------------------------------------
# parametric (binormal) ROC model by maximum likelihood: rocreg, probit ml
# ---------------------------------------------------------------------------

def adjust_split(theta: np.ndarray, kx: int, kz: int):
    """theta = (casecov: X then _cons, casesd, ctrlcov: Z then _cons, ctrlsd)."""
    return (theta[:kx + 1], theta[kx + 1], theta[kx + 2:kx + 3 + kz], theta[kx + 3 + kz])


def adjust_derivatives(theta: np.ndarray, y: np.ndarray, d: np.ndarray, X: np.ndarray, Z: np.ndarray):
    """The covariate control adjustment model of rocreg, probit ml: the controls are
    N(Z1 gamma, ctrlsd^2) and the cases N(Z1 gamma + X1 delta, casesd^2), with X1 and Z1 the ROC
    and control covariates followed by a constant. Returns the log likelihood of each observation,
    the scores (one row per observation) and the Hessian of the total."""
    kx, kz = X.shape[1], Z.shape[1]
    delta, sd_d, gamma, sd_c = adjust_split(theta, kx, kz)
    n, k = len(y), len(theta)
    X1 = np.column_stack([X, np.ones(n)])
    Z1 = np.column_stack([Z, np.ones(n)])
    mu = Z1 @ gamma + np.where(d, X1 @ delta, 0.0)
    sd = np.where(d, sd_d, sd_c)
    if sd_d <= 0 or sd_c <= 0:
        return np.full(n, -np.inf), np.zeros((n, k)), np.zeros((k, k))
    r = y - mu
    ll = -0.5 * math.log(2 * math.pi) - np.log(sd) - r * r / (2 * sd * sd)
    # derivative of the mean (W) and of the standard deviation (Sg) with respect to theta
    W = np.zeros((n, k))
    W[:, :kx + 1] = X1 * d[:, None]
    W[:, kx + 2:kx + 3 + kz] = Z1
    Sg = np.zeros((n, k))
    Sg[d, kx + 1] = 1.0
    Sg[~d, kx + 3 + kz] = 1.0
    dmu = r / sd ** 2
    dsd = -1 / sd + r * r / sd ** 3
    S = W * dmu[:, None] + Sg * dsd[:, None]
    H = -(W * (1 / sd ** 2)[:, None]).T @ W
    cross = (W * (-2 * r / sd ** 3)[:, None]).T @ Sg
    H += cross + cross.T
    H += (Sg * (1 / sd ** 2 - 3 * r * r / sd ** 4)[:, None]).T @ Sg
    return ll, S, H


def adjust_start(y: np.ndarray, d: np.ndarray, kx: int, kz: int) -> np.ndarray:
    """Starting values: no covariate effects, the controls' mean, the difference of the means
    and the standard deviations of the two groups (n - 1 divisor; Stata's iteration 0, compat
    1912)."""
    yc, yd = y[~d], y[d]
    return np.concatenate([np.zeros(kx), [yd.mean() - yc.mean(), yd.std(ddof=1)],
                           np.zeros(kz), [yc.mean(), yc.std(ddof=1)]])


def adjust_ml(y: np.ndarray, d: np.ndarray, X: np.ndarray, Z: np.ndarray):
    """Newton-Raphson fit of the adjustment model. Returns the ml.MLResult."""
    from . import ml

    def fun(t, todo):
        ll, S, H = adjust_derivatives(t, y, d, X, Z)
        v = float(ll.sum())
        return (v if np.isfinite(v) else -np.inf), S.sum(axis=0), H
    return ml.maximize(fun, adjust_start(y, d, X.shape[1], Z.shape[1]))


def roc_params(theta: np.ndarray, kx: int, kz: int) -> np.ndarray:
    """The ROC model Phi^-1 ROC(t) = i_cons + x'b + s_cons Phi^-1(t) of the adjustment model:
    (i_cons, b, s_cons), and the AUC when there are no ROC covariates."""
    delta, sd_d, _, sd_c = adjust_split(theta, kx, kz)
    i, b, sl = delta[-1] / sd_d, delta[:-1] / sd_d, sd_c / sd_d
    out = np.concatenate([[i], b, [sl]])
    if not kx:
        out = np.append(out, st.norm.cdf(i / math.sqrt(1 + sl * sl)))
    return out


def roc_jacobian(theta: np.ndarray, kx: int, kz: int) -> np.ndarray:
    """Derivatives of roc_params by central differences."""
    base = roc_params(theta, kx, kz)
    J = np.zeros((len(base), len(theta)))
    for j in range(len(theta)):
        h = 1e-6 * max(1.0, abs(theta[j]))
        e = np.zeros(len(theta))
        e[j] = h
        J[:, j] = (roc_params(theta + e, kx, kz) - roc_params(theta - e, kx, kz)) / (2 * h)
    return J

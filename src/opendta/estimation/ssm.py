"""Linear Gaussian state-space models (Kalman filter) shared by ucm, sspace and dfactor.

    y_t = D z_t + d_t + e_t,        e_t ~ N(0, H)
    z_t = A z_{t-1} + c + w_t,      w_t ~ N(0, S)

Initial state: mean a0 and covariance P0. States flagged diffuse get a very large variance and the
observations they make uninformative are left out of the likelihood (the limit of De Jong's diffuse
filter)."""

from __future__ import annotations

import numpy as np
from scipy import linalg

KAPPA = 1e7


def stationary_cov(A: np.ndarray, S: np.ndarray) -> np.ndarray | None:
    """Solution of P = A P A' + S, or None if A has eigenvalues on or outside the unit circle."""
    if A.size and np.max(np.abs(np.linalg.eigvals(A))) >= 1 - 1e-9:
        return None
    try:
        return linalg.solve_discrete_lyapunov(A, S)
    except (np.linalg.LinAlgError, ValueError):
        return None


def initial_state(A: np.ndarray, S: np.ndarray, diffuse: np.ndarray | None = None):
    """(a0, P0): diffuse states get variance KAPPA; the rest the stationary covariance (zero
    mean). If nothing is declared diffuse and A is not stable, every state is diffuse."""
    r = A.shape[0]
    mask = np.zeros(r, dtype=bool) if diffuse is None else np.asarray(diffuse, dtype=bool)
    P0 = np.zeros((r, r))
    st = ~mask
    if st.any():
        Pst = stationary_cov(A[np.ix_(st, st)], S[np.ix_(st, st)])
        if Pst is None:
            mask[:] = True
            st = ~mask
        else:
            P0[np.ix_(st, st)] = Pst
    P0[np.ix_(mask, mask)] = KAPPA * np.eye(int(mask.sum()))
    return np.zeros(r), P0, mask


def kalman_loglik(y: np.ndarray, D: np.ndarray, A: np.ndarray, S: np.ndarray, H: np.ndarray,
                  c: np.ndarray | None = None, d: np.ndarray | None = None,
                  diffuse: np.ndarray | None = None, return_states: bool = False):
    """Log-likelihood contributions (length T; zero for the missing observations, only the constant for the diffuse ones).
    y: (T, m) with NaN for missing. Returns None when the filter breaks down."""
    T, m = y.shape
    r = A.shape[0]
    a, P, mask = initial_state(A, S, diffuse)
    c = np.zeros(r) if c is None else c
    d = np.zeros((T, m)) if d is None else d
    out = np.zeros(T)
    fa = np.zeros((T, r)) if return_states else None
    fP = np.zeros((T, r, r)) if return_states else None
    log2pi = np.log(2 * np.pi)
    for t in range(T):
        a = A @ a + c
        P = A @ P @ A.T + S
        obs = np.isfinite(y[t])
        if obs.any():
            Dt = D[obs]
            v = y[t, obs] - Dt @ a - d[t, obs]
            F = Dt @ P @ Dt.T + H[np.ix_(obs, obs)]
            try:
                Fi = np.linalg.inv(F)
            except np.linalg.LinAlgError:
                return None
            sign, ld = np.linalg.slogdet(F)
            if sign <= 0:
                return None
            if not (mask.any() and np.max(np.diag(Dt @ P @ Dt.T)) >= KAPPA / 100):
                out[t] = -0.5 * (obs.sum() * log2pi + ld + float(v @ Fi @ v))
            else:
                out[t] = -0.5 * obs.sum() * log2pi      # VERIFY (case 1904): a diffuse observation keeps its constant
            K = P @ Dt.T @ Fi
            a = a + K @ v
            P = P - K @ F @ K.T
            P = (P + P.T) / 2
        if return_states:
            fa[t], fP[t] = a, P
    if not np.all(np.isfinite(out)):
        return None
    return (out, fa, fP) if return_states else out


def kalman_onestep(y: np.ndarray, D: np.ndarray, A: np.ndarray, S: np.ndarray, H: np.ndarray,
                   c: np.ndarray | None = None, d: np.ndarray | None = None,
                   diffuse: np.ndarray | None = None):
    """One-step-ahead predictions of y (T, m) and the variances of the prediction errors, from the
    same filter as kalman_loglik; missing observations are predicted but do not update the
    states."""
    T, m = y.shape
    r = A.shape[0]
    a, P, _mask = initial_state(A, S, diffuse)
    c = np.zeros(r) if c is None else c
    d = np.zeros((T, m)) if d is None else d
    yhat = np.full((T, m), np.nan)
    fvar = np.full((T, m), np.nan)
    for t in range(T):
        a = A @ a + c
        P = A @ P @ A.T + S
        yhat[t] = D @ a + d[t]
        Fall = D @ P @ D.T + H
        fvar[t] = np.diag(Fall)
        obs = np.isfinite(y[t])
        if obs.any():
            Dt = D[obs]
            v = y[t, obs] - Dt @ a - d[t, obs]
            F = Fall[np.ix_(obs, obs)]
            K = P @ Dt.T @ np.linalg.inv(F)
            a = a + K @ v
            P = P - K @ F @ K.T
            P = (P + P.T) / 2
    return yhat, fvar

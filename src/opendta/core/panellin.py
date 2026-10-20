"""Linear panel estimators: feasible GLS and panel-corrected standard errors ([XT] xtgls,
[XT] xtpcse), AR(1) disturbances ([XT] xtregar), random coefficients ([XT] xtrc),
Hausman-Taylor ([XT] xthtaylor) and panel instrumental variables ([XT] xtivreg).

Plain functions on arrays: no Session, no variables, no printing. Observations are rows; `inv`
is the group number (0..G-1) of every row and `tidx` the integer time index inside the panel."""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# panel structure
# ---------------------------------------------------------------------------

class PanelIndex:
    """Group structure and the consecutive-period lag of every observation."""

    def __init__(self, inv: np.ndarray, tidx: np.ndarray):
        self.inv = np.asarray(inv)
        self.tidx = np.asarray(tidx)
        self.G = int(self.inv.max()) + 1
        self.counts = np.bincount(self.inv, minlength=self.G)
        n = len(self.inv)
        order = np.lexsort((self.tidx, self.inv))
        self.order = order
        lag = np.full(n, -1, dtype=int)
        gi, ti = self.inv[order], self.tidx[order]
        same = (gi[1:] == gi[:-1]) & (ti[1:] - ti[:-1] == 1)
        lag[order[1:][same]] = order[:-1][same]
        self.lag = lag
        self.has_lag = lag >= 0

    def lagged(self, v: np.ndarray) -> np.ndarray:
        out = np.zeros_like(v, dtype=float)
        m = self.has_lag
        out[m] = v[self.lag[m]]
        return out

    def group_sums(self, v: np.ndarray) -> np.ndarray:
        if v.ndim == 1:
            return np.bincount(self.inv, weights=v, minlength=self.G)
        out = np.zeros((self.G, v.shape[1]))
        np.add.at(out, self.inv, v)
        return out

    def group_means(self, v: np.ndarray) -> np.ndarray:
        s = self.group_sums(v)
        return s / (self.counts if v.ndim == 1 else self.counts[:, None])

    def within(self, v: np.ndarray) -> np.ndarray:
        return v - self.group_means(v)[self.inv]

    @property
    def balanced(self) -> bool:
        return bool(self.counts.min() == self.counts.max())


def prais_winsten(v: np.ndarray, pi: PanelIndex, rho_obs: np.ndarray) -> np.ndarray:
    """z*_t = z_t - rho z_{t-1}; the first observation of a panel (or after a gap) is scaled by
    sqrt(1 - rho^2). v may be 1-D or 2-D."""
    lagv = np.zeros_like(v, dtype=float)
    m = pi.has_lag
    lagv[m] = v[pi.lag[m]]
    r = np.clip(rho_obs, -0.9999, 0.9999)      # VERIFY: |rho| >= 1 is truncated
    r = r if v.ndim == 1 else r[:, None]
    mm = m if v.ndim == 1 else m[:, None]
    return np.where(mm, v - r * lagv, np.sqrt(1 - r ** 2) * v)


# ---------------------------------------------------------------------------
# autocorrelation estimators
# ---------------------------------------------------------------------------

RHOTYPES = ("regress", "freg", "tscorr", "dw", "theil", "nagar")


def rho_from_residuals(e: np.ndarray, pi: PanelIndex, kind: str = "regress", k: int = 0) -> np.ndarray:
    """AR(1) coefficient of every panel from its residuals; `kind` as rhotype()."""
    m = pi.has_lag
    a = np.zeros_like(e)
    a[m] = e[pi.lag[m]]
    b = np.where(m, e, 0.0)
    ab = pi.group_sums(a * b)
    a2 = pi.group_sums(a * a)
    b2 = pi.group_sums(b * b)
    e2 = pi.group_sums(e * e)
    d2 = pi.group_sums((b - a) ** 2 * m)
    with np.errstate(divide="ignore", invalid="ignore"):
        if kind == "regress":
            r = ab / a2
        elif kind == "freg":
            r = ab / b2
        elif kind == "tscorr":
            r = ab / e2
        elif kind == "dw":
            r = 1 - d2 / e2 / 2
        elif kind == "theil":
            T = pi.counts.astype(float)
            r = ab / e2 * (T - k) / T
        elif kind == "nagar":
            T = pi.counts.astype(float)
            rdw = 1 - d2 / e2 / 2
            r = (rdw * T ** 2 + k ** 2) / (T ** 2 - k ** 2)
        else:
            raise ValueError(kind)
    return np.nan_to_num(r)


def common_rho(rho_i: np.ndarray, pi: PanelIndex, np1: bool = False) -> float:
    """Average of the panel autocorrelations weighted by T_i - 1 (T_i with np1)."""
    w = pi.counts.astype(float) - (0 if np1 else 1)
    return float((w * rho_i).sum() / w.sum()) if w.sum() > 0 else 0.0


# ---------------------------------------------------------------------------
# FGLS (xtgls) and PCSE (xtpcse)
# ---------------------------------------------------------------------------

def _present_times(pi: PanelIndex):
    """For every time index: the observation of every panel present (rows of obs, panels)."""
    tmax = int(pi.tidx.max()) + 1
    obs_by_t = [[] for _ in range(tmax)]
    pan_by_t = [[] for _ in range(tmax)]
    for i in range(len(pi.inv)):
        obs_by_t[pi.tidx[i]].append(i)
        pan_by_t[pi.tidx[i]].append(pi.inv[i])
    return [(np.array(o, dtype=int), np.array(p, dtype=int)) for o, p in zip(obs_by_t, pan_by_t) if len(o)]


def sigma_correlated(e: np.ndarray, pi: PanelIndex, pairwise: bool = False) -> tuple[np.ndarray, int]:
    """Panel covariance matrix S_ij = e_i'e_j / T (casewise: the periods in which every panel is
    observed; pairwise: the periods panels i and j share). Returns S and the number of periods."""
    G = pi.G
    times = _present_times(pi)
    E = np.full((len(times), G), np.nan)
    for r, (obs, pan) in enumerate(times):
        E[r, pan] = e[obs]
    full = ~np.isnan(E).any(axis=1)
    if not pairwise:
        Ef = E[full]
        T = int(full.sum())
        return (Ef.T @ Ef) / max(T, 1), T
    S = np.zeros((G, G))
    for i in range(G):
        for j in range(i, G):
            ok = ~np.isnan(E[:, i]) & ~np.isnan(E[:, j])
            S[i, j] = S[j, i] = float((E[ok, i] * E[ok, j]).sum() / max(int(ok.sum()), 1))
    return S, int(full.sum())


def gls_fit(y: np.ndarray, X: np.ndarray, pi: PanelIndex, panels: str, e: np.ndarray,
            pairwise: bool = False, nmk: bool = False) -> dict:
    """One FGLS step from the residuals `e`: Sigma is estimated according to `panels`
    ('iid', 'hetero', 'correlated') and the GLS estimator and its covariance computed.
    Returns b, V, Sigma, the number of estimated covariances."""
    n, k = X.shape
    scale = n / (n - k) if nmk else 1.0
    if panels == "iid":
        s2 = float(e @ e) / n * scale
        A = np.linalg.inv(X.T @ X)
        b = A @ (X.T @ y)
        return {"b": b, "V": s2 * A, "sigma": np.array([[s2]]), "ncov": 1, "T": None}
    if panels == "hetero":
        s2 = pi.group_sums(e * e) / pi.counts * scale
        w = 1 / s2[pi.inv]
        A = np.linalg.inv((X * w[:, None]).T @ X)
        b = A @ ((X * w[:, None]).T @ y)
        return {"b": b, "V": A, "sigma": np.diag(s2), "ncov": pi.G, "T": None}
    S, T = sigma_correlated(e, pi, pairwise)
    S = S * scale
    G = pi.G
    if T < G:
        from .errors import StataError
        raise StataError(498, "panels(correlated) requires at least as many time periods as panels")   # VERIFY
    XtWX = np.zeros((k, k))
    XtWy = np.zeros(k)
    for obs, pan in _present_times(pi):
        Si = np.linalg.inv(S[np.ix_(pan, pan)])
        Xt = X[obs]
        XtWX += Xt.T @ Si @ Xt
        XtWy += Xt.T @ Si @ y[obs]
    A = np.linalg.inv(XtWX)
    return {"b": A @ XtWy, "V": A, "sigma": S, "ncov": G * (G + 1) // 2, "T": T}


def gls_loglik(e: np.ndarray, pi: PanelIndex, panels: str, sigma: np.ndarray, rho_g: np.ndarray | None) -> float:
    """Gaussian log likelihood at the final residuals and covariances."""
    n = len(e)
    ll = -0.5 * n * np.log(2 * np.pi)
    if panels == "iid":
        s2 = float(sigma[0, 0])
        ll += -0.5 * (n * np.log(s2) + float(e @ e) / s2)
    elif panels == "hetero":
        s2 = np.diag(sigma)
        ll += -0.5 * float((pi.counts * np.log(s2)).sum() + (pi.group_sums(e * e) / s2).sum())
    else:
        tot = 0.0
        for obs, pan in _present_times(pi):
            Ssub = sigma[np.ix_(pan, pan)]
            tot += np.linalg.slogdet(Ssub)[1] + float(e[obs] @ np.linalg.solve(Ssub, e[obs]))
        ll += -0.5 * tot
    if rho_g is not None:
        ll += 0.5 * float(np.log(1 - np.clip(rho_g, -0.9999, 0.9999) ** 2).sum())
    return float(ll)


def pcse_covariance(X: np.ndarray, e: np.ndarray, pi: PanelIndex, mode: str, pairwise: bool = False):
    """(X'X)^-1 X' Omega X (X'X)^-1 with Omega from `mode`: 'correlated' (full panel covariance),
    'hetonly' (panel variances) or 'independent' (one variance). Returns V and the number of
    covariances estimated."""
    n, k = X.shape
    A = np.linalg.inv(X.T @ X)
    if mode == "independent":
        s2 = float(e @ e) / n
        return A * s2, 1
    if mode == "hetonly":
        s2 = pi.group_sums(e * e) / pi.counts
        meat = (X * s2[pi.inv][:, None]).T @ X
        return A @ meat @ A, pi.G
    S, T = sigma_correlated(e, pi, pairwise)
    meat = np.zeros((k, k))
    for obs, pan in _present_times(pi):
        Xt = X[obs]
        meat += Xt.T @ S[np.ix_(pan, pan)] @ Xt
    return A @ meat @ A, pi.G * (pi.G + 1) // 2


# ---------------------------------------------------------------------------
# xtrc: Swamy random coefficients
# ---------------------------------------------------------------------------

def swamy_rc(y: np.ndarray, X: np.ndarray, pi: PanelIndex) -> dict:
    """Swamy's random-coefficients estimator. Panels with at most k observations are not
    estimable (the caller drops them). Returns b, V, Delta, the group estimates and the
    parameter-constancy statistic."""
    n, k = X.shape
    G = pi.G
    bi = np.zeros((G, k))
    Vi = np.zeros((G, k, k))
    for g in range(G):
        ix = np.flatnonzero(pi.inv == g)
        Xg, yg = X[ix], y[ix]
        A = np.linalg.inv(Xg.T @ Xg)
        b = A @ (Xg.T @ yg)
        r = yg - Xg @ b
        s2 = float(r @ r) / (len(ix) - k)
        bi[g], Vi[g] = b, s2 * A
    bbar = bi.mean(axis=0)
    dev = bi - bbar
    first = dev.T @ dev / (G - 1)
    delta = first - Vi.mean(axis=0)
    # the second term is dropped when the difference is not nonnegative definite
    if np.min(np.linalg.eigvalsh((delta + delta.T) / 2)) < -1e-12:
        delta = first
    Ws = np.array([np.linalg.inv(delta + Vi[g]) for g in range(G)])
    Wsum = Ws.sum(axis=0)
    V = np.linalg.inv(Wsum)
    b = V @ sum(Ws[g] @ bi[g] for g in range(G))
    # test of parameter constancy
    Pi = np.array([np.linalg.inv(Vi[g]) for g in range(G)])
    bstar = np.linalg.solve(Pi.sum(axis=0), sum(Pi[g] @ bi[g] for g in range(G)))
    chi2 = float(sum((bi[g] - bstar) @ Pi[g] @ (bi[g] - bstar) for g in range(G)))
    return {"b": b, "V": (V + V.T) / 2, "delta": delta, "bi": bi, "Vi": Vi, "chi2": chi2, "df": k * (G - 1)}


# ---------------------------------------------------------------------------
# two-stage least squares on transformed data
# ---------------------------------------------------------------------------

def tsls(y: np.ndarray, X: np.ndarray, Z: np.ndarray) -> dict:
    """b = (X'Pz X)^-1 X'Pz y; also returns (X'Pz X)^-1 and the residuals y - Xb."""
    ZtZ = Z.T @ Z
    PzX = Z @ np.linalg.solve(ZtZ, Z.T @ X)
    A = np.linalg.inv(PzX.T @ X)
    b = A @ (PzX.T @ y)
    return {"b": b, "A": A, "e": y - X @ b, "PzX": PzX}


def variance_components_iv(y, X, Z, pi: PanelIndex, k_between: int, nosa: bool = False) -> tuple[float, float]:
    """sigma_e^2 and sigma_u^2 for the random-effects IV estimators: the first from the within
    2SLS, the second from the between 2SLS (adapted Swamy-Arora). X and Z carry no constant for
    the within step; the between step adds it. Returns (sigma_e2, sigma_u2)."""
    n = len(y)
    G = pi.G
    yw, Xw, Zw = pi.within(y), pi.within(X), pi.within(Z)
    kx = X.shape[1]
    fit = tsls(yw, Xw, Zw)
    # residuals use the actual regressors, with the constant of the within fit absorbed
    e_w = fit["e"]
    sigma_e2 = float(e_w @ e_w) / (n - G - kx)
    yb, Xb, Zb = pi.group_means(y), pi.group_means(X), pi.group_means(Z)
    Xb1 = np.column_stack([Xb, np.ones(G)])
    Zb1 = np.column_stack([Zb, np.ones(G)])
    fb = tsls(yb, Xb1, Zb1)
    e_b = fb["e"]
    Tbar = G / float(np.sum(1.0 / pi.counts))
    s_b2 = float(e_b @ e_b) / (G - k_between)
    sigma_u2 = max(s_b2 - sigma_e2 / Tbar, 0.0)
    return sigma_e2, sigma_u2


# ---------------------------------------------------------------------------
# Hausman-Taylor
# ---------------------------------------------------------------------------

def hausman_taylor(y: np.ndarray, X1: np.ndarray, X2: np.ndarray, Z1: np.ndarray, Z2: np.ndarray,
                   pi: PanelIndex, amacurdy: bool = False) -> dict:
    """Hausman-Taylor estimator. X1/X2: time-varying exogenous/endogenous; Z1/Z2: time-invariant
    exogenous/endogenous (without the constant, which is always included). Returns the
    coefficients in the order X1, X2, Z1, Z2, cons, the covariance, sigma_u, sigma_e."""
    n, G = len(y), pi.G
    X = np.column_stack([X1, X2]) if X1.shape[1] + X2.shape[1] else np.zeros((n, 0))
    Z = np.column_stack([Z1, Z2]) if Z1.shape[1] + Z2.shape[1] else np.zeros((n, 0))
    kx, kz = X.shape[1], Z.shape[1]
    yw, Xw = pi.within(y), pi.within(X)
    bw = np.linalg.solve(Xw.T @ Xw, Xw.T @ yw)
    e_w = yw - Xw @ bw
    sigma_e2 = float(e_w @ e_w) / (n - G - kx)
    # between equation: d_i = ybar_i - xbar_i bw on [Z, 1] with the instruments [X1bar, Z1, 1]
    d = pi.group_means(y) - pi.group_means(X) @ bw
    # one row per group: time-invariant variables take the first row of every group
    fi = np.zeros(G, dtype=int)
    fi[pi.inv[::-1]] = np.arange(n)[::-1]
    Zg = Z[fi] if kz else np.zeros((G, 0))
    X1b = pi.group_means(X1) if X1.shape[1] else np.zeros((G, 0))
    Z1g = Z1[fi] if Z1.shape[1] else np.zeros((G, 0))
    W = np.column_stack([Zg, np.ones(G)])
    Q = np.column_stack([X1b, Z1g, np.ones(G)])
    fit = tsls(d, W, Q)
    rss = float(fit["e"] @ fit["e"])
    Tbar = G / float(np.sum(1.0 / pi.counts))
    sigma_b2 = rss / G
    sigma_u2 = max(sigma_b2 - sigma_e2 / Tbar, 0.0)
    theta = 1 - np.sqrt(sigma_e2 / (sigma_e2 + pi.counts * sigma_u2))
    th = theta[pi.inv]
    ones = np.ones(n)

    def tr(v):
        if v.ndim == 1:
            return v - th * pi.group_means(v)[pi.inv]
        return v - th[:, None] * pi.group_means(v)[pi.inv]
    R = np.column_stack([X, Z, ones])
    Rs = tr(R)
    ys = tr(y)
    # instruments: within deviations of every time-varying variable, group means of X1, Z1, constant
    ins = [Xw, X1b[pi.inv] if X1.shape[1] else np.zeros((n, 0)), Z1, ones]
    if amacurdy:
        # Amemiya-MaCurdy: the values of X1 in every period of the panel (balanced panels)
        order = np.lexsort((pi.tidx, pi.inv))
        T = int(pi.counts[0])
        X1o = X1[order].reshape(G, T, X1.shape[1])
        wide = X1o.reshape(G, T * X1.shape[1])
        ins = [Xw, wide[pi.inv], Z1, ones]
    I = np.column_stack(ins)
    f2 = tsls(ys, Rs, I)
    return {"b": f2["b"], "V": f2["A"], "sigma_u": float(np.sqrt(sigma_u2)), "sigma_e": float(np.sqrt(sigma_e2)),
            "e": f2["e"], "theta": theta, "bw": bw}


# ---------------------------------------------------------------------------
# panel instrumental variables ([XT] xtivreg)
# ---------------------------------------------------------------------------

def panel_iv(y: np.ndarray, X: np.ndarray, Z: np.ndarray, pi: PanelIndex, method: str,
             nosa: bool = False) -> dict:
    """Panel 2SLS. X: regressors (endogenous and exogenous) and Z: instruments (the exogenous
    regressors and the excluded instruments), both without the constant. method: g2sls and
    ec2sls (random effects), fe (within), be (between), fd (first differences).

    Returns b (the constant last, except in fd), the matrix A with V = s2 * A, the residuals of
    the fitted (transformed) equation, s2 and the variance components. Rows of the equation are
    the rows of the transformed data (`rows` indexes the original ones, `groups` their panel)."""
    n, G, k = len(y), pi.G, X.shape[1]
    if method == "fd":
        m = pi.has_lag
        dy = (y - pi.lagged(y))[m]
        dX = (X - np.column_stack([pi.lagged(X[:, j]) for j in range(k)]))[m]
        dZ = (Z - np.column_stack([pi.lagged(Z[:, j]) for j in range(Z.shape[1])]))[m]
        f = tsls(dy, dX, dZ)
        ne = int(m.sum())
        s2 = float(f["e"] @ f["e"]) / (ne - k)
        return dict(b=f["b"], A=f["A"], e=f["e"], s2=s2, rows=np.flatnonzero(m), df=ne - k, n=ne)
    if method == "be":
        yb, Xb, Zb = pi.group_means(y), pi.group_means(X), pi.group_means(Z)
        X1, Z1 = np.column_stack([Xb, np.ones(G)]), np.column_stack([Zb, np.ones(G)])
        f = tsls(yb, X1, Z1)
        s2 = float(f["e"] @ f["e"]) / (G - k - 1)
        return dict(b=f["b"], A=f["A"], e=f["e"], s2=s2, df=G - k - 1, n=G)
    yw, Xw, Zw = pi.within(y), pi.within(X), pi.within(Z)
    fw = tsls(yw, Xw, Zw)
    s2w = float(fw["e"] @ fw["e"]) / (n - G - k)
    if method == "fe":
        xbar = X.mean(axis=0)
        cons = float(y.mean() - xbar @ fw["b"])
        return dict(b=np.append(fw["b"], cons), A=fw["A"], e=fw["e"], s2=s2w, df=n - G - k, xbar=xbar, cons=cons,
                    n=n, sigma_e2=s2w)
    se2, su2 = variance_components_iv(y, X, Z, pi, k + 1, nosa)
    theta = 1 - np.sqrt(se2 / (se2 + pi.counts * su2))
    th = theta[pi.inv]

    def tr(v):
        gm = pi.group_means(v)[pi.inv]
        return v - (th * gm if v.ndim == 1 else th[:, None] * gm)
    ys = tr(y)
    Xs = np.column_stack([tr(X), 1 - th])
    if method == "g2sls":
        Zs = np.column_stack([tr(Z), 1 - th])
    else:
        Zs = np.column_stack([Zw, pi.group_means(Z)[pi.inv], np.ones(n)])
    f = tsls(ys, Xs, Zs)
    return dict(b=f["b"], A=f["A"], e=f["e"], s2=se2, df=n - k - 1, sigma_e2=se2, sigma_u2=su2, theta=theta, n=n)

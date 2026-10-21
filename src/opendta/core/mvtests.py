"""Computations behind Hotelling's test, MANOVA, canonical correlations and the mvtest family
([MV] hotelling, [MV] manova, [MV] canon, [MV] mvtest).

Plain functions on arrays and numbers: no Session, no variables, no printing."""

from __future__ import annotations

import math

import numpy as np

from ..lazy import lazy

st = lazy("scipy.stats")


# ---------------------------------------------------------------------------
# Hotelling's T-squared
# ---------------------------------------------------------------------------

def hotelling_one(xbar: np.ndarray, S: np.ndarray, n: float, mu0: np.ndarray | None = None) -> dict:
    """One-group T2 = n (xbar - mu0)' S^-1 (xbar - mu0) and its F transformation."""
    p = len(xbar)
    d = xbar - (0.0 if mu0 is None else mu0)
    t2 = float(n * d @ np.linalg.solve(S, d))
    F = (n - p) / ((n - 1.0) * p) * t2
    return {"T2": t2, "F": F, "df1": p, "df2": n - p, "p": float(st.f.sf(F, p, n - p))}


def hotelling_two(x1: np.ndarray, S1: np.ndarray, n1: float, x2: np.ndarray, S2: np.ndarray, n2: float) -> dict:
    """Two-group T2 with the pooled covariance matrix."""
    p = len(x1)
    Sp = ((n1 - 1.0) * S1 + (n2 - 1.0) * S2) / (n1 + n2 - 2.0)
    d = x1 - x2
    t2 = float(n1 * n2 / (n1 + n2) * d @ np.linalg.solve(Sp, d))
    F = (n1 + n2 - p - 1.0) / ((n1 + n2 - 2.0) * p) * t2
    return {"T2": t2, "F": F, "df1": p, "df2": n1 + n2 - p - 1.0, "p": float(st.f.sf(F, p, n1 + n2 - p - 1.0))}


# ---------------------------------------------------------------------------
# multivariate test statistics from the eigenvalues of E^-1 H
# ---------------------------------------------------------------------------

def mv_statistics(lam: np.ndarray, p: int, q: float, v: float) -> list[dict]:
    """Wilks' lambda, Pillai's trace, Lawley-Hotelling trace and Roy's largest root with their F
    approximations. `lam`: nonzero eigenvalues of E^-1 H (descending); p: number of dependent
    variables; q: hypothesis degrees of freedom; v: error degrees of freedom. Each result has
    the statistic, the F test (df1, df2, F, p) and a flag: e (exact), a (approximate) or u
    (upper bound on F)."""
    lam = np.asarray(lam, dtype=float)
    lam = lam[lam > 0] if len(lam) else lam
    s = min(p, q)
    m = (abs(p - q) - 1.0) / 2.0
    n = (v - p - 1.0) / 2.0
    out = []
    wilks = float(np.prod(1.0 / (1.0 + lam))) if len(lam) else 1.0
    if p * p + q * q - 5.0 > 0:
        t = math.sqrt((p * p * q * q - 4.0) / (p * p + q * q - 5.0))
    else:
        t = 1.0
    df1 = p * q
    df2 = t * (v - (p - q + 1.0) / 2.0) - (p * q - 2.0) / 2.0
    root = wilks ** (1.0 / t) if wilks > 0 else 0.0
    F = (1.0 - root) / root * df2 / df1 if root > 0 else float("inf")
    out.append({"name": "W", "stat": wilks, "df1": float(df1), "df2": df2, "F": F,
                "p": float(st.f.sf(F, df1, df2)) if np.isfinite(F) and df2 > 0 else float("nan"),
                "flag": "e" if (p <= 2 or q <= 2) else "a"})
    pil = float(np.sum(lam / (1.0 + lam)))
    df1 = s * (2.0 * m + s + 1.0)
    df2 = s * (2.0 * n + s + 1.0)
    F = (2.0 * n + s + 1.0) / (2.0 * m + s + 1.0) * pil / (s - pil) if s - pil > 0 else float("inf")
    out.append({"name": "P", "stat": pil, "df1": df1, "df2": df2, "F": F,
                "p": float(st.f.sf(F, df1, df2)) if np.isfinite(F) and df2 > 0 else float("nan"),
                "flag": "e" if s == 1 else "a"})
    lh = float(lam.sum())
    df1 = s * (2.0 * m + s + 1.0)
    df2 = 2.0 * (s * n + 1.0)
    F = df2 * lh / (s * s * (2.0 * m + s + 1.0))
    out.append({"name": "L", "stat": lh, "df1": df1, "df2": df2, "F": F,
                "p": float(st.f.sf(F, df1, df2)) if df2 > 0 else float("nan"), "flag": "e" if s == 1 else "a"})
    l1 = float(lam[0]) if len(lam) else 0.0
    d = max(p, q)
    df1 = float(d)
    df2 = v - d + q
    F = l1 * df2 / d
    out.append({"name": "R", "stat": l1, "df1": df1, "df2": df2, "F": F,
                "p": float(st.f.sf(F, df1, df2)) if df2 > 0 else float("nan"), "flag": "e" if s == 1 else "u"})
    return out


def eig_ratio(E: np.ndarray, H: np.ndarray) -> np.ndarray:
    """Eigenvalues of E^-1 H in descending order (real parts, clipped at zero)."""
    L = np.linalg.cholesky(E)
    Li = np.linalg.inv(L)
    A = Li @ H @ Li.T
    w = np.linalg.eigvalsh((A + A.T) / 2)
    return np.sort(np.maximum(w, 0.0))[::-1]


# ---------------------------------------------------------------------------
# canonical correlations
# ---------------------------------------------------------------------------

def inv_sqrt(S: np.ndarray) -> np.ndarray:
    w, V = np.linalg.eigh((S + S.T) / 2)
    return (V / np.sqrt(w)) @ V.T


def canonical(S11: np.ndarray, S22: np.ndarray, S12: np.ndarray) -> dict:
    """Canonical correlations and raw coefficients (canonical variates with unit variance) of two
    sets of variables from their covariance blocks."""
    A = inv_sqrt(S11)
    B = inv_sqrt(S22)
    U, d, Vt = np.linalg.svd(A @ S12 @ B)
    k = min(S11.shape[0], S22.shape[0])
    d = np.clip(d[:k], 0.0, 1.0)
    return {"r": d, "A1": A @ U[:, :k], "A2": B @ Vt.T[:, :k]}


# ---------------------------------------------------------------------------
# covariance matrices
# ---------------------------------------------------------------------------

def box_m(Ss: list[np.ndarray], ns: list[float]) -> dict:
    """Modified LR statistic (Box's M with n-1 weights) for equality of covariance matrices across
    groups, with its chi2 and F approximations."""
    g = len(Ss)
    p = Ss[0].shape[0]
    dfs = [n - 1.0 for n in ns]
    tot = sum(dfs)
    Sp = sum(d * S for d, S in zip(dfs, Ss)) / tot
    M = tot * np.linalg.slogdet(Sp)[1] - sum(d * np.linalg.slogdet(S)[1] for d, S in zip(dfs, Ss))
    c1 = (sum(1.0 / d for d in dfs) - 1.0 / tot) * (2.0 * p * p + 3.0 * p - 1.0) / (6.0 * (p + 1.0) * (g - 1.0))
    c2 = (sum(1.0 / d ** 2 for d in dfs) - 1.0 / tot ** 2) * (p - 1.0) * (p + 2.0) / (6.0 * (g - 1.0))
    df = p * (p + 1) * (g - 1) / 2.0
    chi2 = float(M * (1.0 - c1))
    if c2 > c1 * c1:
        df2 = (df + 2.0) / (c2 - c1 * c1)
        F = M * (1.0 - c1 - df / df2) / df
    else:
        df2 = (df + 2.0) / (c1 * c1 - c2)
        b = df2 / (1.0 - c1 + 2.0 / df2)
        F = df2 * M / (df * (b - M))
    return {"M": float(M), "chi2": chi2, "df": df, "p": float(st.chi2.sf(chi2, df)),
            "F": float(F), "df2": float(df2), "p_F": float(st.f.sf(F, df, df2))}


def lawley(R: np.ndarray, n: float) -> dict:
    """Lawley's chi2 test that all correlations are equal (compound symmetry of the correlation matrix)."""
    p = R.shape[0]
    iu = np.triu_indices(p, 1)
    rbar = float(R[iu].mean())
    rk = (R.sum(1) - 1.0) / (p - 1.0)
    gam = (p - 1.0) ** 2 * (1.0 - (1.0 - rbar) ** 2) / (p - (p - 2.0) * (1.0 - rbar) ** 2)
    chi2 = float((n - 1.0) / (1.0 - rbar) ** 2 * (((R[iu] - rbar) ** 2).sum() - gam * ((rk - rbar) ** 2).sum()))
    df = p * (p - 1) / 2.0 - 1.0
    return {"chi2": chi2, "df": df, "p": float(st.chi2.sf(chi2, df))}


def lr_diagonal(R: np.ndarray, n: float) -> dict:
    """LR test that the covariance matrix is diagonal (correlations are all zero), Bartlett-corrected."""
    p = R.shape[0]
    chi2 = float(-(n - 1.0 - (2.0 * p + 5.0) / 6.0) * np.linalg.slogdet(R)[1])
    df = p * (p - 1) / 2.0
    return {"chi2": chi2, "df": df, "p": float(st.chi2.sf(chi2, df))}


def lr_spherical(S: np.ndarray, n: float) -> dict:
    """LR test that the covariance matrix is proportional to the identity (sphericity of the
    matrix itself), with Bartlett's correction."""
    p = S.shape[0]
    ld = np.linalg.slogdet(S)[1]
    stat = ld - p * math.log(np.trace(S) / p)
    chi2 = float(-(n - 1.0 - (2.0 * p * p + p + 2.0) / (6.0 * p)) * stat)
    df = p * (p + 1) / 2.0 - 1.0
    return {"chi2": chi2, "df": df, "p": float(st.chi2.sf(chi2, df))}


def lr_equal_correlations(Rs: list[np.ndarray], ns: list[float]) -> dict:
    """Likelihood-ratio test that the correlation matrices of the groups are equal (Jennrich's
    test is used by Stata for the general case; here the LR test with the pooled correlation)."""
    g = len(Rs)
    p = Rs[0].shape[0]
    dfs = [n - 1.0 for n in ns]
    Rp = sum(d * R for d, R in zip(dfs, Rs)) / sum(dfs)
    stat = sum(dfs) * np.linalg.slogdet(Rp)[1] - sum(d * np.linalg.slogdet(R)[1] for d, R in zip(dfs, Rs))
    df = (g - 1) * p * (p - 1) / 2.0
    return {"chi2": float(stat), "df": df, "p": float(st.chi2.sf(stat, df))}


# ---------------------------------------------------------------------------
# multivariate normality
# ---------------------------------------------------------------------------

def mardia(X: np.ndarray) -> dict:
    """Mardia's skewness and kurtosis tests of multivariate normality."""
    n, p = X.shape
    Xc = X - X.mean(0)
    S = Xc.T @ Xc / n
    Si = np.linalg.inv(S)
    G = Xc @ Si @ Xc.T
    b1 = float((G ** 3).sum() / (n * n))
    b2 = float((np.diag(G) ** 2).sum() / n)
    chi_s = n * b1 / 6.0
    df_s = p * (p + 1) * (p + 2) / 6.0
    z_k = (b2 - p * (p + 2.0)) / math.sqrt(8.0 * p * (p + 2.0) / n)
    return {"b1": b1, "b2": b2, "chi2_skew": chi_s, "df_skew": df_s, "p_skew": float(st.chi2.sf(chi_s, df_s)),
            "z_kurt": z_k, "p_kurt": float(2 * st.norm.sf(abs(z_k))),
            "chi2_joint": chi_s + z_k * z_k, "p_joint": float(st.chi2.sf(chi_s + z_k * z_k, df_s + 1))}


def doornik_hansen(X: np.ndarray) -> dict:
    """Doornik-Hansen omnibus test of multivariate normality: the data are orthogonalised with the
    correlation matrix, the skewness and kurtosis of each column are transformed to normal deviates
    and their squares are added up (chi2 with 2p degrees of freedom). `U` has one row per variable:
    skewness, kurtosis and the two transformed deviates."""
    n, p = X.shape
    S = np.cov(X.T, ddof=1).reshape(p, p)
    V = np.diag(1.0 / np.sqrt(np.diag(S)))
    lam, H = np.linalg.eigh(V @ S @ V)
    Z = (X - X.mean(0)) @ V @ (H @ np.diag(lam ** -0.5) @ H.T)
    beta = 3.0 * (n * n + 27.0 * n - 70.0) * (n + 1.0) * (n + 3.0) / ((n - 2.0) * (n + 5.0) * (n + 7.0) * (n + 9.0))
    w2 = -1.0 + math.sqrt(2.0 * (beta - 1.0))
    delta = 1.0 / math.sqrt(math.log(math.sqrt(w2)))
    d = (n - 3.0) * (n + 1.0) * (n * n + 15.0 * n - 4.0)
    a = (n - 2.0) * (n + 5.0) * (n + 7.0) * (n * n + 27.0 * n - 70.0) / (6.0 * d)
    c = (n - 7.0) * (n + 5.0) * (n + 7.0) * (n * n + 2.0 * n - 5.0) / (6.0 * d)
    k = (n + 5.0) * (n + 7.0) * (n ** 3 + 37.0 * n * n + 11.0 * n - 313.0) / (12.0 * d)
    U = np.zeros((p, 4))
    for j in range(p):
        z = Z[:, j] - Z[:, j].mean()
        m2, m3, m4 = (z ** 2).mean(), (z ** 3).mean(), (z ** 4).mean()
        sb1, b2 = m3 / m2 ** 1.5, m4 / m2 ** 2
        y = sb1 * math.sqrt((w2 - 1.0) * (n + 1.0) * (n + 3.0) / (12.0 * (n - 2.0)))
        z1 = delta * math.log(y + math.sqrt(y * y + 1.0))
        al = a + sb1 ** 2 * c
        z2 = (np.cbrt((b2 - 1.0 - sb1 ** 2) * 2.0 * k / (2.0 * al)) - 1.0 + 1.0 / (9.0 * al)) * math.sqrt(9.0 * al)
        U[j] = (sb1, b2, z1, z2)
    chi2 = float((U[:, 2:] ** 2).sum())
    df = 2.0 * p
    return {"chi2": chi2, "df": df, "p": float(st.chi2.sf(chi2, df)), "U": U}


def henze_zirkler(X: np.ndarray) -> dict:
    """Henze-Zirkler test of multivariate normality (lognormal approximation of the null)."""
    n, p = X.shape
    Xc = X - X.mean(0)
    S = Xc.T @ Xc / n
    Si = np.linalg.inv(S)
    G = Xc @ Si @ Xc.T
    D = np.diag(G)
    Dij = D[:, None] + D[None, :] - 2.0 * G
    b = (1.0 / math.sqrt(2.0)) * ((2.0 * p + 1.0) / 4.0) ** (1.0 / (p + 4.0)) * n ** (1.0 / (p + 4.0))
    hz = float(np.exp(-0.5 * b * b * Dij).sum() / n
               - 2.0 * (1.0 + b * b) ** (-p / 2.0) * np.exp(-0.5 * b * b / (1.0 + b * b) * D).sum()
               + n * (1.0 + 2.0 * b * b) ** (-p / 2.0))
    b2 = b * b
    a = 1.0 + 2.0 * b2
    mu = 1.0 - a ** (-p / 2.0) * (1.0 + p * b2 / a + p * (p + 2.0) * b2 * b2 / (2.0 * a * a))
    # lognormal moments of the null distribution (Henze and Zirkler, 1990)
    w = (1.0 + b2) * (1.0 + 3.0 * b2)
    s2 = (2.0 * (1.0 + 4.0 * b2) ** (-p / 2.0)
          + 2.0 * a ** (-p) * (1.0 + 2.0 * p * b2 * b2 / (a * a) + 3.0 * p * (p + 2.0) * b2 ** 4 / (4.0 * a ** 4))
          - 4.0 * w ** (-p / 2.0) * (1.0 + 3.0 * p * b2 * b2 / (2.0 * w) + p * (p + 2.0) * b2 ** 4 / (2.0 * w * w)))
    pmu = math.log(math.sqrt(mu ** 4 / (s2 + mu * mu)))
    psig = math.sqrt(math.log((s2 + mu * mu) / (mu * mu)))
    z = (math.log(hz) - pmu) / psig if hz > 0 else float("-inf")
    return {"hz": hz, "z": z, "p": float(st.norm.sf(z))}

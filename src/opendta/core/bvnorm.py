"""Bivariate normal distribution function and its derivatives, vectorized.

    F(a, b; rho) = P(X <= a, Y <= b)   with X, Y standard normal, corr(X, Y) = rho

Plackett's identity, F = Phi(a) Phi(b) + int_0^rho phi2(a, b; r) dr, is integrated by
Gauss-Legendre quadrature (the integrand is smooth for |rho| < 1). The derivatives are exact:

    dF/da = phi(a) Phi((b - rho a) / s)        dF/db = phi(b) Phi((a - rho b) / s)
    dF/drho = phi2(a, b; rho)                  with s = sqrt(1 - rho^2)

Used by the bivariate probit models and the multinomial probit models.
"""

from __future__ import annotations

import numpy as np

from ..lazy import lazy

st = lazy("scipy.stats")

_RHO_MAX = 1 - 1e-9
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(64)


def _clip_rho(rho: np.ndarray) -> np.ndarray:
    return np.clip(rho, -_RHO_MAX, _RHO_MAX)


def bvn_pdf(a, b, rho):
    """Density of the standard bivariate normal."""
    a, b, rho = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float), _clip_rho(np.asarray(rho, float)))
    s2 = 1 - rho * rho
    return np.exp(-(a * a - 2 * rho * a * b + b * b) / (2 * s2)) / (2 * np.pi * np.sqrt(s2))


def bvn_cdf(a, b, rho):
    """P(X <= a, Y <= b) for standard normals with correlation rho (arrays of one shape)."""
    a, b, rho = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float), _clip_rho(np.asarray(rho, float)))
    r = rho[..., None] * (1 + _NODES) / 2
    s2 = 1 - r * r
    f = np.exp(-(a[..., None] ** 2 - 2 * r * a[..., None] * b[..., None] + b[..., None] ** 2) / (2 * s2)) \
        / (2 * np.pi * np.sqrt(s2))
    integral = rho / 2 * (f @ _WEIGHTS)
    out = st.norm.cdf(a) * st.norm.cdf(b) + integral
    return np.clip(out, 0.0, 1.0)


def bvn_grad(a, b, rho):
    """(F, dF/da, dF/db, dF/drho)."""
    a, b, rho = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float), _clip_rho(np.asarray(rho, float)))
    s = np.sqrt(1 - rho * rho)
    F = bvn_cdf(a, b, rho)
    Fa = st.norm.pdf(a) * st.norm.cdf((b - rho * a) / s)
    Fb = st.norm.pdf(b) * st.norm.cdf((a - rho * b) / s)
    return F, Fa, Fb, bvn_pdf(a, b, rho)

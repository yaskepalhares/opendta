"""Log likelihoods of the binary, count, censored and selection models, per observation.

Every function takes arrays (the outcome and the linear indexes of the equations, plus the
ancillary parameters on the scale where they are estimated) and returns the log likelihood of
each observation and its derivatives with respect to those indexes and parameters, in the order
of the arguments. No Session, no variables, no printing: the commands build the indexes and the
engine in `mleq.py` assembles scores and Hessians.
"""

from __future__ import annotations

import numpy as np

from ..core.bvnorm import bvn_grad
from ..lazy import lazy

special = lazy("scipy.special")
st = lazy("scipy.stats")


def _lam(x):
    """phi(x) / Phi(x), stable in the lower tail."""
    return np.exp(st.norm.logpdf(x) - st.norm.logcdf(x))


def _lam_up(x):
    """phi(x) / (1 - Phi(x))."""
    return np.exp(st.norm.logpdf(x) - st.norm.logsf(x))


# ---------------------------------------------------------------------------
# binary
# ---------------------------------------------------------------------------

def cloglog(y, eta):
    """P(y = 1) = 1 - exp(-exp(eta))."""
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        u = np.exp(np.clip(eta, -700, 700))
        pos = y > 0
        ll = np.where(pos, np.log(-np.expm1(-u)), -u)
        d = np.where(pos, u / np.expm1(u), -u)
    d = np.where(np.isfinite(d), d, 0.0)
    return ll, d


def scobit(y, eta, lnalpha):
    """P(y = 1) = 1 - (1 + exp(eta))^-alpha, alpha = exp(lnalpha) (skewed logistic)."""
    alpha = np.exp(lnalpha)
    sp = np.logaddexp(0, eta)
    u = alpha * sp
    sig = special.expit(eta)
    pos = y > 0
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        em = np.expm1(u)
        ll = np.where(pos, np.log(-np.expm1(-u)), -u)
        d_eta = np.where(pos, alpha * sig / em, -alpha * sig)
        d_la = np.where(pos, u / em, -u)
    d_eta = np.where(np.isfinite(d_eta), d_eta, 0.0)
    d_la = np.where(np.isfinite(d_la), d_la, 0.0)
    return ll, d_eta, d_la


def hetprobit(y, eta, zeta):
    """P(y = 1) = Phi(eta / exp(zeta))."""
    t = eta * np.exp(-zeta)
    q = np.where(y > 0, 1.0, -1.0)
    ll = st.norm.logcdf(q * t)
    d = q * _lam(q * t)
    return ll, d * np.exp(-zeta), -d * t


def fractional(y, eta, link):
    """Bernoulli quasi-log-likelihood for y in [0, 1]: y ln F + (1 - y) ln (1 - F)."""
    if link == "logit":
        ll = -y * np.logaddexp(0, -eta) - (1 - y) * np.logaddexp(0, eta)
        return ll, y - special.expit(eta)
    if link == "probit":
        lp, lq = st.norm.logcdf(eta), st.norm.logcdf(-eta)
        ll = y * lp + (1 - y) * lq
        return ll, y * np.exp(st.norm.logpdf(eta) - lp) - (1 - y) * np.exp(st.norm.logpdf(eta) - lq)
    # cloglog: ln(1 - F) = -u
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        u = np.exp(np.clip(eta, -700, 700))
        ll = np.where(y > 0, y * np.log(-np.expm1(-u)), 0.0) - (1 - y) * u
        d = np.where(y > 0, y * u / np.expm1(u), 0.0) - (1 - y) * u
    return ll, np.where(np.isfinite(d), d, 0.0)


# ---------------------------------------------------------------------------
# count
# ---------------------------------------------------------------------------

def _nb_terms(y, mu, m):
    """ln NB2(y; mu, alpha = 1/m) and the derivatives of it with respect to eta = ln mu and
    to ln alpha."""
    lg = special.gammaln
    ll = lg(y + m) - lg(m) - lg(y + 1) + m * np.log(m / (m + mu)) + y * np.log(mu / (m + mu))
    d_eta = (y - mu) / (1 + mu / m)
    dm = special.digamma(y + m) - special.digamma(m) + np.log(m / (m + mu)) + 1 - (m + y) / (m + mu)
    return ll, d_eta, -m * dm


def _infl_prob(zeta, link):
    """(pi, ln pi, ln(1 - pi), dpi/dzeta) of the inflation equation."""
    if link == "probit":
        return (st.norm.cdf(zeta), st.norm.logcdf(zeta), st.norm.logcdf(-zeta), st.norm.pdf(zeta))
    pi = special.expit(zeta)
    return pi, -np.logaddexp(0, -zeta), -np.logaddexp(0, zeta), pi * (1 - pi)


def zip_(y, eta, zeta, link="logit"):
    """Zero-inflated Poisson: mean exp(eta), inflation probability F(zeta)."""
    mu = np.exp(np.clip(eta, -700, 700))
    pi, lpi, lq, dpi = _infl_prob(zeta, link)
    zero = y == 0
    lp0 = np.logaddexp(lpi, lq - mu)
    ll = np.where(zero, lp0, lq - mu + y * np.log(mu) - special.gammaln(y + 1))
    p0 = np.exp(lp0)
    d_eta = np.where(zero, (1 - pi) * np.exp(-mu) * (-mu) / p0, y - mu)
    d_zeta = np.where(zero, dpi * (1 - np.exp(-mu)) / p0, -dpi / (1 - pi))
    return ll, d_eta, d_zeta


def zinb(y, eta, zeta, lnalpha, link="logit"):
    """Zero-inflated negative binomial (NB2)."""
    mu = np.exp(np.clip(eta, -700, 700))
    m = np.exp(-lnalpha)
    pi, lpi, lq, dpi = _infl_prob(zeta, link)
    zero = y == 0
    lnb, dnb_eta, dnb_la = _nb_terms(y, mu, m)
    lnb0, dnb0_eta, dnb0_la = _nb_terms(np.zeros_like(y), mu, m)
    lp0 = np.logaddexp(lpi, lq + lnb0)
    ll = np.where(zero, lp0, lq + lnb)
    p0 = np.exp(lp0)
    w0 = (1 - pi) * np.exp(lnb0) / p0
    d_eta = np.where(zero, w0 * dnb0_eta, dnb_eta)
    d_la = np.where(zero, w0 * dnb0_la, dnb_la)
    d_zeta = np.where(zero, dpi * (1 - np.exp(lnb0)) / p0, -dpi / (1 - pi))
    return ll, d_eta, d_zeta, d_la


def tpoisson(y, eta, lower):
    """Poisson truncated below: only y > lower is observed."""
    mu = np.exp(np.clip(eta, -700, 700))
    lsf = st.poisson.logsf(lower, mu)
    ll = -mu + y * np.log(mu) - special.gammaln(y + 1) - lsf
    d_eta = (y - mu) - mu * np.exp(st.poisson.logpmf(lower, mu) - lsf)
    return ll, d_eta


def tnbreg(y, eta, lnalpha, lower):
    """Negative binomial (NB2) truncated below: only y > lower is observed. `lnalpha` may be a
    scalar or an array."""
    mu = np.exp(np.clip(eta, -700, 700))
    m = np.exp(-lnalpha) * np.ones_like(mu)
    lnb, d_eta, d_la = _nb_terms(y, mu, m)
    cdf = np.zeros_like(mu)
    dc_eta = np.zeros_like(mu)
    dc_la = np.zeros_like(mu)
    for j in range(int(np.floor(lower)) + 1):
        lj, de, dl = _nb_terms(np.full_like(mu, j), mu, m)
        pj = np.exp(lj)
        cdf += pj
        dc_eta += pj * de
        dc_la += pj * dl
    sf = np.maximum(1 - cdf, 1e-300)
    return lnb - np.log(sf), d_eta + dc_eta / sf, d_la + dc_la / sf


def gnbreg(y, eta, zeta):
    """Negative binomial with ln alpha = zeta, an index of its own."""
    mu = np.exp(np.clip(eta, -700, 700))
    m = np.exp(-zeta)
    return _nb_terms(y, mu, m)


# ---------------------------------------------------------------------------
# censored, truncated, interval
# ---------------------------------------------------------------------------

def intreg(lo, hi, eta, lnsigma):
    """Interval regression: lo or hi missing (nan) means left or right censoring; lo == hi an
    exact observation."""
    n = len(eta)
    lns = lnsigma * np.ones(n)
    sg = np.exp(lns)
    ll = np.zeros(n)
    d_eta = np.zeros(n)
    d_ls = np.zeros(n)
    exact = (lo == hi)
    left = np.isnan(lo) & ~np.isnan(hi)
    right = np.isnan(hi) & ~np.isnan(lo)
    inter = ~exact & ~np.isnan(lo) & ~np.isnan(hi)
    e = (lo[exact] - eta[exact]) / sg[exact]
    ll[exact] = st.norm.logpdf(e) - lns[exact]
    d_eta[exact] = e / sg[exact]
    d_ls[exact] = e * e - 1
    z = (hi[left] - eta[left]) / sg[left]
    lam = _lam(z)
    ll[left] = st.norm.logcdf(z)
    d_eta[left] = -lam / sg[left]
    d_ls[left] = -lam * z
    z = (lo[right] - eta[right]) / sg[right]
    lam = _lam(-z)
    ll[right] = st.norm.logcdf(-z)
    d_eta[right] = lam / sg[right]
    d_ls[right] = lam * z
    zl = (lo[inter] - eta[inter]) / sg[inter]
    zu = (hi[inter] - eta[inter]) / sg[inter]
    # Phi(zu) - Phi(zl) in the tail where it is smaller
    flip = (zl + zu) > 0
    a = np.where(flip, -zu, zl)
    b = np.where(flip, -zl, zu)
    p = np.maximum(st.norm.cdf(b) - st.norm.cdf(a), 1e-300)
    ll[inter] = np.log(p)
    fa, fb = st.norm.pdf(zl), st.norm.pdf(zu)
    d_eta[inter] = (fa - fb) / sg[inter] / p
    d_ls[inter] = (zl * fa - zu * fb) / p
    return ll, d_eta, d_ls


def truncreg(y, eta, lnsigma, lower, upper):
    """Normal regression truncated to (lower, upper); infinite limits are +-inf."""
    lns = lnsigma * np.ones(len(eta))
    sg = np.exp(lns)
    e = (y - eta) / sg
    zl = (lower - eta) / sg
    zu = (upper - eta) / sg
    flip = (zl + zu) > 0
    a = np.where(flip, -zu, zl)
    b = np.where(flip, -zl, zu)
    p = np.maximum(st.norm.cdf(b) - st.norm.cdf(a), 1e-300)
    fl = np.where(np.isfinite(zl), st.norm.pdf(np.where(np.isfinite(zl), zl, 0)), 0.0)
    fu = np.where(np.isfinite(zu), st.norm.pdf(np.where(np.isfinite(zu), zu, 0)), 0.0)
    zlf = np.where(np.isfinite(zl), zl, 0.0)
    zuf = np.where(np.isfinite(zu), zu, 0.0)
    ll = st.norm.logpdf(e) - lns - np.log(p)
    d_eta = e / sg - (fl - fu) / sg / p
    d_ls = e * e - 1 - (zlf * fl - zuf * fu) / p
    return ll, d_eta, d_ls


# ---------------------------------------------------------------------------
# selection and bivariate probit
# ---------------------------------------------------------------------------

def heckman(sel, y, eta, zeta, athrho, lnsigma):
    """Normal regression observed when sel = 1; selection probit with index zeta."""
    rho = np.tanh(athrho) * np.ones(len(eta))
    lns = lnsigma * np.ones(len(eta))
    sg = np.exp(lns)
    r = np.sqrt(1 - rho * rho)
    obs = sel > 0
    e = np.where(obs, (np.where(obs, y, 0.0) - eta) / sg, 0.0)
    a = (zeta + rho * e) / r
    lam = _lam(a)
    ll = np.where(obs, st.norm.logpdf(e) - lns + st.norm.logcdf(a), st.norm.logcdf(-zeta))
    lam0 = _lam(-zeta)
    d_eta = np.where(obs, e / sg - lam * rho / (sg * r), 0.0)
    d_zeta = np.where(obs, lam / r, -lam0)
    d_ls = np.where(obs, e * e - 1 - lam * rho * e / r, 0.0)
    d_ath = np.where(obs, lam * (e + rho * zeta) / r, 0.0)
    return ll, d_eta, d_zeta, d_ath, d_ls


def etreg(t, y, eta, zeta, athrho, lnsigma):
    """Normal outcome with an endogenous binary treatment: y = eta + e, t = 1(zeta + u > 0),
    e ~ N(0, sigma^2), corr(e, u) = rho. Returns ll and the derivatives with respect to eta,
    zeta, athrho and lnsigma."""
    n = len(y)
    rho = np.tanh(athrho) * np.ones(n)
    lns = lnsigma * np.ones(n)
    sg = np.exp(lns)
    c = np.sqrt(1 - rho * rho)
    q = np.where(t > 0, 1.0, -1.0)
    u = (y - eta) / sg
    a = q * (zeta + rho * u) / c
    lam = _lam(a)
    ll = st.norm.logpdf(u) - lns + st.norm.logcdf(a)
    d_eta = (u - lam * q * rho / c) / sg
    d_zeta = lam * q / c
    d_ls = u * u - 1 - u * lam * q * rho / c
    da_drho = q * (u / c + (zeta + rho * u) * rho / c ** 3)
    d_ath = lam * da_drho * c * c
    return ll, d_eta, d_zeta, d_ath, d_ls


def etpois(t, y, eta, zeta, athrho, lnsigma, nodes=32):
    """Poisson outcome with an endogenous binary treatment: y ~ Poisson(exp(eta + e)), e ~ N(0, sigma^2),
    t = 1(zeta + u > 0), corr(e, u) = rho. The integral over e is a Gauss-Hermite quadrature.
    Returns ll and the derivatives with respect to eta, zeta, athrho and lnsigma."""
    from scipy.special import gammaln, logsumexp
    x, wq = np.polynomial.hermite.hermgauss(nodes)
    sk = np.sqrt(2.0) * x
    lw = np.log(wq / np.sqrt(np.pi))
    rho = np.tanh(athrho)
    c = np.sqrt(1 - rho * rho)
    sg = np.exp(lnsigma)
    q = np.where(t > 0, 1.0, -1.0)[:, None]
    mu = np.exp(np.clip(eta[:, None] + sg * sk[None, :], -50, 50))
    a = q * (zeta[:, None] + rho * sk[None, :]) / c
    lf = y[:, None] * np.log(mu) - mu - gammaln(y + 1)[:, None] + st.norm.logcdf(a) + lw[None, :]
    ll = logsumexp(lf, axis=1)
    pi = np.exp(lf - ll[:, None])
    lam = np.exp(st.norm.logpdf(a) - st.norm.logcdf(a))
    r = y[:, None] - mu
    d_eta = (pi * r).sum(axis=1)
    d_zeta = (pi * lam * q / c).sum(axis=1)
    d_ls = (pi * r * sg * sk[None, :]).sum(axis=1)
    da = q * (sk[None, :] / c + (zeta[:, None] + rho * sk[None, :]) * rho / c ** 3)
    d_ath = (pi * lam * da).sum(axis=1) * c * c
    return ll, d_eta, d_zeta, d_ath, d_ls


def biprobit(y1, y2, eta1, eta2, athrho):
    """Bivariate probit: P(y1, y2) = Phi2(q1 eta1, q2 eta2; q1 q2 rho)."""
    rho = np.tanh(athrho)
    q1 = np.where(y1 > 0, 1.0, -1.0)
    q2 = np.where(y2 > 0, 1.0, -1.0)
    F, Fa, Fb, Fr = bvn_grad(q1 * eta1, q2 * eta2, q1 * q2 * rho)
    F = np.maximum(F, 1e-300)
    return (np.log(F), q1 * Fa / F, q2 * Fb / F, q1 * q2 * Fr * (1 - rho * rho) / F)


def heckprobit(sel, y, eta, zeta, athrho):
    """Probit for y observed when sel = 1 (selection probit with index zeta)."""
    rho = np.tanh(athrho)
    obs = sel > 0
    q = np.where(y > 0, 1.0, -1.0)
    # selected: P = Phi2(zeta, q eta; q rho)
    F, Fa, Fb, Fr = bvn_grad(zeta, q * eta, q * rho)
    F = np.maximum(F, 1e-300)
    ll = np.where(obs, np.log(F), st.norm.logcdf(-zeta))
    d_eta = np.where(obs, q * Fb / F, 0.0)
    d_zeta = np.where(obs, Fa / F, -_lam(-zeta))
    d_ath = np.where(obs, q * Fr * (1 - rho * rho) / F, 0.0)
    return ll, d_eta, d_zeta, d_ath


def heckoprobit(sel, yi, ncat, eta, zeta, cuts, athrho):
    """Ordered probit for y (categories 0..ncat-1) observed when sel = 1. Returns the log
    likelihood and the derivatives with respect to eta, zeta, athrho and the cutpoints (n x
    (ncat - 1))."""
    rho = np.tanh(athrho)
    obs = sel > 0
    n = len(eta)
    k = ncat - 1
    big = 1e10
    cu = np.concatenate([[-big], cuts, [big]])
    hi = cu[np.where(obs, yi, 0) + 1]
    lo = cu[np.where(obs, yi, 0)]
    # P(s = 1, y = j) = Phi2(zeta, hi - eta; -rho) - Phi2(zeta, lo - eta; -rho)
    Fh, Fha, Fhb, Fhr = bvn_grad(zeta, hi - eta, -rho)
    Fl, Fla, Flb, Flr = bvn_grad(zeta, lo - eta, -rho)
    P = np.maximum(Fh - Fl, 1e-300)
    ll = np.where(obs, np.log(P), st.norm.logcdf(-zeta))
    d_zeta = np.where(obs, (Fha - Fla) / P, -_lam(-zeta))
    d_eta = np.where(obs, -(Fhb - Flb) / P, 0.0)
    d_ath = np.where(obs, -(Fhr - Flr) * (1 - rho * rho) / P, 0.0)
    d_cut = np.zeros((n, k))
    for j in range(k):
        up = obs & (yi == j)
        dn = obs & (yi == j + 1)
        d_cut[:, j] = np.where(up, Fhb / P, 0.0) - np.where(dn, Flb / P, 0.0)
    return ll, d_eta, d_zeta, d_ath, d_cut


# ---------------------------------------------------------------------------
# endogenous regressor: ivprobit, ivtobit (y2 = eta2 + v, v ~ N(0, sigma_v^2), corr(u, v) = rho)
# ---------------------------------------------------------------------------

def ivprobit(y, y2, eta1, eta2, athrho, lnsigma):
    """ll_obs and derivatives wrt eta1, eta2, athrho and lnsigma (sigma = sd of y2's error);
    the probit error has variance one."""
    n = len(y)
    r = np.tanh(athrho) * np.ones(n)
    c = np.sqrt(1 - r * r)
    sv = np.exp(lnsigma) * np.ones(n)
    v = y2 - eta2
    u = v / sv
    q = 2 * y - 1
    t = (eta1 + r * u) / c
    z = q * t
    lcdf = np.where(z > -30, special.log_ndtr(z), -0.5 * z * z - np.log(-z) - 0.5 * np.log(2 * np.pi))
    ll = -0.5 * u * u - lnsigma - 0.5 * np.log(2 * np.pi) + lcdf
    lam = q * np.exp(-0.5 * z * z - 0.5 * np.log(2 * np.pi) - lcdf)
    d1 = lam / c
    d2 = u / sv - lam * r / (sv * c)
    da = lam * (c * u + t * r)
    dls = u * u - 1 - lam * r * u / c
    return ll, d1, d2, da, dls


def ivtobit(y, y2, eta1, eta2, alpha, lns, lnv, lower, upper):
    """Tobit with an endogenous regressor, in Stata's parameters: y* = eta1 + alpha*v + e with
    v = y2 - eta2 ~ N(0, exp(lnv)^2) and e ~ N(0, exp(lns)^2). Derivatives wrt eta1, eta2, alpha,
    lns and lnv; lower/upper are the censoring limits (+-inf when absent)."""
    n = len(y)
    s_ = np.exp(lns) * np.ones(n)
    sv = np.exp(lnv) * np.ones(n)
    v = y2 - eta2
    u = v / sv
    m = eta1 + alpha * v
    left = y <= lower
    right = y >= upper
    unc = ~(left | right)
    zz = np.where(left, (lower - m) / s_, np.where(right, (m - upper) / s_, (y - m) / s_))
    lcdf = np.where(zz > -30, special.log_ndtr(zz), -0.5 * zz * zz - np.log(np.abs(zz) + 1e-300) - 0.5 * np.log(2 * np.pi))
    h = np.exp(-0.5 * zz * zz - 0.5 * np.log(2 * np.pi) - lcdf)
    llc = np.where(unc, -0.5 * zz * zz - np.log(s_) - 0.5 * np.log(2 * np.pi), lcdf)
    Dm = np.where(unc, zz / s_, np.where(left, -h / s_, h / s_))
    DS = np.where(unc, (zz * zz - 1) / s_, -h * zz / s_)
    ll = llc - 0.5 * u * u - lnv - 0.5 * np.log(2 * np.pi)
    d1 = Dm
    d2 = -alpha * Dm + u / sv
    da = Dm * v
    dls = DS * s_
    dlv = u * u - 1
    return ll, d1, d2, da, dls, dlv

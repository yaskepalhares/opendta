"""Mixed-effects generalized linear models: meglm and its special cases melogit, meprobit,
mecloglog, mepoisson, menbreg, meologit and meoprobit ([ME] meglm and the others).

Two-level models with random intercepts and slopes. The likelihood of
each group is integrated by mean-variance adaptive Gauss-Hermite
quadrature (intmethod(mvaghermite), intpoints(7) by default; with q
random effects, the product rule with 7^q points); intmethod(ghermite)
gives the non-adaptive rule and intmethod(laplace) a single point. The
adaptive parameters (posterior means and Cholesky factors) are updated at
each iteration and held fixed within it; the gradient is analytic and the
Hessian is the numerical derivative of the gradient.

Stata's starting values ("Refining starting values" with a grid search)
are not documented, so the iteration logs may differ; the final
estimates should agree to the precision of the quadrature (VERIFY:
compat/do/1002_melogit.do).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..core.stats import chi2_tail
from ..lang.syntax import find_top
from ..lazy import lazy
from . import fvars
from .mixed import (
    Level,
    _p_chibar,
    _right,
    build_levels,
    estat_group,
    estat_icc,
    estat_recov,
    group_header,
    lr_line,
    split_levels,
)
from .ml import maximize, print_log
from .ml import options as ml_options
from .results import CoefRow, Estimates, coef_table, post, rows_from_columns

# imported on first use (start-up time; see opendta/lazy.py)
special = lazy("scipy.special")
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _g(x: float, fmt: str = "%9.0g") -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


# ---------------------------------------------------------------------------
# families
# ---------------------------------------------------------------------------

def _family(kind: str, y: np.ndarray):
    """(ln f(y|eta), d ln f / d eta, d2 ln f / d eta2) vectorized over eta (any shape,
    with y broadcast on the first axis). The optional second argument holds the auxiliary
    parameters of the family (unused by logit, probit and poisson)."""
    yc = y[:, None]
    if kind == "logit":
        def f(eta, a=None):
            lf = yc * eta - np.logaddexp(0, eta)
            p = special.expit(eta)
            return lf, yc - p, -p * (1 - p)
    elif kind == "probit":
        q = 2 * yc - 1

        def f(eta, a=None):
            lc = special.log_ndtr(q * eta)
            lam = q * np.exp(st.norm.logpdf(q * eta) - lc)
            return lc, lam, -lam * (lam + eta)
    else:
        lgy = special.gammaln(yc + 1)

        def f(eta, a=None):
            mu = np.exp(np.clip(eta, -700, 700))
            return yc * eta - mu - lgy, yc - mu, -mu
    return f


def family_of(fam: str, link: str, y: np.ndarray, *, nbin: np.ndarray | None = None,
              disp: str = "mean"):
    """The conditional distribution of the other families of meglm: f(eta, a) with the
    auxiliary parameters a (ln variance of the gaussian, ln scale of the gamma, ln dispersion
    of the negative binomial, the cutpoints of the ordinal models) and the number of them."""
    yc = y[:, None]
    ln2pi = math.log(2 * math.pi)
    if fam == "gaussian" and link == "log":
        def f(eta, a):
            v = math.exp(a[0])
            mu = np.exp(np.clip(eta, -700, 700))
            r = yc - mu
            return -0.5 * (ln2pi + a[0]) - r * r / (2 * v), r * mu / v, (r * mu - mu * mu) / v
        return f, 1
    if fam == "gaussian":
        def f(eta, a):
            v = math.exp(a[0])
            r = yc - eta
            return -0.5 * (ln2pi + a[0]) - r * r / (2 * v), r / v, np.full(eta.shape, -1.0 / v)
        return f, 1
    if fam == "gamma":
        lny = np.log(yc)

        def f(eta, a):
            lphi = 2.0 * a[0]            # /logs is the log of the square root of the scale (compat 1902)
            phi = math.exp(lphi)
            k = 1.0 / phi
            r = yc * np.exp(-np.clip(eta, -700, 700)) / phi
            return (k - 1) * lny - r - k * (lphi + eta) - special.gammaln(k), r - k, -r
        return f, 1
    if fam == "nbinomial":
        lgy = special.gammaln(yc + 1)
        if disp == "mean":
            def f(eta, a):
                al = math.exp(a[0])
                ia = 1.0 / al
                mu = np.exp(np.clip(eta, -700, 700))
                lf = (special.gammaln(yc + ia) - special.gammaln(ia) - lgy + yc * (a[0] + eta)
                      - (yc + ia) * np.log1p(al * mu))
                return lf, (yc - mu) / (1 + al * mu), -mu * (1 + al * yc) / (1 + al * mu) ** 2
        else:
            def f(eta, a):
                dl = math.exp(a[0])
                mu = np.exp(np.clip(eta, -700, 700))
                rr = mu / dl
                lf = (special.gammaln(yc + rr) - special.gammaln(rr) - lgy - rr * math.log1p(dl)
                      + yc * (a[0] - math.log1p(dl)))
                s1 = rr * (special.digamma(yc + rr) - special.digamma(rr) - math.log1p(dl))
                return lf, s1, s1 + rr * rr * (special.polygamma(1, yc + rr) - special.polygamma(1, rr))
        return f, 1
    if fam == "binomial":
        nc = nbin[:, None]
        lch = special.gammaln(nc + 1) - special.gammaln(yc + 1) - special.gammaln(nc - yc + 1)
        if link == "logit":
            def f(eta, a=None):
                p = special.expit(eta)
                return yc * eta - nc * np.logaddexp(0, eta) + lch, yc - nc * p, -nc * p * (1 - p)
        elif link == "probit":
            def f(eta, a=None):
                lp, lq = special.log_ndtr(eta), special.log_ndtr(-eta)
                l1 = np.exp(st.norm.logpdf(eta) - lp)
                l0 = np.exp(st.norm.logpdf(eta) - lq)
                return (yc * lp + (nc - yc) * lq + lch, yc * l1 - (nc - yc) * l0,
                        -yc * l1 * (l1 + eta) - (nc - yc) * l0 * (l0 - eta))
        else:
            f = _cloglog_binomial(yc, nc, lch)
        return f, 0
    raise StataError(198, f"family({fam}) is not available yet in OpenDTA")   # OpenDTA


def _cloglog_binomial(yc, nc, lch):
    """p = 1 - exp(-exp(eta)); with m = exp(eta), d ln p / d eta = m / (e^m - 1) = h and
    d2 ln p / d eta2 = h - (m/2 / sinh(m/2))^2, d ln(1-p) / d eta = -m."""
    def f(eta, a=None):
        m = np.exp(np.clip(eta, -700, 6.5))
        h = m / np.expm1(m)
        w = (m / 2) / np.sinh(m / 2)
        lp = np.log(-np.expm1(-m))
        return (yc * lp - (nc - yc) * m + lch, yc * h - (nc - yc) * m,
                yc * (h - w * w) - (nc - yc) * m)
    return f


def ordinal_family(link: str, yi: np.ndarray):
    """Ordinal outcome 0..J-1 with P(y = j) = F(k_j - eta) - F(k_(j-1) - eta); a holds the J-1
    cutpoints (in increasing order)."""
    if link == "logit":
        F = special.expit
        pdf = lambda z: F(z) * (1 - F(z))                       # noqa: E731
        dpdf = lambda z: pdf(z) * (1 - 2 * F(z))                 # noqa: E731
    else:
        F = special.ndtr
        pdf = st.norm.pdf
        dpdf = lambda z: np.where(np.isfinite(z), -np.where(np.isfinite(z), z, 0) * st.norm.pdf(z), 0.0)   # noqa: E731

    def f(eta, a):
        kk = np.concatenate([[-np.inf], a, [np.inf]])
        up = kk[yi + 1][:, None] - eta
        lo = kk[yi][:, None] - eta
        p = np.maximum(F(up) - F(lo), 1e-300)
        s1 = -(pdf(up) - pdf(lo)) / p
        return np.log(p), s1, (dpdf(up) - dpdf(lo)) / p - s1 * s1
    return f, int(yi.max())


def _gh_nodes(n: int, q: int):
    """Product rule for N(0, I_q): nodes (K x q) and log weights (K)."""
    x, w = np.polynomial.hermite_e.hermegauss(n)
    w = w / math.sqrt(2 * math.pi)
    grids = np.meshgrid(*([x] * q), indexing="ij")
    nodes = np.column_stack([g.ravel() for g in grids])
    wg = np.meshgrid(*([w] * q), indexing="ij")
    lw = np.sum([np.log(g.ravel()) for g in wg], axis=0)
    return nodes, lw


def _chol_rows(cov, old):
    """Cholesky factors of a stack of covariance matrices; a matrix that is not numerically
    positive definite keeps its old factor."""
    try:
        return np.linalg.cholesky(cov)
    except np.linalg.LinAlgError:
        out = old.copy()
        for i in range(cov.shape[0]):
            try:
                out[i] = np.linalg.cholesky(cov[i])
            except np.linalg.LinAlgError:
                pass
        return out


class GLMM:
    def __init__(self, kind, y, X, offset, level: Level, npoints: int, method: str, fam=None,
                 naux: int = 0, ordered: bool = False, wobs=None, wgrp=None):
        self.kind = kind
        self.wobs = wobs                  # observation weights (None: all 1)
        self.wgrp = wgrp                  # group weights (None: all 1)
        self.naux = naux                  # auxiliary parameters between beta and the variances
        self.ordered = ordered            # the auxiliary parameters are cutpoints
        self.y, self.X, self.off = y, X, offset
        self.lv = level
        self.q = level.q
        self.gid = level.gid
        self.G = int(self.gid.max()) + 1
        self.Z = level.Z
        self.fam = fam if fam is not None else _family(kind, y)
        if method == "laplace":
            npoints = 1
        self.nodes, self.lw = _gh_nodes(npoints, self.q)
        self.lphi = st.norm.logpdf(self.nodes).sum(axis=1)       # ln phi(t_k)
        self.adaptive = method in ("mvaghermite", "mcaghermite", "laplace")
        self.mu = np.zeros((self.G, self.q))
        self.C = None                                           # G x q x q
        self.npar_th = level.npar()

    @property
    def mus(self):
        return [self.mu]

    @property
    def gtop(self):
        return self.gid

    def reset(self):
        self.C = None
        self.mu[:] = 0

    def _wg(self):
        return np.ones(self.G) if self.wgrp is None else self.wgrp

    def split(self, b):
        p = self.X.shape[1]
        return b[:p], b[p + self.naux:]

    def aux(self, b):
        p = self.X.shape[1]
        return b[p:p + self.naux]

    def _u(self):
        """Quadrature points of each group: G x K x q."""
        return self.mu[:, None, :] + np.einsum("gij,kj->gki", self.C, self.nodes)

    def _terms(self, b):
        beta, th = self.split(b)
        S = self.lv.sigma(th)
        try:
            Sinv = np.linalg.inv(S)
            _s, ldS = np.linalg.slogdet(S)
        except np.linalg.LinAlgError:
            return None
        if self.C is None:
            L = np.linalg.cholesky(S + 1e-12 * np.eye(self.q))
            self.C = np.broadcast_to(L, (self.G, self.q, self.q)).copy()
        if not self.adaptive:
            L = np.linalg.cholesky(S)
            self.C = np.broadcast_to(L, (self.G, self.q, self.q)).copy()
            self.mu[:] = 0
        U = self._u()                                           # G x K x q
        Uo = U[self.gid]                                        # N x K x q
        eta = (self.X @ beta + self.off)[:, None] + np.einsum("nq,nkq->nk", self.Z, Uo)
        lf, s1, s2 = self.fam(eta, self.aux(b))
        if self.wobs is not None:
            w_ = self.wobs[:, None]
            lf, s1, s2 = lf * w_, s1 * w_, s2 * w_
        K = self.nodes.shape[0]
        sumlf = np.zeros((self.G, K))
        np.add.at(sumlf, self.gid, lf)
        quad = np.einsum("gki,ij,gkj->gk", U, Sinv, U)
        lprior = -0.5 * (self.q * math.log(2 * math.pi) + ldS + quad)
        _s2, ldC = np.linalg.slogdet(self.C)
        lg = sumlf + lprior - self.lphi[None, :] + ldC[:, None] + self.lw[None, :]
        lL = special.logsumexp(lg, axis=1)
        return lL, lg, s1, s2, U, Sinv, S, eta

    def loglik(self, b) -> float:
        if self.ordered and np.any(np.diff(self.aux(b)) <= 0):
            return -np.inf
        r = self._terms(b)
        if r is None:
            return -np.inf
        ll = float(r[0] @ self._wg())
        return ll if np.isfinite(ll) else -np.inf

    def contrib(self, b) -> np.ndarray:
        """Weighted log likelihood of each group (the adaptive parameters held fixed, except
        with the Laplace approximation, whose mode moves with b)."""
        if self.nodes.shape[0] == 1:
            self._mode_curvature(b)
        r = self._terms(b)
        return r[0] * self._wg()

    def scores(self, b) -> np.ndarray:
        """G x k: derivatives of the weighted log likelihood of each group (central differences;
        for the robust variance)."""
        k = b.size
        Sc = np.zeros((self.G, k))
        mu0, C0 = self.mu.copy(), self.C.copy()
        for i in range(k):
            h = 1e-5 * (abs(b[i]) + 1)
            e = np.zeros(k)
            e[i] = h
            up = self.contrib(b + e)
            self.mu, self.C = mu0.copy(), C0.copy()
            dn = self.contrib(b - e)
            self.mu, self.C = mu0.copy(), C0.copy()
            Sc[:, i] = (up - dn) / (2 * h)
        return Sc

    def adapt(self, b, steps: int = 20) -> None:
        """Posterior means and covariances of the random effects by quadrature, iterated."""
        if not self.adaptive:
            return
        for _ in range(steps):
            r = self._terms(b)
            if r is None:
                return
            lL, lg, *_rest = r
            U = r[4]
            pw = np.exp(lg - lL[:, None])                         # G x K
            mu = np.einsum("gk,gkq->gq", pw, U)
            dev = U - mu[:, None, :]
            cov = np.einsum("gk,gki,gkj->gij", pw, dev, dev)
            if self.nodes.shape[0] == 1:
                return self._mode_curvature(b)
            fin = np.all(np.isfinite(cov), axis=(1, 2))
            ok = fin & np.all(np.linalg.eigvalsh(np.where(fin[:, None, None], cov, 0)) > 1e-12, axis=1)
            newC = self.C.copy()
            newC[ok] = _chol_rows(cov[ok], newC[ok])
            if not ok.all():
                # points too wide for the prior all but one fall in its tails and the
                # posterior variance comes out zero: restart those groups from the prior
                newC[~ok] = np.linalg.cholesky(r[6] + 1e-300 * np.eye(self.q))
            change = np.max(np.abs(mu - self.mu)) + np.max(np.abs(newC - self.C))
            self.mu = mu
            self.C = newC
            if change < 1e-8:
                break

    def _mode_curvature(self, b):
        """Laplace: mode and curvature of each group's integrand (Newton)."""
        beta, th = self.split(b)
        S = self.lv.sigma(th)
        Sinv = np.linalg.inv(S)
        mu = self.mu.copy()
        for _ in range(50):
            eta = self.X @ beta + self.off + np.einsum("nq,nq->n", self.Z, mu[self.gid])
            _lf, s1, s2 = self.fam(eta[:, None], self.aux(b))
            if self.wobs is not None:
                s1, s2 = s1 * self.wobs[:, None], s2 * self.wobs[:, None]
            g = np.zeros((self.G, self.q))
            np.add.at(g, self.gid, self.Z * s1[:, 0][:, None])
            g -= mu @ Sinv
            H = np.zeros((self.G, self.q, self.q))
            np.add.at(H, self.gid, self.Z[:, :, None] * self.Z[:, None, :] * s2[:, 0][:, None, None])
            H -= Sinv
            step = np.linalg.solve(H, g[..., None])[..., 0]
            mu = mu - step
            if np.max(np.abs(step)) < 1e-10:
                break
        self.mu = mu
        self.C = np.linalg.cholesky(np.linalg.inv(-H))

    def gradient(self, b) -> np.ndarray:
        r = self._terms(b)
        lL, lg, s1, s2, U, Sinv, S, eta = r
        pw = np.exp(lg - lL[:, None]) * self._wg()[:, None]
        pwo = pw[self.gid]                                      # N x K
        gb = self.X.T @ np.sum(pwo * s1, axis=1)
        ga = np.zeros(self.naux)
        a = self.aux(b)
        wo = 1.0 if self.wobs is None else self.wobs[:, None]
        for m in range(self.naux):
            h = 1e-5 * (abs(a[m]) + 1)
            e = np.zeros(self.naux)
            e[m] = h
            ga[m] = float(np.sum(pwo * wo * (self.fam(eta, a + e)[0] - self.fam(eta, a - e)[0])) / (2 * h))
        # variance parameters: d/dth ln phi_S(u) = -1/2 tr(S^-1 dS) + 1/2 u' S^-1 dS S^-1 u
        beta, th = self.split(b)
        gth = np.zeros(th.size)
        h = 1e-6
        for m in range(th.size):
            e = np.zeros(th.size)
            e[m] = h
            if not self.adaptive:
                # non-adaptive points u = C t move with C = chol(S): the prior terms cancel and
                # only the linear predictor depends on theta, through dC
                dC = (np.linalg.cholesky(self.lv.sigma(th + e)) - np.linalg.cholesky(self.lv.sigma(th - e))) / (2 * h)
                dU = np.einsum("ij,kj->ki", dC, self.nodes)                    # K x q
                dEta = np.einsum("nq,kq->nk", self.Z, dU)                      # N x K
                contrib = np.zeros((self.G, dEta.shape[1]))
                np.add.at(contrib, self.gid, s1 * dEta)
                gth[m] = float(np.sum(pw * contrib))
                continue
            dS = (self.lv.sigma(th + e) - self.lv.sigma(th - e)) / (2 * h)
            A = Sinv @ dS @ Sinv
            quad = np.einsum("gki,ij,gkj->gk", U, A, U)
            gth[m] = float(np.sum(pw * (-0.5 * np.trace(Sinv @ dS) + 0.5 * quad)))
        return np.concatenate([gb, ga, gth])

    def _laplace_ll(self, b) -> float:
        self._mode_curvature(b)
        return self.loglik(b)

    def _laplace_fun(self, b, todo):
        """The Laplace approximation depends on the parameters through the mode and the curvature
        of each group, so its derivatives are numerical (the analytic gradient of the other
        methods holds the quadrature points fixed)."""
        ll = self._laplace_ll(b)
        if todo == 0 or not np.isfinite(ll):
            return ll, None, None
        k = b.size

        def grad(x):
            g = np.zeros(k)
            for i in range(k):
                h = 1e-5 * (abs(x[i]) + 1)
                e = np.zeros(k)
                e[i] = h
                g[i] = (self._laplace_ll(x + e) - self._laplace_ll(x - e)) / (2 * h)
            return g
        g = grad(b)
        H = np.zeros((k, k))
        for i in range(k):
            h = 1e-4 * (abs(b[i]) + 1)
            e = np.zeros(k)
            e[i] = h
            H[:, i] = (grad(b + e) - grad(b - e)) / (2 * h)
        self._mode_curvature(b)
        return ll, g, (H + H.T) / 2

    def fun(self, b, todo):
        if self.nodes.shape[0] == 1:
            return self._laplace_fun(b, todo)
        # the quadrature points follow b also in the line search: points adapted to a larger
        # variance overstate the likelihood of a smaller one (compat 1916)
        self.adapt(b)
        ll = self.loglik(b)
        if todo == 0 or not np.isfinite(ll):
            return (ll if np.isfinite(ll) else -np.inf), None, None
        if not self.adaptive:
            g = self.gradient(b)
            k = b.size
            H = np.zeros((k, k))
            for i in range(k):
                hi = 1e-5 * (abs(b[i]) + 1)
                e = np.zeros(k)
                e[i] = hi
                H[:, i] = (self.gradient(b + e) - self.gradient(b - e)) / (2 * hi)
            return ll, g, (H + H.T) / 2
        # derivatives of the likelihood whose points follow the parameters: the estimates and
        # standard errors of Stata (compat 1916)
        state = (self.mu.copy(), self.C.copy())

        def f(x):
            self.mu, self.C = state[0].copy(), state[1].copy()
            self.adapt(x)
            return self.loglik(x)
        k = b.size
        hs = 1e-5 * (np.abs(b) + 1)
        g = np.zeros(k)
        for i in range(k):
            e = np.zeros(k)
            e[i] = hs[i]
            g[i] = (f(b + e) - f(b - e)) / (2 * hs[i])
        hh = 1e-4 * (np.abs(b) + 1)
        H = np.zeros((k, k))
        for i in range(k):
            for j in range(i, k):
                ei = np.zeros(k)
                ej = np.zeros(k)
                ei[i] += hh[i]
                ej[j] += hh[j]
                H[i, j] = H[j, i] = (f(b + ei + ej) - f(b + ei - ej) - f(b - ei + ej) + f(b - ei - ej)) \
                    / (4 * hh[i] * hh[j])
        f(b)
        return ll, g, H


class MGLMM:
    """Nested random effects at two or more levels (the first level is the highest). The
    likelihood of a highest-level group is the nested integral over the effects of each level,
    computed by product quadrature: every level has its own points, adapted to the posterior
    mean and variance of the effects of each of its groups (computed with the joint posterior
    weights of the levels above), held fixed within an iteration. The derivatives are numerical
    (central differences with the adaptive points fixed). VERIFY: Stata adapts the points of a
    lower level to each point of the levels above; the estimates agree as the points increase."""

    def __init__(self, X, offset, levels, npoints: int, method: str, fam, naux: int = 0,
                 ordered: bool = False, wobs=None, wgrp=None):
        self.X, self.off, self.lvs = X, offset, levels
        self.fam, self.naux, self.ordered = fam, naux, ordered
        self.m = len(levels)
        self.N = X.shape[0]
        self.G = [int(lv.gid.max()) + 1 for lv in levels]
        self.parent = [None]
        for li in range(1, self.m):
            par = np.zeros(self.G[li], dtype=np.int64)
            par[levels[li].gid] = levels[li - 1].gid
            self.parent.append(par)
        self.nodes, self.lw, self.lphi = [], [], []
        for lv in levels:
            nd, lw = _gh_nodes(npoints, lv.q)
            self.nodes.append(nd)
            self.lw.append(lw)
            self.lphi.append(st.norm.logpdf(nd).sum(axis=1))
        self.adaptive = method in ("mvaghermite", "mcaghermite")
        self.wobs = wobs
        self.wgrp = wgrp if wgrp is not None else [None] * self.m
        self.npar = [lv.npar() for lv in levels]
        self.reset()

    @property
    def gtop(self):
        return self.lvs[0].gid

    def reset(self):
        self.mus = [np.zeros((g, lv.q)) for g, lv in zip(self.G, self.lvs)]
        self.Cs = None
        self.Bs = [None] + [np.zeros((self.G[li], self.lvs[li].q, self.lvs[li - 1].q)) for li in range(1, self.m)]

    # the attributes the single-level model exposes to _fit_me
    @property
    def C(self):
        return self.Cs

    @C.setter
    def C(self, v):
        self.Cs = v

    @property
    def mu(self):
        return self.mus

    @mu.setter
    def mu(self, v):
        self.mus = v

    def split(self, b):
        p = self.X.shape[1]
        return b[:p], b[p + self.naux:]

    def aux(self, b):
        p = self.X.shape[1]
        return b[p:p + self.naux]

    def _pieces(self, th):
        out, k = [], 0
        for n in self.npar:
            out.append(th[k:k + n])
            k += n
        return out

    def _points(self, li, U_par):
        """Quadrature points of level li: G_li x K_0 x ... x K_li x q. The points of a group
        move with the effects of its parent group, by the regression of the posterior (the
        conditional posterior of a Gaussian model)."""
        base = self.mus[li][:, None, :] + np.einsum("gij,kj->gki", self.Cs[li], self.nodes[li])
        if li == 0:
            return base
        par = self.parent[li]
        dev = U_par[par] - self.mus[li - 1][par].reshape((-1,) + (1,) * li + (U_par.shape[-1],))
        shift = np.einsum("gij,g...j->g...i", self.Bs[li], dev)
        G, K, q = base.shape
        return shift[..., None, :] + base.reshape((G,) + (1,) * li + (K, q))

    def _eval(self, b, keep: bool = False):
        beta, th = self.split(b)
        Ss = [lv.sigma(t) for lv, t in zip(self.lvs, self._pieces(th))]
        try:
            Sinvs = [np.linalg.inv(S) for S in Ss]
            ldSs = [np.linalg.slogdet(S)[1] for S in Ss]
            chols = [np.linalg.cholesky(S + 1e-12 * np.eye(S.shape[0])) for S in Ss]
        except np.linalg.LinAlgError:
            return None
        if self.Cs is None or not self.adaptive:
            if not self.adaptive:
                self.reset()
            self.Cs = [np.broadcast_to(L, (g, L.shape[0], L.shape[0])).copy() for L, g in zip(chols, self.G)]
        m, N = self.m, self.N
        eta = (self.X @ beta + self.off).reshape((N,) + (1,) * m)
        Us = []
        for li, lv in enumerate(self.lvs):
            U = self._points(li, Us[-1] if li else None)
            Us.append(U)
            E = np.einsum("nq,n...q->n...", lv.Z, U[lv.gid])
            eta = eta + E.reshape(E.shape + (1,) * (m - 1 - li))
        full = eta.shape
        lf = self.fam(eta.reshape(N, -1), self.aux(b))[0].reshape(full)
        if self.wobs is not None:
            lf = lf * self.wobs.reshape((N,) + (1,) * m)
        cur, idx = lf, self.lvs[-1].gid
        store = [None] * m
        for li in reversed(range(m)):
            S = np.zeros((self.G[li],) + cur.shape[1:])
            np.add.at(S, idx, cur)
            U = Us[li]
            quad = np.einsum("g...i,ij,g...j->g...", U, Sinvs[li], U)
            q = self.lvs[li].q
            lprior = -0.5 * (q * math.log(2 * math.pi) + ldSs[li] + quad)
            ldC = np.linalg.slogdet(self.Cs[li])[1]
            T = lprior + (self.lw[li] - self.lphi[li]) + ldC.reshape((-1,) + (1,) * (li + 1))
            S = S + T
            R = special.logsumexp(S, axis=-1)
            if keep:
                store[li] = (S, R, U)
            if self.wgrp[li] is not None:
                R = R * self.wgrp[li].reshape((-1,) + (1,) * li)
            cur = R
            if li > 0:
                idx = self.parent[li]
        return cur, store

    def contrib(self, b) -> np.ndarray:
        r = self._eval(b)
        return r[0] if r is not None else np.full(self.G[0], -np.inf)

    def loglik(self, b) -> float:
        if self.ordered and np.any(np.diff(self.aux(b)) <= 0):
            return -np.inf
        r = self._eval(b)
        if r is None:
            return -np.inf
        ll = float(r[0].sum())
        return ll if np.isfinite(ll) else -np.inf

    def adapt(self, b, steps: int = 20) -> None:
        if not self.adaptive:
            return
        for _ in range(steps):
            r = self._eval(b, keep=True)
            if r is None:
                return
            store = r[1]
            J = None
            change = 0.0
            newmus, newCs = [], []
            for li in range(self.m):
                S, R, U = store[li]
                cond = np.exp(S - R[..., None])
                J = cond if li == 0 else J[self.parent[li]][..., None] * cond      # joint weights
                G_ = J.shape[0]
                Jf, Uf = J.reshape(G_, -1), U.reshape(G_, -1, U.shape[-1])
                mu = np.einsum("gk,gkq->gq", Jf, Uf)
                dev = U - mu.reshape((-1,) + (1,) * (li + 1) + (mu.shape[1],))
                devf = dev.reshape(Uf.shape)
                cov = np.einsum("gk,gki,gkj->gij", Jf, devf, devf)
                C = self.Cs[li].copy()
                if li > 0:
                    Up = store[li - 1][2][self.parent[li]]
                    Jp = J.sum(axis=-1).reshape(G_, -1)
                    Upf = Up.reshape(G_, -1, Up.shape[-1])
                    mup = np.einsum("gk,gkq->gq", Jp, Upf)
                    dp = Up - mup.reshape((-1,) + (1,) * li + (mup.shape[1],))
                    dpf = dp.reshape(Upf.shape)
                    Vp = np.einsum("gk,gki,gkj->gij", Jp, dpf, dpf)
                    dpx = np.broadcast_to(dp[..., None, :], dev.shape[:-1] + (dp.shape[-1],)).reshape(G_, -1, dp.shape[-1])
                    X_ = np.einsum("gk,gki,gkj->gij", Jf, devf, dpx)
                    try:
                        Bn = np.einsum("gij,gjk->gik", X_, np.linalg.pinv(Vp, hermitian=True))
                    except np.linalg.LinAlgError:
                        Bn = self.Bs[li]
                    if np.all(np.isfinite(Bn)):
                        cov = cov - np.einsum("gij,gkj->gik", Bn, X_)
                        change = max(change, float(np.max(np.abs(Bn - self.Bs[li]))))
                        self.Bs[li] = Bn
                fin = np.all(np.isfinite(cov), axis=(1, 2))
                ok = fin & np.all(np.linalg.eigvalsh(np.where(fin[:, None, None], cov, 0)) > 1e-12, axis=1)
                C[ok] = _chol_rows(cov[ok], C[ok])
                if not ok.all():        # see GLMM.adapt
                    S_ = self.lvs[li].sigma(self._pieces(self.split(b)[1])[li])
                    C[~ok] = np.linalg.cholesky(S_ + 1e-300 * np.eye(S_.shape[0]))
                change = max(change, float(np.max(np.abs(mu - self.mus[li]))), float(np.max(np.abs(C - self.Cs[li]))))
                newmus.append(mu)
                newCs.append(C)
            self.mus, self.Cs = newmus, newCs
            if change < 1e-8:
                break

    def scores(self, b) -> np.ndarray:
        k = b.size
        Sc = np.zeros((self.G[0], k))
        for i in range(k):
            h = 1e-5 * (abs(b[i]) + 1)
            e = np.zeros(k)
            e[i] = h
            Sc[:, i] = (self.contrib(b + e) - self.contrib(b - e)) / (2 * h)
        return Sc

    def _ll_moving(self, x, state) -> float:
        """The log likelihood with the quadrature points adapted to x, from the state given."""
        self.mus = [a.copy() for a in state[0]]
        self.Cs = [a.copy() for a in state[1]]
        self.Bs = [None if a is None else a.copy() for a in state[2]]
        self.adapt(x)
        return self.loglik(x)

    def fun(self, b, todo):
        # the quadrature points follow b also in the line search: points adapted to a larger
        # variance overstate the likelihood of a smaller one (compat 1916)
        self.adapt(b)
        ll = self.loglik(b)
        if todo == 0 or not np.isfinite(ll):
            return (ll if np.isfinite(ll) else -np.inf), None, None
        if not self.adaptive:
            f = self.loglik
        else:
            # derivatives of the likelihood whose points follow the parameters: the standard
            # errors of Stata (compat 1916)
            state = ([a.copy() for a in self.mus], [a.copy() for a in self.Cs],
                     [None if a is None else a.copy() for a in self.Bs])
            f = lambda x: self._ll_moving(x, state)        # noqa: E731
        k = b.size
        hs = 1e-5 * (np.abs(b) + 1)
        g = np.zeros(k)
        for i in range(k):
            e = np.zeros(k)
            e[i] = hs[i]
            g[i] = (f(b + e) - f(b - e)) / (2 * hs[i])
        hh = 1e-4 * (np.abs(b) + 1)
        H = np.zeros((k, k))
        for i in range(k):
            for j in range(i, k):
                ei = np.zeros(k)
                ej = np.zeros(k)
                ei[i] += hh[i]
                ej[j] += hh[j]
                H[i, j] = H[j, i] = (f(b + ei + ej) - f(b + ei - ej) - f(b - ei + ej) + f(b - ei - ej)) \
                    / (4 * hh[i] * hh[j])
        if self.adaptive:
            self._ll_moving(b, state)
        return ll, g, H


# ---------------------------------------------------------------------------
# the commands
# ---------------------------------------------------------------------------

# command -> (family, link, title, model, reference model of the LR test)
_CMDS = {
    "melogit": ("bernoulli", "logit", "Mixed-effects logistic regression", "logistic", "logistic model"),
    "meprobit": ("bernoulli", "probit", "Mixed-effects probit regression", "probit", "probit model"),
    "mecloglog": ("bernoulli", "cloglog", "Mixed-effects cloglog regression", "cloglog", "cloglog model"),
    "mepoisson": ("poisson", "log", "Mixed-effects Poisson regression", "Poisson", "Poisson model"),
    "menbreg": ("nbinomial", "log", "Mixed-effects nbinomial regression", "nbinomial", "nbinomial model"),
    "meologit": ("ordinal", "logit", "Mixed-effects ologit regression", "ologit", "ologit model"),
    "meoprobit": ("ordinal", "probit", "Mixed-effects oprobit regression", "oprobit", "oprobit model"),
    "meglm": (None, None, "Mixed-effects GLM", "", ""),
    "meqrlogit": ("bernoulli", "logit", "Mixed-effects logistic regression", "logistic", "logistic model"),
    "meqrpoisson": ("poisson", "log", "Mixed-effects Poisson regression", "Poisson", "Poisson model"),
}
_FAMILIES = ("gaussian", "bernoulli", "binomial", "gamma", "nbinomial", "ordinal", "poisson")
_DEFAULT_LINK = {"gaussian": "identity", "bernoulli": "logit", "binomial": "logit", "gamma": "log",
                 "nbinomial": "log", "ordinal": "logit", "poisson": "log"}
_LINKS_OK = {"gaussian": ("identity", "log"), "bernoulli": ("logit", "probit", "cloglog"),
             "binomial": ("logit", "probit", "cloglog"), "gamma": ("log",), "nbinomial": ("log",),
             "ordinal": ("logit", "probit"), "poisson": ("log",)}
# reference model of the LR test of meglm (VERIFY the names)
_VS_GLM = {("gaussian", "identity"): "linear model", ("gaussian", "log"): "linear model", ("bernoulli", "logit"): "logistic model",
           ("bernoulli", "probit"): "probit model", ("bernoulli", "cloglog"): "cloglog model",
           ("binomial", "logit"): "logistic model", ("binomial", "probit"): "probit model",
           ("binomial", "cloglog"): "cloglog model", ("gamma", "log"): "gamma model",
           ("nbinomial", "log"): "nbinomial model", ("ordinal", "logit"): "ologit model",
           ("ordinal", "probit"): "oprobit model", ("poisson", "log"): "Poisson model"}
_FAM_TITLE = {"gaussian": "Gaussian", "bernoulli": "Bernoulli", "binomial": "binomial",
              "gamma": "gamma", "nbinomial": "negative binomial", "ordinal": "ordinal",
              "poisson": "Poisson"}                                  # compat 1902
_LINK_TITLE = {"identity": "identity", "logit": "logit", "probit": "probit", "log": "log",
               "cloglog": "complementary log-log"}                   # compat 1902
_GLOBAL = {"or": 2, "irr": 3, "intpoints": 4, "intmethod": 7, "nolog": 4, "log": 3, "iterate": 4,
           "tolerance": 3, "ltolerance": 4, "nrtolerance": 4, "level": 1, "noheader": 6, "nogroup": 5,
           "nolrtest": 5, "noconstant": 4, "offset": 3, "exposure": 3, "binomial": 3, "covariance": 3,
           "startgrid": 6, "nofetable": 5, "noretable": 5, "estmetric": 4, "family": 3, "link": 4,
           "dispersion": 4, "asis": 4, "laplace": 3, "stddeviations": 3, "variance": 3, "refineopts": 6}


def _word(text: str, choices, what: str) -> str:
    t = text.strip().lower()
    hit = [c for c in choices if c.startswith(t)] if t else []
    if len(hit) != 1 and t in choices:
        hit = [t]
    if len(hit) != 1:
        raise StataError(198, f"{what}({text.strip()}) invalid")   # VERIFY
    return hit[0]


def _resolve(cmdname: str, o: dict) -> dict:
    """Family, link and titles of the command, from its name and the family()/link() options of
    meglm."""
    fam, link, title, model, vs = _CMDS[cmdname]
    disp = ""
    nb = ""
    if cmdname == "meglm":
        words = str(o.get("family") or "gaussian").split()
        fam = _word(words[0], _FAMILIES, "family")
        rest = words[1:]
        if fam == "nbinomial":
            disp = _word(rest[0], ("mean", "constant"), "family") if rest else "mean"
        elif fam == "binomial":
            nb = rest[0] if rest else "1"
        elif rest:
            raise StataError(198, f"family({' '.join(words)}) invalid")   # VERIFY
        link = _word(str(o["link"]), ("identity", "log", "logit", "probit", "cloglog"), "link") \
            if o.get("link") else _DEFAULT_LINK[fam]
        if link not in _LINKS_OK[fam]:
            raise StataError(198, f"link({link}) with family({fam}) is not available yet in OpenDTA")   # OpenDTA
        vs = _VS_GLM[(fam, link)]
        model = vs[:-6] if vs.endswith(" model") else ""
    elif cmdname == "menbreg":
        disp = _word(str(o["dispersion"]), ("mean", "constant"), "dispersion") if o.get("dispersion") else "mean"
    elif o.get("dispersion"):
        raise StataError(198, "option dispersion() not allowed")
    if o.get("binomial") is not None and o.get("binomial") is not False:
        if cmdname not in ("melogit", "meprobit", "mecloglog"):
            raise StataError(198, "option binomial() not allowed")
        fam, nb = "binomial", str(o["binomial"]).strip()
    if fam == "binomial" and not nb:
        nb = "1"
    return {"cmd": cmdname, "fam": fam, "link": link, "title": title, "model": model, "vs": vs,
            "disp": disp, "nbin": nb}


def _mean_fn(sp: dict):
    """Inverse link of the family: mean of the outcome given the linear predictor."""
    link = sp["link"]
    if link == "logit":
        return special.expit
    if link == "probit":
        return st.norm.cdf
    if link == "cloglog":
        return lambda e: -np.expm1(-np.exp(np.clip(e, -700, 6.5)))
    if link == "log":
        return lambda e: np.exp(np.clip(e, -700, 700))
    return lambda e: e


def _fixed_start(sp, y, off, p, constant, n_nb):
    """Starting values (beta, auxiliary) of the fixed-effects model."""
    fam, link = sp["fam"], sp["link"]
    b0 = np.zeros(p)
    ybar = float(y.mean())
    if fam == "gaussian":
        c = ybar - float(off.mean())
        a0 = [math.log(max(float(y.var()), 1e-8))]
    elif fam == "gamma":
        c = math.log(max(ybar, 1e-10)) - float(off.mean())
        a0 = [0.5 * math.log(max(float(y.var()) / max(ybar * ybar, 1e-10), 1e-2))]
    elif fam == "nbinomial":
        c = math.log(max(ybar, 1e-10)) - float(off.mean())
        over = max(float(y.var()) / max(ybar, 1e-10) - 1.0, 0.05)
        a0 = [math.log(over / max(ybar, 1e-10)) if sp["disp"] == "mean" else math.log(over)]
    elif fam == "poisson":
        c = math.log(max(ybar, 1e-10)) - float(off.mean())
        a0 = []
    else:           # bernoulli, binomial
        pr = min(max(float((y / n_nb).mean()), 1e-4), 1 - 1e-4)
        c = {"logit": math.log(pr / (1 - pr)), "probit": float(st.norm.ppf(pr)),
             "cloglog": math.log(-math.log(1 - pr))}[link]
        a0 = []
    if constant:
        b0[-1] = c
    return b0, np.array(a0, dtype=float)


def _num_gh(gfun, b):
    """Hessian by central differences of the gradient."""
    k = b.size
    H = np.zeros((k, k))
    for i in range(k):
        hi = 1e-5 * (abs(b[i]) + 1)
        e = np.zeros(k)
        e[i] = hi
        H[:, i] = (gfun(b + e) - gfun(b - e)) / (2 * hi)
    return (H + H.T) / 2


@dataclass
class Fit:
    """What a command hands to _fit_me once it has read its data: design, offset, level, the
    conditional distribution of the outcome with its auxiliary parameters, and the details that
    differ from one command to the next."""
    smp: object
    o: dict
    sp: dict
    X: np.ndarray
    off: np.ndarray
    lv: Level
    fam_f: object
    naux: int
    b0: np.ndarray                    # starting values of the fixed-effects model (with auxiliary)
    npts: int
    method: str
    family: str
    link: str
    aux_eq: list
    off_name: str = ""
    ordered: bool = False
    analytic: bool = False            # analytic Hessian of the fixed-effects model
    aux_in_table: bool = False
    aux_rows: object = None           # (aux, Va, z) -> rows of the coefficient table
    resid: object = None              # latent residual variance of the outcome (float or f(aux))
    nbin: object = None
    nlev: int = 0
    qr: bool = False
    eform: str = ""
    predict: object = None
    predict_multi: dict = field(default_factory=dict)
    macro_edit: object = None
    cmdline: str = ""                 # e(cmdline) when Stata rebuilds it (meglm, mestreg)
    late_error: str = ""              # an error that Stata only raises after the fit
    predict_opts: dict | None = None
    cats: list = field(default_factory=list)      # outcome values of the ordinal models
    lvs: list = field(default_factory=list)       # all the levels, the highest first (default [lv])
    wobs: object = None               # observation weights ([fweight], [pweight], [iweight])
    wgrp: object = None               # group weights of each level (fweight(), pweight()) or None
    level_opts: list = field(default_factory=list)   # the options of each random-effects equation


def cmd_me(s: "Session", args: str, cmdname: str) -> None:
    from .postest import replay
    from .regress import Sample, cmdline
    if not args.strip() or args.strip().startswith(","):
        replay(s, cmdname, args)
        return
    fixed, levels_spec, glob = split_levels(args)
    import re
    m_ = re.search(r"\b(exposure|offset)\s*\(\s*\w+\s*\)", glob)
    if m_ and not cmdname.startswith("meqr"):      # compat 1902: these belong to the fixed-effects equation (before the ||)
        raise StataError(198, f"option {re.sub(chr(32), '', m_.group(0))} allowed only in the fixed-effects equation")
    if not levels_spec:
        raise StataError(198, "random-effects equation required")   # VERIFY
    qr_ = cmdname.startswith("meqr")
    c = find_top(fixed, ",")
    fixed_body, fe_opts = (fixed, "") if c == -1 else (fixed[:c], fixed[c + 1:])
    opts_text = " ".join(x for x in (fe_opts.strip(), glob.strip()) if x)
    extra = " ".join([lv for lv, _r, _o in levels_spec] + [r for _l, r, _o in levels_spec if r])
    for _l, _r, _o in levels_spec:
        for _wk in ("fweight", "pweight"):
            if isinstance(_o.get(_wk), str):
                extra += " " + _o[_wk].strip()
    try:
        smp = Sample(s, fixed_body + (f", {opts_text}" if opts_text else ""), opts_spec=_GLOBAL,
                     weights=() if qr_ else ("fweight", "iweight", "pweight"), extra_vars=extra)
    except StataError as e_:
        if e_.message == "aweights not allowed":
            raise StataError(e_.rc, "aweight not allowed") from None     # compat 1916
        raise
    o = smp.o
    sp = _resolve(cmdname, o)
    fam, link = sp["fam"], sp["link"]
    qr = cmdname.startswith("meqr")          # the older commands: same model, layout of mixed
    for _opt in ("stddeviations", "variance"):
        if o.get(_opt) and not qr:      # compat 1902: only the meqr commands accept them
            raise StataError(198, f"option {_opt} not allowed")
    if fam == "ordinal":
        smp.constant = False
    out = s.output
    ds = s.data
    y = smp.y.astype(np.float64)
    n = smp.n
    nbin = np.ones(n)
    if fam == "bernoulli":
        y = (y != 0).astype(float)
    elif fam == "binomial":
        t = sp["nbin"]
        try:
            nbin = np.full(n, float(t))
        except ValueError:
            nbin = ds.get(t).data[smp.mask].astype(float)
        if np.any(y < 0) or np.any(y > nbin):
            raise StataError(459, f"{smp.depname} must be between 0 and the number of trials")   # VERIFY
    elif fam == "gamma":
        if np.any(y <= 0):
            raise StataError(459, f"{smp.depname} must be greater than zero")   # VERIFY
    elif fam in ("poisson", "nbinomial"):
        if np.any(y < 0):
            raise StataError(459, f"{smp.depname} must be greater than or equal to zero")   # VERIFY
    yi = None
    if fam == "ordinal":
        cats = sorted(set(y.tolist()))
        if len(cats) < 2:
            raise StataError(2000, "outcome does not vary")   # VERIFY
        if len(cats) > 50:
            raise StataError(149, f"{smp.depname}:  too many outcomes")   # VERIFY
        yi = np.searchsorted(cats, y)
    fvars.drop_collinear(smp.cols, smp.constant)
    est_cols = [col for col in smp.cols if col.values is not None]
    X = np.column_stack([col.values for col in est_cols] + ([np.ones(n)] if smp.constant else [])) \
        if (est_cols or smp.constant) else np.zeros((n, 0))
    off = np.zeros(n)
    off_name = ""
    if o.get("offset") or o.get("exposure"):
        if o.get("exposure") and cmdname in ("melogit", "meprobit", "mecloglog", "meologit", "meoprobit"):
            raise StataError(198, "option exposure() not allowed")
        nm = str(o.get("offset") or o.get("exposure")).strip()
        v = ds.get(nm).data[smp.mask].astype(float)
        off = np.log(v) if o.get("exposure") else v
        off_name = f"ln({nm})" if o.get("exposure") else nm
    lvspec = levels_spec[-1]
    if lvspec[2].get("covariance") is None and o.get("covariance"):
        lvspec[2]["covariance"] = o["covariance"]
    if any(lv_ == "_all" for lv_, _r, _o in levels_spec):
        raise StataError(198, "crossed random effects (_all:) are not available for this command")   # OpenDTA
    levels = build_levels(s, smp.mask, levels_spec, n)
    for lv_, (_l, _r, ro_) in zip(levels, levels_spec):
        if lv_.cov == "identity" and lv_.q > 1 and not ro_.get("covariance"):
            lv_.cov = "independent"
    lv = levels[0]
    wobs, wgrp = level_weights(s, smp, levels, levels_spec)
    npts = int(o["intpoints"]) if o.get("intpoints") else 7
    method = str(o.get("intmethod") or "mvaghermite").strip().lower()
    method = {"mv": "mvaghermite", "mvaghermite": "mvaghermite", "mc": "mcaghermite",
              "mcaghermite": "mcaghermite", "gh": "ghermite", "ghermite": "ghermite",
              "laplace": "laplace"}.get(method, method)
    if qr and o.get("laplace"):
        method = "laplace"
    # the conditional distribution of the outcome
    if fam == "bernoulli" and link in ("logit", "probit"):
        fam_f, naux = _family(link, y), 0
    elif fam == "poisson":
        fam_f, naux = _family("poisson", y), 0
    elif fam == "ordinal":
        fam_f, naux = ordinal_family(link, yi)
    else:
        fam_f, naux = family_of("binomial" if fam == "bernoulli" else fam, link, y, nbin=nbin,
                                disp=sp["disp"])
    p = X.shape[1]
    analytic = naux == 0 and fam != "binomial" and not (fam == "bernoulli" and link == "cloglog")
    if analytic:
        b0 = _fixed_start(sp, y, off, p, smp.constant, nbin)[0]
    elif fam == "ordinal":
        cum = np.array([np.mean(yi <= j) for j in range(len(cats) - 1)])
        b0 = np.concatenate([np.zeros(p), (special.logit if link == "logit" else st.norm.ppf)(cum)])
    else:
        b0b, a0 = _fixed_start(sp, y, off, p, smp.constant, nbin)
        if fam == "gaussian" and link == "log":
            b0b = np.zeros(p)
            if smp.constant:
                b0b[p - 1] = math.log(max(float(y.mean()), 0.1))
            a0 = np.array([math.log(max(float(np.var(y)), 1e-8))])
        elif fam == "gaussian" and p:
            b0b = np.linalg.lstsq(X, y - off, rcond=None)[0]
            a0 = np.array([math.log(max(float(np.mean((y - off - X @ b0b) ** 2)), 1e-8))])
        b0 = np.concatenate([b0b, a0])
    if fam == "ordinal":
        aux_eq = [f"/cut{j + 1}" for j in range(naux)]
    elif fam == "nbinomial":
        aux_eq = ["/lnalpha" if sp["disp"] == "mean" else "/lndelta"]
    elif fam == "gamma":
        aux_eq = ["/logs"]                                      # compat 1902
    elif fam == "gaussian":
        aux_eq = [f"/var(e.{smp.depname})"]
    else:
        aux_eq = []
    resid = {("bernoulli", "logit"): math.pi ** 2 / 3, ("bernoulli", "probit"): 1.0,
             ("bernoulli", "cloglog"): math.pi ** 2 / 6, ("binomial", "logit"): math.pi ** 2 / 3,
             ("binomial", "probit"): 1.0, ("binomial", "cloglog"): math.pi ** 2 / 6,
             ("ordinal", "logit"): math.pi ** 2 / 3, ("ordinal", "probit"): 1.0}.get((fam, link))
    if fam == "gaussian":
        resid = lambda a: float(math.exp(a[0]))                 # noqa: E731
    eform = "Odds Ratio" if (link == "logit" and o.get("or")) else \
        ("IRR" if (fam in ("poisson", "nbinomial") and o.get("irr")) else "")
    F = Fit(smp=smp, o=o, sp=sp, X=X, off=off, off_name=off_name, lv=lv, lvs=levels, wobs=wobs, wgrp=wgrp,
            level_opts=[ro_ for _l, _r, ro_ in levels_spec],
            fam_f=fam_f, naux=naux,
            ordered=fam == "ordinal", b0=b0, analytic=analytic, aux_eq=aux_eq,
            aux_in_table=fam in ("ordinal", "nbinomial", "gamma"), resid=resid, nbin=nbin,
            nlev=(len(cats) if fam == "ordinal" else 0), qr=qr, eform=eform, npts=npts, method=method,
            family=fam, link=link)
    if fam == "ordinal":
        F.predict_multi = {"": _pr_multi, "pr": _pr_multi}
        F.cats = [float(c) for c in cats]
    if cmdname == "meglm":      # Stata rebuilds the command line (compat 1902; VERIFY for other forms)
        F.cmdline = (f"meglm {fixed_body.strip()}   ,  "
                     + "".join(f"|| {lvn}: {rev_.strip()}{' ' if rev_.strip() else ''}" for lvn, rev_, _ro in levels_spec)
                     + ", " + "    ".join(_split_opts(glob)))
    _fit_me(s, args, cmdname, F)


def level_weights(s: "Session", smp, levels, levels_spec):
    """Observation weights (from [fweight], [pweight] or [iweight]) and the group weights of
    each level (its fweight() or pweight() option, which must be constant within a group)."""
    wobs = None
    if smp.wtype:
        wobs = np.asarray(smp.w, dtype=float)[smp.mask]
    wgrp = []
    kinds = {smp.wtype} if smp.wtype else set()
    for lv, (_l, _r, ro) in zip(levels, levels_spec):
        nm = next(((k, ro[k]) for k in ("fweight", "pweight") if isinstance(ro.get(k), str)), None)
        if nm is None:
            wgrp.append(None)
            continue
        kinds.add(nm[0])
        from ..core.varlist import resolve_name
        v = s.data.get(resolve_name(s.data, nm[1].strip())).data[smp.mask].astype(float)
        G = int(lv.gid.max()) + 1
        first = np.full(G, np.nan)
        first[lv.gid[::-1]] = v[::-1]
        if np.any(v != first[lv.gid]):
            raise StataError(198, f"weights in variable {nm[1].strip()} not constant within groups defined by: {lv.var}")   # compat 1916
        if np.any(first < 0):
            raise StataError(402, "negative weights encountered")   # VERIFY
        wgrp.append(first)
    if len(kinds - {"iweight"}) > 1:
        raise StataError(198, "frequency weights and sampling weights may not be combined")   # VERIFY
    if all(w is None for w in wgrp):
        wgrp = None
    return wobs, wgrp


def _split_opts(text: str) -> list[str]:
    """The options of a command line as separate words (parentheses kept together)."""
    out, depth, cur = [], 0, ""
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch.isspace() and depth == 0:
            if cur:
                out.append(cur)
            cur = ""
        else:
            cur += ch
    if cur:
        out.append(cur)
    return out


def _me_vce(F: "Fit"):
    """vce(oim), vce(robust) or vce(cluster clustvar); sampling weights imply vce(robust)."""
    smp = F.smp
    vce = getattr(smp, "vce", "ols")
    clustvar = getattr(smp, "clustvar", "") or ""
    if vce not in ("ols", "oim", "robust", "cluster"):
        raise StataError(198, f"vcetype '{vce}' not allowed")   # VERIFY
    pw = getattr(smp, "wtype", "") == "pweight" or any(
        isinstance(ro.get("pweight"), str) for ro in (F.level_opts or []))
    if vce in ("ols", "oim"):
        vce = "robust" if pw else "oim"
    if vce != "oim" and F.qr:
        raise StataError(198, f"option vce({vce}) not allowed")   # VERIFY
    if clustvar:
        gtop = (F.lvs or [F.lv])[0].gid
        cl = np.asarray(smp.cluster)
        first = np.zeros(int(gtop.max()) + 1, dtype=cl.dtype)
        first[gtop[::-1]] = cl[::-1]
        if np.any(first[gtop] != cl):
            raise StataError(198, f"highest-level groups are not nested within {clustvar}")   # VERIFY
    return vce, clustvar


def _robust_v(model, b, V, F: "Fit", clustvar: str):
    """Sandwich variance with the scores of the highest-level groups summed within the clusters
    (each highest-level group is a cluster under vce(robust)); multiplier G/(G-1) (VERIFY)."""
    Sc = model.scores(b)
    gtop = model.gtop
    if clustvar:
        cl = np.asarray(F.smp.cluster)
        G = Sc.shape[0]
        first = np.zeros(G, dtype=cl.dtype)
        first[gtop[::-1]] = cl[::-1]
        _, cid = np.unique(first, return_inverse=True)
        Sc = np.array([Sc[cid == c].sum(axis=0) for c in range(cid.max() + 1)])
    g = Sc.shape[0]
    M_ = Sc.T @ Sc * (g / (g - 1) if g > 1 else 1.0)
    Vr = V @ M_ @ V
    return (Vr + Vr.T) / 2, g


def _obs_weight(F: "Fit"):
    """The weight of each observation in the model without random effects: its own weight times
    the weights of its groups."""
    w = None if F.wobs is None else np.asarray(F.wobs, dtype=float).copy()
    for lv, wg in zip(F.lvs or [F.lv], F.wgrp or []):
        if wg is not None:
            w = (np.ones(F.X.shape[0]) if w is None else w) * wg[lv.gid]
    return w


def _fit_fixed(F: "Fit", kw: dict):
    """The model without random effects: its log likelihood is the reference of the LR test and
    its estimates start the full model."""
    X, off, naux, b0 = F.X, F.off, F.naux, F.b0
    p = X.shape[1]
    wt = _obs_weight(F)
    if wt is None:
        fam_f = F.fam_f
    else:
        def fam_f(eta, a=None):
            r = F.fam_f(eta, a) if a is not None else F.fam_f(eta)
            return tuple(x * wt[:, None] for x in r)
    if F.analytic:
        def fe_fun(b, todo):
            eta = (X @ b + off)[:, None]
            lf, s1, s2 = fam_f(eta)
            ll = float(lf.sum())
            if todo == 0:
                return ll, None, None
            return ll, X.T @ s1[:, 0], (X * s2[:, 0][:, None]).T @ X
    else:
        def fe_grad(b):
            eta = (X @ b[:p] + off)[:, None]
            lf, s1, s2 = fam_f(eta, b[p:])
            g = np.zeros(b.size)
            g[:p] = X.T @ s1[:, 0]
            for m in range(naux):
                h = 1e-5 * (abs(b[p + m]) + 1)
                e = np.zeros(naux)
                e[m] = h
                g[p + m] = float(np.sum(fam_f(eta, b[p:] + e)[0] - fam_f(eta, b[p:] - e)[0])) / (2 * h)
            return g

        def fe_fun(b, todo):
            if F.ordered and np.any(np.diff(b[p:]) <= 0):
                return -np.inf, None, None
            eta = (X @ b[:p] + off)[:, None]
            ll = float(fam_f(eta, b[p:])[0].sum())
            if todo == 0:
                return ll, None, None
            return ll, fe_grad(b), _num_gh(fe_grad, b)
    return maximize(fe_fun, b0, **kw)


def _fit_me(s: "Session", args: str, cmdname: str, F: "Fit") -> None:
    """Fits the mixed-effects model with the conditional distribution F.fam_f and posts and
    shows the results (the part that meglm and mestreg have in common)."""
    from .regress import cmdline
    smp, o, sp, X, off, lv = F.smp, F.o, F.sp, F.X, F.off, F.lv
    fam, link, qr, naux = F.family, F.link, F.qr, F.naux
    off_name, aux_eq, npts, method = F.off_name, F.aux_eq, F.npts, F.method
    out = s.output
    ds = s.data
    n = smp.n
    lvs = F.lvs or [lv]
    if len(lvs) == 1:
        model = GLMM(fam, np.zeros(n), X, off, lv, npts, method, fam=F.fam_f, naux=naux, ordered=F.ordered,
                     wobs=F.wobs, wgrp=(F.wgrp[0] if F.wgrp else None))
    else:
        if method == "laplace":
            raise StataError(198, "intmethod(laplace) with more than one level is not available yet in OpenDTA")   # OpenDTA
        model = MGLMM(X, off, lvs, npts, method, F.fam_f, naux=naux, ordered=F.ordered, wobs=F.wobs,
                      wgrp=F.wgrp)
    vce, clustvar = _me_vce(F)
    robust = vce in ("robust", "cluster")
    kw = ml_options(o)
    p = X.shape[1]
    est_cols = [col for col in smp.cols if col.values is not None]
    rfe = _fit_fixed(F, kw)
    # starting values for the variance components: best of a small grid (VERIFY)
    best = None
    for ls in (math.log(0.5), 0.0, math.log(0.25), math.log(1.5)):
        ths = []
        for lv_ in lvs:
            th_ = np.zeros(lv_.npar())
            th_[:min(lv_.q, th_.size) if lv_.cov != "exchangeable" else 1] = ls
            ths.append(th_)
        model.reset()
        bb = np.concatenate([rfe.b] + ths)
        model.adapt(bb)
        ll = model.loglik(bb)
        if best is None or ll > best[0]:
            best = (ll, bb, model.mu.copy(), model.C.copy())
    model.mu, model.C = best[2], best[3]
    lab = "log pseudolikelihood" if robust else "log likelihood"
    if not o.get("nolog"):
        if qr:
            # no fixed-effects fit is shown; the first stage refines the starting values (VERIFY)
            import re as _re
            if not _re.search(r"iterate\(\s*0\s*\)", str(o.get("refineopts") or "")):   # compat 1902
                out.write("\nRefining starting values:\n\n", "text")
                out.write(f"{'Iteration 0:':<14} {lab} = {_g(best[0], '%10.0g'):>10}\n", "text")
            out.write("\nPerforming gradient-based optimization:\n\n", "text")
        else:
            # the fixed-effects fit and the grid show "log likelihood" even with vce(robust) or
            # sampling weights (compat 1916)
            out.write("\nFitting fixed-effects model:\n\n", "text")
            print_log(out, rfe)
            out.write("\nRefining starting values:\n\n", "text")
            out.write(f"{'Grid node 0:':<14} log likelihood = {_g(best[0], '%10.0g'):>10}\n", "text")
            out.write("\nFitting full model:\n\n", "text")
    r = maximize(model.fun, best[1], **kw)
    if not o.get("nolog"):
        print_log(out, r, lab)
    b = r.b
    try:
        V = np.linalg.inv(-r.H)
        V = (V + V.T) / 2
    except np.linalg.LinAlgError:
        V = np.full((b.size, b.size), np.nan)
    n_clust = 0
    V_model = V
    if robust:
        V, n_clust = _robust_v(model, b, V, F, clustvar)
        # the note and e(N_clust) count the highest-level groups also with vce(cluster) (compat
        # 1916: "30 clusters in cl" for 6 clusters of 5 groups; VERIFY)
        n_clust = int(lvs[0].gid.max()) + 1
    beta, th = model.split(b)
    aux = model.aux(b)
    Vb = V[:p, :p]
    Va = V[p:p + naux, p:p + naux]
    Vth = V[p + naux:, p + naux:]
    slopes = list(range(len(est_cols)))
    try:
        wchi2 = float(beta[slopes] @ np.linalg.solve(Vb[np.ix_(slopes, slopes)], beta[slopes])) if slopes else SYS
    except np.linalg.LinAlgError:
        wchi2 = SYS
    wp = float(chi2_tail(wchi2, len(slopes))) if slopes and wchi2 < SYS else SYS
    chi2_c = max(2 * (r.ll - rfe.ll), 0.0)
    names = [col.name if not col.omitted else fvars.omitted_name(col.name) for col in smp.cols]
    if smp.constant:
        names.append("_cons")
    kk = len(names)
    pos = [i for i, col in enumerate(smp.cols) if col.values is not None] + ([kk - 1] if smp.constant else [])
    bf = np.zeros(kk)
    Vf = np.zeros((kk, kk))
    bf[pos] = beta
    Vf[np.ix_(pos, pos)] = Vb
    bmap = {col.name: bf[i] for i, col in enumerate(smp.cols)}
    semap = {col.name: float(np.sqrt(max(Vf[i, i], 0))) for i, col in enumerate(smp.cols)}
    if smp.constant:
        bmap["_cons"], semap["_cons"] = bf[-1], float(np.sqrt(max(Vf[-1, -1], 0)))
    if len(lvs) == 1:
        th_names = [f"{t}:_cons" for t in lv.labels()]
    else:       # VERIFY the names of e(b) with several levels
        th_names = [f"{t.replace('_', str(li + 1) + '_', 1)}:_cons" for li, lv_ in enumerate(lvs) for t in lv_.labels()]
    # auxiliary parameters: names in e(b) and the rows of the table are in F.aux_eq
    ball = np.concatenate([bf, aux, th])
    Vall = np.zeros((ball.size, ball.size))
    Vall[:kk, :kk] = Vf
    Vall[kk:kk + naux, kk:kk + naux] = Va
    Vall[kk:kk + naux, kk + naux:] = V[p:p + naux, p + naux:]
    Vall[kk + naux:, kk:kk + naux] = V[p + naux:, p:p + naux]
    Vall[kk + naux:, kk + naux:] = Vth
    bidx = pos + list(range(kk, ball.size))
    V_model_all = np.zeros_like(Vall)
    V_model_all[np.ix_(bidx, bidx)] = V_model
    est = Estimates(cmdname, smp.depname, names + ["_cons"] * (naux + th.size), ball, Vall, n, "z", None)
    est.eqnames = [smp.depname] * kk + aux_eq + [f"/{t.split(':')[0]}" for t in th_names]
    est.rows = rows_from_columns(smp.cols, bmap, semap, ds, constant=smp.constant)
    if off_name:
        est.rows.append(CoefRow("fixed", off_name, eq="(exposure)" if o.get("exposure") else "(offset)"))
    if F.aux_in_table:
        est.rows.append(CoefRow("sep"))
        if F.aux_rows is not None:
            est.rows += F.aux_rows(aux, Va, st.norm.ppf(1 - (1 - smp.level / 100) / 2))
        else:
            for j, nm_ in enumerate(aux_eq):
                se_a = float(np.sqrt(max(Va[j, j], 0)))
                est.rows.append(CoefRow("coef", nm_, float(aux[j]), se_a, name="_cons", noeform=True))
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = "Robust" if robust else ""
    est.robust = robust
    if robust:
        est.n_clust = n_clust
        est.clustvar = clustvar or lvs[0].var
    est.offset_name = off_name
    resid = F.resid(aux) if callable(F.resid) else F.resid
    est.glmm = {"model": model, "theta": th, "Vth": Vth, "beta": beta, "levels": lvs, "aux": aux,
                "Va": Va, "aux_eq": aux_eq, "resid_var": resid, "mask": smp.mask.copy(), "nbin": F.nbin,
                "nlev": F.nlev}
    # frequency weights: the observations of a group are counted with their [fweight] and the
    # groups (and the observations) with the fweight() of their level and of the levels above
    # (compat 1916)
    fwo = np.asarray(F.wobs, dtype=float) if (F.wobs is not None and getattr(smp, "wtype", "") == "fweight") \
        else np.ones(n)
    fwg = np.ones(n)
    groups = []
    for li_, lv_ in enumerate(lvs):
        wg_ = F.wgrp[li_] if (F.wgrp and isinstance(F.level_opts[li_].get("fweight"), str)) else None
        if wg_ is not None:
            fwg = fwg * wg_[lv_.gid]
        cnt = np.bincount(lv_.gid, weights=fwo)
        ng = np.bincount(lv_.gid, weights=fwg) / np.bincount(lv_.gid)    # the group's weight
        groups.append((lv_.var, int(round(ng.sum())), int(round(cnt.min())), float(np.sum(cnt * ng) / ng.sum()),
                       int(round(cnt.max()))))
    nvc = sum(lv_.npar() for lv_ in lvs)
    N_ = float(np.sum(fwo * fwg))
    est.extra = {"N": N_, "ll": r.ll, "chi2": wchi2, "p": wp, "df_m": len(slopes), "ll_c": rfe.ll,
                 "chi2_c": chi2_c, "nvc": nvc, "spec": sp, "method": method, "qr": qr,
                 "sd": bool(o.get("stddeviations")), "robust": robust, "n_clust": n_clust,
                 "clustvar": clustvar or lvs[0].var, "groups": groups}
    est.extra["npts"] = 1 if method == "laplace" else npts
    est.display = _display_qr if qr else display_me
    est.predict = F.predict or _predict_me
    est.predict_opts = F.predict_opts if F.predict_opts is not None else \
        {"mu": 2, "eta": 3, "reffects": 4, "fixedonly": 5, "pr": 2, "outcome": 3, "conditional": 4}
    est.predict_multi = {"reffects": _reffects_me, **F.predict_multi}
    est.estat = {"icc": (3, estat_icc), "group": (3, estat_group)}
    if qr:
        est.estat["recovariance"] = (5, estat_recov)       # compat 1902: unrecognized after the others
    est.icc_ok = cmdname in ("melogit", "meprobit", "meqrlogit", "meqrpoisson") or \
        (cmdname == "meglm" and fam in ("bernoulli", "gaussian"))   # VERIFY for meglm
    est.eform_default = F.eform
    # e() in the order Stata lists it (compat 1002 and 1902)
    nvar_p = sum({"identity": 1, "exchangeable": 1}.get(lv_.cov, lv_.q) for lv_ in lvs)
    ncov_p = nvc - nvar_p
    from ..commands.matrix import Matrix
    from .fvars import base_vars
    from .ml import ilog
    from .results import datasignature
    covs = [c.name for c in smp.cols if c.values is not None]
    family = fam
    names_all = list(est.names)
    grad = np.zeros((1, len(names_all)))
    if r.g is not None and len(r.g) == len(names_all):
        grad[0] = r.g
    ivars = " ".join(lv_.var for lv_ in lvs)
    gmats = [(nm, Matrix(np.array([[g_[j] for g_ in groups]], dtype=float), ["r1"], [lv_.var for lv_ in lvs]))
             for nm, j in (("g_max", 4), ("g_avg", 3), ("g_min", 2), ("N_g", 1))]
    cmdline_ = F.cmdline or cmdline(cmdname, args)
    if qr:
        # meqrlogit and meqrpoisson: the e() of the xtmelogit family. The plain fit lists it in the
        # order below and a fit with refineopts() in the reverse order (compat 1902; VERIFY)
        sigvars = [smp.depname] + [v for v in base_vars(smp.terms) if v != smp.depname] + [lv.var] \
            + [nm for nm in lv.names if nm != "_cons"]
        scalars = [("converged", 1.0 if r.converged else 0.0), ("rc", 0.0), ("p", wp), ("chi2", wchi2),
                   ("df_m", float(len(slopes))), ("N", float(n)), ("p_c", _p_chibar(chi2_c, nvc)),
                   ("chi2_c", chi2_c), ("df_c", float(nvc)), ("ll_c", rfe.ll), ("ll", r.ll),
                   ("k_rc", float(ncov_p)), ("k_rs", float(nvar_p)), ("k_r", float(th.size)),
                   ("k_f", float(len(slopes) + int(smp.constant))), ("k", float(ball.size)),
                   ("reparm_rc", 0.0), ("rank", float(ball.size))]
        macros = [("cmdline", cmdline_), ("datasignaturevars", " ".join(sigvars)),
                  ("datasignature", datasignature(s, sigvars, smp.mask)), ("chi2type", "Wald"),
                  ("n_quad", str(int(est.extra["npts"]))), ("redim", " ".join(str(l_.q) for l_ in lvs)),
                  ("ivars", " ".join(l_.var for l_ in lvs)),   # nested levels: VERIFY
                  ("vartypes", " ".join(_VT[l_.cov] for l_ in lvs)),
                  ("revars", " ".join(" ".join(l_.names) for l_ in lvs)), ("title", sp["title"]),
                  ("technique", "nr"), ("ml_method", "d0"), ("opt", "moptimize"), ("model", sp["model"]),
                  ("estat_cmd", f"{cmdname}_estat"), ("predict", f"{cmdname}_p"), ("cmd", cmdname),
                  ("marginsdefault", "predict(xb)"),
                  ("marginsnotok", "stdp Pearson Deviance Anscombe REFfects RESEs"), ("method", "ML"),
                  ("depvar", smp.depname), ("properties", "b V")]
        extra_mats = list(gmats)
        if o.get("refineopts"):
            scalars.reverse()
            macros.reverse()
            extra_mats.reverse()
        if method == "laplace":
            macros.append(("laplace", "laplace"))
    else:
        # the level variables from the lowest level up (compat 1916)
        sigvars = [lv_.var for lv_ in reversed(lvs)] + [v for v in base_vars(smp.terms) if v != smp.depname] + [smp.depname]
        ordinal = fam == "ordinal"
        scalars = [("rank", float(ball.size)), ("N", N_), ("ic", float(r.iterations)),
                   ("converged", 1.0 if r.converged else 0.0), ("rc", 0.0), ("ll", r.ll),
                   ("k", float(ball.size)), ("k_eq", float(1 + naux + nvc))]
        if ordinal:
            scalars.append(("k_cat", float(F.nlev)))
        scalars += [("k_eq_model", 1.0), ("ll_c", rfe.ll), ("k_rs", float(nvar_p)), ("k_rc", float(ncov_p)),
                    ("k_dv", 1.0)]
        if len(lvs) > 1:        # compat 1916: with nested levels this part comes in reverse order
            scalars.reverse()
        scalars += [("k_r", float(sum(lv_.q for lv_ in lvs))),
                    ("k_f", float(len(slopes)))]   # VERIFY: k_f without _cons
        scalars += [("p_c", _p_chibar(chi2_c, nvc)), ("df_c", float(nvc)), ("chi2_c", chi2_c)]
        scalars += [("p", wp), ("chi2", wchi2), ("df_m", float(len(slopes)))]
        if robust:
            scalars += [("N_clust", float(n_clust))]
        common = [("n_quad", str(int(est.extra["npts"]))), ("ivars", ivars), ("chi2type", "Wald"),
                  ("datasignaturevars", " ".join(sigvars)),
                  ("datasignature", datasignature(s, sigvars, smp.mask)), ("intmethod", method),
                  ("depvar", smp.depname)] + \
            ([("clustvar", clustvar or lvs[0].var)] if robust else []) + \
            [("marginsnotok", "LATent1(passthru) LATent se(passthru)"), ("covariates", " ".join(covs))]
        tail = [("opt", "moptimize")] + ([("vcetype", "Robust")] if robust else []) + [("vce", vce)] + \
            [("user", "_gsem_eval_chol()"), ("ml_method", "gf2"),
             ("technique", "nr"), ("which", "max"), ("properties", "b V")]   # VERIFY vcetype/clustvar
        if getattr(smp, "wtype", ""):
            from .regress import wexp_macros
            tail[1:1] = wexp_macros(smp)   # VERIFY the place
        if ordinal:
            tail.insert(0, ("marginsdefault", " ".join(f"predict(pr outcome({c:g}))" for c in F.cats)))
        if sp["cmd"] == "meglm":
            macros = [("cmdline", cmdline_), ("cmd2", "meglm"), ("cmd", "gsem"), ("title", sp["title"]),
                      ("family", family), ("model", sp["model"]), ("predict", "meglm_p"),
                      ("estat_cmd", "meglm_estat")] + common + [("link", link)] + tail
        else:
            macros = [("cmd", "meglm"), ("estat_cmd", f"{cmdname}_estat"), ("predict", f"{cmdname}_p"),
                      ("cmdline", cmdline_), ("cmd2", cmdname), ("model", sp["model"]), ("title", sp["title"]),
                      ("dispersion", sp["disp"]), ("link", link), ("family", family)] + common + tail
            if cmdname == "meoprobit":      # compat 1902
                d_ = dict(macros)
                macros = [(k_, d_[k_]) for k_ in ("cmd", "cmdline", "cmd2", "estat_cmd", "predict", "model",
                                                  "title", "link", "family")] + \
                         [m for m in macros if m[0] not in ("cmd", "cmdline", "cmd2", "estat_cmd", "predict",
                                                            "model", "title", "link", "family")]
        extra_mats = list(gmats)
        if ordinal:
            extra_mats.append(("cat", Matrix(np.array([F.cats], dtype=float), ["r1"],
                                             [f"c{i + 1}" for i in range(len(F.cats))])))
        if robust:
            extra_mats.append(("V_modelbased", Matrix(V_model_all, names_all, names_all)))
        extra_mats += [("gradient", Matrix(grad, ["r1"], names_all)),
                       ("ilog", Matrix(ilog(r), ["r1"], [f"c{i + 1}" for i in range(20)]))]
    if F.macro_edit is not None:
        macros = F.macro_edit(macros)
    post(s, est, scalars, macros, smp.mask, extra_mats=extra_mats)
    if F.late_error:
        raise StataError(198, F.late_error)
    if qr:
        _display_qr(s, est, header=not o.get("noheader"), group=not o.get("nogroup"),
                    lrtest=not o.get("nolrtest"), table=not o.get("nofetable"),
                    retable=not o.get("noretable"))
    else:
        display_me(s, est, header=not o.get("noheader"), group=not o.get("nogroup"),
                   lrtest=not o.get("nolrtest"))


def display_me(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
               level: float | None = None, group: bool = True, lrtest: bool = True, eform=None, **_kw) -> None:
    out = s.output
    x = est.extra
    sp = x["spec"]
    lev = level or est.level
    if header:
        ll = f"Log {'pseudolikelihood' if x.get('robust') else 'likelihood'} = {_g(x['ll'], '%10.0g'):>10}"
        meth = x["method"]
        mid = []
        if sp["cmd"] == "meglm":
            mid = [f"{'Family:':<7}{_FAM_TITLE[sp['fam']]:>24}", f"{'Link:':<5}{_LINK_TITLE[sp['link']]:>26}"]
        if sp["fam"] == "nbinomial":
            mid.append(f"{'Overdispersion:':<15}{sp['disp']:>16}")      # VERIFY
        if sp["fam"] == "binomial":
            nb_ = str(sp["nbin"])
            lab_ = "Binomial trials:" if nb_.replace(".", "").isdigit() else "Binomial variable:"
            mid.append(f"{lab_}{nb_:>{31 - len(lab_)}}")
        meth_line = (f"{'Integration method:':<19}{meth:>12}", "", "") if meth == "laplace" else \
            (f"{'Integration method:':<19}{meth:>12}", "Integration pts.", f"{int(x['npts'])}")
        group_header(out, sp["title"], x["groups"] if group else [], x["N"],
                     [("", "", ""),
                      meth_line,
                      ("", "", ""),
                      ("", f"Wald chi2({int(x['df_m'])})", _g(x["chi2"], "%10.2f")),
                      (ll, "Prob > chi2", _g(x["p"], "%10.4f"))], gv_right=True, mid=mid)
    if x.get("robust") and header:
        note = f"(Std. Err. adjusted for {int(x['n_clust']):,} clusters in {x['clustvar']})"   # VERIFY
        out.write(f"{note:>78}\n", "text")
    if not table:
        return
    z = st.norm.ppf(1 - (1 - lev / 100) / 2)
    ef = est.eform_default if eform is None else eform
    rows = list(est.rows)
    G = est.glmm
    from .mixed import _vc_rows
    sd = bool(x.get("sd"))
    k0 = 0
    for li, lv in enumerate(G["levels"]):
        npar = lv.npar()
        t_ = G["theta"][k0:k0 + npar]
        V_ = G["Vth"][k0:k0 + npar, k0:k0 + npar]
        k0 += npar
        # the label of a nested level joins the levels above it: school>class (VERIFY)
        eqlab = ">".join(l_.var for l_ in G["levels"][:li + 1])
        vc = _vc_rows(lv, t_, V_, sd, z)
        rows.append(CoefRow("sep"))
        rows.append(CoefRow("eq", eqlab))
        for lab, v, se, lo, hi in vc:
            if not lab.startswith("cov("):      # stddeviations: sd() and corr() rows, VERIFY the layout
                rows.append(CoefRow("aux", lab, v, se, eq=f"{float(lo)!r} {float(hi)!r}", noeform=True,
                                    tight=True))
        covs = [r_ for r_ in vc if r_[0].startswith("cov(")]
        if covs:
            # covariances: a section of their own, with z tests, and the pair written in reverse
            # order of e(revars) (compat 1002: cov(_cons,x) after var(x) and var(_cons))
            rows.append(CoefRow("sep"))
            rows.append(CoefRow("eq", eqlab))
            for lab, v, se, _lo, _hi in covs:
                a, b_ = lab[4:-1].split(",", 1) if "," in lab else (lab[4:-1], "")
                lab2 = f"cov({b_},{a})" if b_ else lab
                rows.append(CoefRow("coef", lab2, v, se, noeform=True, tight=True))
    if sp["fam"] == "gaussian":
        # residual variance, after the variances of the levels (VERIFY the layout)
        a_, se_a = float(G["aux"][0]), float(np.sqrt(max(G["Va"][0, 0], 0)))
        h_ = 0.5 if sd else 1.0
        v_ = math.exp(h_ * a_)
        rows.append(CoefRow("sep"))
        rows.append(CoefRow("aux", f"{'sd' if sd else 'var'}(e.{est.depvar})", v_, v_ * h_ * se_a,
                            eq=f"{math.exp(h_ * (a_ - z * se_a))!r} {math.exp(h_ * (a_ + z * se_a))!r}",
                            noeform=True, tight=True))
    coef_table(s, est.depvar, rows, stat="z", df=None, level=lev,
               coef_title=ef or "Coef.", eform=bool(ef), vcetype="Robust" if x.get("robust") else "")
    if lrtest and not x.get("robust"):
        lr_line(out, sp["vs"], x["chi2_c"], x["nvc"])


_VT = {"identity": "Identity", "independent": "Independent", "exchangeable": "Exchangeable",
       "unstructured": "Unstructured"}


def _display_qr(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                level: float | None = None, group: bool = True, lrtest: bool = True, eform=None,
                retable: bool = True, **_kw) -> None:
    """meqrlogit and meqrpoisson: the layout of mixed (header with the integration points, table
    of random-effects parameters); VERIFY against the Stata run."""
    from .mixed import _re_row
    out = s.output
    x = est.extra
    sp = x["spec"]
    lev = level or est.level
    if header:
        ll = f"Log likelihood = {_g(x['ll'], '%10.0g'):>10}"
        if len(x["groups"]) > 1 and group:
            # nested levels: the group table has a column of integration points (compat 1916)
            out.write("\n", "text")
            _right(out, sp["title"], "Number of obs", f"{int(x['N']):,}")
            out.write("\n" + "-" * 76 + "\n", "text")
            out.write("                |     No. of       Observations per Group       Integration\n", "text")
            out.write(" Group Variable |     Groups    Minimum    Average    Maximum      Points\n", "text")
            out.write("----------------+" + "-" * 59 + "\n", "text")
            for var, ng, mn, av, mx in x["groups"]:
                out.write(f"{var[:14]:>15} |", "text")
                out.write(f" {ng:>10,} {mn:>10,} {av:>10.1f} {mx:>10,} {int(x['npts']):>11}\n", "result")
            out.write("-" * 76 + "\n", "text")
            for left, lab_, val in [("", "", ""), ("", f"Wald chi2({int(x['df_m'])})", _g(x["chi2"], "%10.2f")),
                                    (ll, "Prob > chi2", _g(x["p"], "%10.4f"))]:
                if lab_:
                    _right(out, left, lab_, val)
                else:
                    out.write(left + "\n", "text")
        else:
            group_header(out, sp["title"], x["groups"] if group else [], x["N"],
                         [("", "", ""),
                          (f"Integration points = {int(x['npts']):>3}", f"Wald chi2({int(x['df_m'])})",
                           _g(x["chi2"], "%10.2f")),
                          (ll, "Prob > chi2", _g(x["p"], "%10.4f"))])
    if table:
        ef = est.eform_default if eform is None else eform
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=lev,
                   coef_title=ef or "Coef.", eform=bool(ef))
    z = st.norm.ppf(1 - (1 - lev / 100) / 2)
    G = est.glmm
    if retable:
        from .mixed import _vc_rows
        out.write("\n" + "-" * 78 + "\n", "text")
        out.write(f"  Random-effects Parameters  |   Estimate   Std. Err.     [{lev:g}% Conf. Interval]\n", "text")
        k0 = 0
        for lv in G["levels"]:
            k1 = k0 + lv.npar()
            out.write("-" * 29 + "+" + "-" * 48 + "\n", "text")
            kind = {"identity": "Identity", "independent": "Independent", "exchangeable": "Exchangeable",
                    "unstructured": "Unstructured"}[lv.cov]
            out.write(f"{lv.var + ': ' + kind:<29}|\n", "text")
            for row in _vc_rows(lv, G["theta"][k0:k1], G["Vth"][k0:k1, k0:k1], x["sd"], z):
                _re_row(out, *row)
            k0 = k1
        out.write("-" * 78 + "\n", "text")
    if lrtest:
        lr_line(out, sp["vs"], x["chi2_c"], x["nvc"])
    if x["method"] == "laplace":      # compat 1902
        out.write("\nNote: Log-likelihood calculations are based on the Laplacian approximation.\n", "text")


def _posterior_eta(s: "Session", est: Estimates, fixedonly: bool = False):
    """Linear predictor for every observation of the data: fixed part (with the offset) plus the
    posterior means of the random effects of the group (zero when the group is not in the
    estimation sample)."""
    from .postest import xb_all
    xb, ok = xb_all(s, est)
    ds = s.data
    G = est.glmm
    model = G["model"]
    if est.offset_name:
        nm = est.offset_name
        v = ds.get(nm[3:-1] if nm.startswith("ln(") else nm).data.astype(float)
        xb = xb + np.where(v < SYS, np.log(np.where(v > 0, v, 1)) if nm.startswith("ln(") else v, 0)
    if fixedonly:
        return xb, ok
    b = np.concatenate([G["beta"], G["aux"], G["theta"]])
    model.adapt(b, steps=50)
    eta = xb.copy()
    for li, lv in enumerate(G["levels"]):
        gidx = _group_index(s, G, li)
        re_all = np.where(gidx[:, None] >= 0, model.mus[li][np.maximum(gidx, 0)], np.nan)
        Zcols = []
        for nm in lv.names:
            Zcols.append(np.ones(ds.nobs) if nm == "_cons" else ds.get(nm).data.astype(float))
        Z = np.column_stack(Zcols)
        eta = eta + np.where(np.isfinite(re_all).all(axis=1), np.einsum("nq,nq->n", Z, np.nan_to_num(re_all)), 0.0)
    return eta, ok


def _group_index(s: "Session", G: dict, li: int) -> np.ndarray:
    """The group of level li of every observation of the data (-1 when the group is not in the
    estimation sample); a nested group is the combination of the values of the levels above."""
    ds = s.data
    rows = np.nonzero(G["mask"])[0]
    lvs = G["levels"][:li + 1]
    cols = [ds.get(lv.var).data for lv in lvs]
    keys = list(zip(*cols))
    key_of = {keys[i]: g for i, g in zip(rows, G["levels"][li].gid)}
    return np.array([key_of.get(k, -1) for k in keys])


def _predict_me(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    sp = est.extra["spec"]
    G = est.glmm
    cond = str((getattr(est, "_pred_opts", None) or {}).get("conditional") or "").strip().lower()
    if cond and not ("fixedonly".startswith(cond) or "ebmeans".startswith(cond)):
        raise StataError(198, f"conditional({cond}) invalid")   # ebmodes: VERIFY, not available
    if cond and "fixedonly".startswith(cond):
        stat = "fixedonly" if stat in ("", "mu") else stat
        if stat == "eta":
            stat = "xb"
    if stat == "xb":
        eta, ok = _posterior_eta(s, est, fixedonly=True)
        return np.where(ok, eta, SYS), ""
    eta, ok = _posterior_eta(s, est, fixedonly=stat == "fixedonly")
    if stat == "eta":
        if est.extra.get("qr") or (cond and "fixedonly".startswith(cond)):
            return np.where(ok, eta, SYS), ""
        return np.where(ok, eta, SYS), f"(predictions based on fixed effects and posterior means of random effects)\n(using {int(est.extra['npts'])} quadrature points)"
    mu = _mean_fn(sp)(eta)
    if sp["fam"] == "binomial":
        t = sp["nbin"]
        try:
            mu = mu * float(t)
        except ValueError:
            mu = mu * s.data.get(t).data.astype(float)
    # notes (compat 1902): fixedonly alone says only that mu was assumed; the meqr commands say
    # "predicted means"; conditional(fixedonly) is silent
    base = "(predictions based on fixed effects and posterior means of random effects)"
    quad = f"(using {int(est.extra['npts'])} quadrature points)"
    lines = []
    if est.extra.get("qr"):
        if stat in ("", "mu"):
            lines = ["(option mu assumed; predicted means)"] if not stat else []
    elif cond and "fixedonly".startswith(cond):
        lines = []
    elif stat == "fixedonly":
        lines = ["(option mu assumed)"]
    else:
        lines = [base] + (["(option mu assumed)"] if stat in ("", "outcome") else []) + [quad]
    note = "\n".join(lines)
    return np.where(ok, mu, SYS), note


def _pr_multi(s: "Session", est: Estimates):
    """predict, pr after meologit/meoprobit: the probability of each outcome (or of outcome(#))
    given the fixed effects and the posterior means of the random effects (VERIFY)."""
    G = est.glmm
    sp = est.extra["spec"]
    J = G["nlev"]
    oc = getattr(est, "_pred_opts", {}).get("outcome")
    nreq = 1        # compat 1902: pr takes one new variable ("too many variables specified" otherwise)
    s.output.write("(predictions based on fixed effects and posterior means of random effects)\n", "text")
    s.output.write(f"(using {int(est.extra['npts'])} quadrature points)\n", "text")

    def fn(mask):
        eta, ok = _posterior_eta(s, est)
        cuts = np.concatenate([[-np.inf], G["aux"], [np.inf]])
        F = special.expit if sp["link"] == "logit" else special.ndtr
        ps = [np.where(ok, F(cuts[j + 1] - eta) - F(cuts[j] - eta), SYS) for j in range(J)]
        if oc:
            t = str(oc).strip().lstrip("#")
            k = int(float(t)) - 1 if str(oc).strip().startswith("#") else None
            if k is None:
                cats = [c for c in range(J)]       # outcome(value): VERIFY with non-consecutive values
                k = int(float(t)) - 1 if 0 <= int(float(t)) - 1 < len(cats) else 0
            return [ps[k]]
        return ps[:1]       # pr with one variable: the first outcome (VERIFY)
    return nreq, fn, "too few variables specified"


def _reffects_me(s: "Session", est: Estimates):
    """predict newvarlist, reffects: posterior means, one variable per random effect in the
    order of the level's effects (VERIFY: the order for random slopes)."""
    G = est.glmm
    nreq = sum(lv.q for lv in G["levels"])
    if not est.extra.get("qr"):
        s.output.write("(calculating posterior means of random effects)\n", "text")

    def fn(mask):
        model = G["model"]
        model.adapt(np.concatenate([G["beta"], G["aux"], G["theta"]]), steps=50)
        out = []
        # the highest level first, within a level the order of its effects (VERIFY)
        for li, lv in enumerate(G["levels"]):
            gidx = _group_index(s, G, li)
            for j in range(lv.q):
                vals = np.where(gidx >= 0, model.mus[li][np.maximum(gidx, 0), j], SYS)
                out.append(np.where(np.isfinite(vals), vals, SYS))
        if not est.extra.get("qr"):
            s.output.write(f"(using {int(est.extra['npts'])} quadrature points)\n", "text")
        return out
    # Stata's message here has a single line (compat 1002)
    return nreq, fn, "too few variables specified"


def cmd_melogit(s, args):
    cmd_me(s, args, "melogit")


def cmd_meprobit(s, args):
    cmd_me(s, args, "meprobit")


def cmd_mepoisson(s, args):
    cmd_me(s, args, "mepoisson")


def cmd_mecloglog(s, args):
    cmd_me(s, args, "mecloglog")


def cmd_menbreg(s, args):
    cmd_me(s, args, "menbreg")


def cmd_meologit(s, args):
    cmd_me(s, args, "meologit")


def cmd_meoprobit(s, args):
    cmd_me(s, args, "meoprobit")


def cmd_meglm(s, args):
    cmd_me(s, args, "meglm")


def cmd_meqrlogit(s, args):
    cmd_me(s, args, "meqrlogit")


def cmd_meqrpoisson(s, args):
    cmd_me(s, args, "meqrpoisson")

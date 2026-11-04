"""The likelihood of the two-level mixed-effects GLMs of meglm and its special cases: the
conditional distributions of the families and the group likelihood integrated by Gauss-Hermite
quadrature (adaptive, non-adaptive or Laplace), with analytic gradient.

GLMM holds a single outcome; MGLMM the multinomial and other models with several linear
predictors. Plain arrays in, log likelihood and derivatives out: no Session, no output.
"""

from __future__ import annotations

import math

import numpy as np

from ..core.errors import StataError
from ..lazy import lazy
from .mixed import Level

# imported on first use (start-up time; see opendta/lazy.py)
special = lazy("scipy.special")
st = lazy("scipy.stats")


def basic_family(kind: str, y: np.ndarray):
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


def gh_nodes(n: int, q: int):
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
        self.fam = fam if fam is not None else basic_family(kind, y)
        if method == "laplace":
            npoints = 1
        self.nodes, self.lw = gh_nodes(npoints, self.q)
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
            nd, lw = gh_nodes(npoints, lv.q)
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

    # the attributes the single-level model exposes to fit_me
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

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
    build_levels,
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
    if fam == "gaussian":
        def f(eta, a):
            v = math.exp(a[0])
            r = yc - eta
            return -0.5 * (ln2pi + a[0]) - r * r / (2 * v), r / v, np.full(eta.shape, -1.0 / v)
        return f, 1
    if fam == "gamma":
        lny = np.log(yc)

        def f(eta, a):
            phi = math.exp(a[0])
            k = 1.0 / phi
            r = yc * np.exp(-np.clip(eta, -700, 700)) / phi
            return (k - 1) * lny - r - k * (a[0] + eta) - special.gammaln(k), r - k, -r
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


class GLMM:
    def __init__(self, kind, y, X, offset, level: Level, npoints: int, method: str, fam=None,
                 naux: int = 0, ordered: bool = False):
        self.kind = kind
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
        ll = float(r[0].sum())
        return ll if np.isfinite(ll) else -np.inf

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
            ok = np.all(np.linalg.eigvalsh(cov) > 1e-12, axis=1)
            newC = self.C.copy()
            newC[ok] = np.linalg.cholesky(cov[ok])
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
        pw = np.exp(lg - lL[:, None])
        pwo = pw[self.gid]                                      # N x K
        gb = self.X.T @ np.sum(pwo * s1, axis=1)
        ga = np.zeros(self.naux)
        a = self.aux(b)
        for m in range(self.naux):
            h = 1e-5 * (abs(a[m]) + 1)
            e = np.zeros(self.naux)
            e[m] = h
            ga[m] = float(np.sum(pwo * (self.fam(eta, a + e)[0] - self.fam(eta, a - e)[0])) / (2 * h))
        # variance parameters: d/dth ln phi_S(u) = -1/2 tr(S^-1 dS) + 1/2 u' S^-1 dS S^-1 u
        beta, th = self.split(b)
        gth = np.zeros(th.size)
        h = 1e-6
        for m in range(th.size):
            e = np.zeros(th.size)
            e[m] = h
            dS = (self.lv.sigma(th + e) - self.lv.sigma(th - e)) / (2 * h)
            A = Sinv @ dS @ Sinv
            quad = np.einsum("gki,ij,gkj->gk", U, A, U)
            gth[m] = float(np.sum(pw * (-0.5 * np.trace(Sinv @ dS) + 0.5 * quad)))
        return np.concatenate([gb, ga, gth])

    def fun(self, b, todo):
        if todo == 0:
            return self.loglik(b), None, None
        self.adapt(b)
        ll = self.loglik(b)
        if not np.isfinite(ll):
            return -np.inf, None, None
        g = self.gradient(b)
        k = b.size
        H = np.zeros((k, k))
        for i in range(k):
            hi = 1e-5 * (abs(b[i]) + 1)
            e = np.zeros(k)
            e[i] = hi
            H[:, i] = (self.gradient(b + e) - self.gradient(b - e)) / (2 * hi)
        H = (H + H.T) / 2
        return ll, g, H


# ---------------------------------------------------------------------------
# the commands
# ---------------------------------------------------------------------------

# command -> (family, link, title, model, reference model of the LR test)
_CMDS = {
    "melogit": ("bernoulli", "logit", "Mixed-effects logistic regression", "logistic", "logistic model"),
    "meprobit": ("bernoulli", "probit", "Mixed-effects probit regression", "probit", "probit model"),
    "mecloglog": ("bernoulli", "cloglog", "Mixed-effects cloglog regression", "cloglog",
                  "complementary log-log model"),
    "mepoisson": ("poisson", "log", "Mixed-effects Poisson regression", "Poisson", "Poisson model"),
    "menbreg": ("nbinomial", "log", "Mixed-effects nbinomial regression", "nbreg", "nbinomial model"),
    "meologit": ("ordinal", "logit", "Mixed-effects ologit regression", "ologit", "ologit model"),
    "meoprobit": ("ordinal", "probit", "Mixed-effects oprobit regression", "oprobit", "oprobit model"),
    "meglm": (None, None, "Mixed-effects GLM", "", ""),
}
_FAMILIES = ("gaussian", "bernoulli", "binomial", "gamma", "nbinomial", "ordinal", "poisson")
_DEFAULT_LINK = {"gaussian": "identity", "bernoulli": "logit", "binomial": "logit", "gamma": "log",
                 "nbinomial": "log", "ordinal": "logit", "poisson": "log"}
_LINKS_OK = {"gaussian": ("identity",), "bernoulli": ("logit", "probit", "cloglog"),
             "binomial": ("logit", "probit", "cloglog"), "gamma": ("log",), "nbinomial": ("log",),
             "ordinal": ("logit", "probit"), "poisson": ("log",)}
# reference model of the LR test of meglm (VERIFY the names)
_VS_GLM = {("gaussian", "identity"): "linear model", ("bernoulli", "logit"): "logistic model",
           ("bernoulli", "probit"): "probit model", ("bernoulli", "cloglog"): "cloglog model",
           ("binomial", "logit"): "logistic model", ("binomial", "probit"): "probit model",
           ("binomial", "cloglog"): "cloglog model", ("gamma", "log"): "gamma model",
           ("nbinomial", "log"): "nbinomial model", ("ordinal", "logit"): "ologit model",
           ("ordinal", "probit"): "oprobit model", ("poisson", "log"): "Poisson model"}
_FAM_TITLE = {"gaussian": "Gaussian", "bernoulli": "Bernoulli", "binomial": "Binomial",
              "gamma": "Gamma", "nbinomial": "Negative binomial", "ordinal": "Ordinal",
              "poisson": "Poisson"}                                  # VERIFY
_LINK_TITLE = {"identity": "Identity", "logit": "Logit", "probit": "Probit", "log": "Log",
               "cloglog": "Complementary log-log"}                   # VERIFY
_GLOBAL = {"or": 2, "irr": 3, "intpoints": 4, "intmethod": 7, "nolog": 4, "log": 3, "iterate": 4,
           "tolerance": 3, "ltolerance": 4, "nrtolerance": 4, "level": 1, "noheader": 6, "nogroup": 5,
           "nolrtest": 5, "noconstant": 4, "offset": 3, "exposure": 3, "binomial": 3, "covariance": 3,
           "startgrid": 6, "nofetable": 5, "noretable": 5, "estmetric": 4, "family": 3, "link": 4,
           "dispersion": 4, "asis": 4}


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
        a0 = [math.log(max(float(y.var()) / max(ybar * ybar, 1e-10), 1e-2))]
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


def cmd_me(s: "Session", args: str, cmdname: str) -> None:
    from .postest import replay
    from .regress import Sample, cmdline
    if not args.strip() or args.strip().startswith(","):
        replay(s, cmdname, args)
        return
    fixed, levels_spec, glob = split_levels(args)
    if len(levels_spec) != 1:
        raise StataError(198, "OpenDTA's me commands support one level of random effects")   # OpenDTA
    c = find_top(fixed, ",")
    fixed_body, fe_opts = (fixed, "") if c == -1 else (fixed[:c], fixed[c + 1:])
    opts_text = " ".join(x for x in (fe_opts.strip(), glob.strip()) if x)
    extra = " ".join([lv for lv, _r, _o in levels_spec] + [r for _l, r, _o in levels_spec if r])
    smp = Sample(s, fixed_body + (f", {opts_text}" if opts_text else ""), opts_spec=_GLOBAL,
                 weights=(), extra_vars=extra)
    o = smp.o
    sp = _resolve(cmdname, o)
    fam, link = sp["fam"], sp["link"]
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
    lvspec = levels_spec[0]
    if lvspec[2].get("covariance") is None and o.get("covariance"):
        lvspec[2]["covariance"] = o["covariance"]
    levels = build_levels(s, smp.mask, [lvspec], n)
    lv = levels[0]
    if lv.cov == "identity" and lv.q > 1 and not lvspec[2].get("covariance"):
        lv.cov = "independent"
    npts = int(o["intpoints"]) if o.get("intpoints") else 7
    method = str(o.get("intmethod") or "mvaghermite").strip().lower()
    method = {"mv": "mvaghermite", "mvaghermite": "mvaghermite", "mc": "mcaghermite",
              "mcaghermite": "mcaghermite", "gh": "ghermite", "ghermite": "ghermite",
              "laplace": "laplace"}.get(method, method)
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
    model = GLMM(fam, y, X, off, lv, npts, method, fam=fam_f, naux=naux, ordered=fam == "ordinal")
    kw = ml_options(o)
    p = X.shape[1]
    # fixed-effects model
    if naux == 0 and fam != "binomial" and not (fam == "bernoulli" and link == "cloglog"):
        def fe_fun(b, todo):
            eta = (X @ b + off)[:, None]
            lf, s1, s2 = fam_f(eta)
            ll = float(lf.sum())
            if todo == 0:
                return ll, None, None
            return ll, X.T @ s1[:, 0], (X * s2[:, 0][:, None]).T @ X
        b0 = _fixed_start(sp, y, off, p, smp.constant, nbin)[0]
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
            if fam == "ordinal" and np.any(np.diff(b[p:]) <= 0):
                return -np.inf, None, None
            eta = (X @ b[:p] + off)[:, None]
            ll = float(fam_f(eta, b[p:])[0].sum())
            if todo == 0:
                return ll, None, None
            return ll, fe_grad(b), _num_gh(fe_grad, b)
        if fam == "ordinal":
            cum = np.array([np.mean(yi <= j) for j in range(len(cats) - 1)])
            b0 = np.concatenate([np.zeros(p), (special.logit if link == "logit" else st.norm.ppf)(cum)])
        else:
            b0b, a0 = _fixed_start(sp, y, off, p, smp.constant, nbin)
            if fam == "gaussian" and p:
                b0b = np.linalg.lstsq(X, y - off, rcond=None)[0]
                a0 = np.array([math.log(max(float(np.mean((y - off - X @ b0b) ** 2)), 1e-8))])
            b0 = np.concatenate([b0b, a0])
    rfe = maximize(fe_fun, b0, **kw)
    # starting values for the variance components: best of a small grid (VERIFY)
    best = None
    for ls in (math.log(0.5), 0.0, math.log(0.25), math.log(1.5)):
        th = np.zeros(lv.npar())
        th[:min(lv.q, th.size) if lv.cov != "exchangeable" else 1] = ls
        model.C = None
        model.mu[:] = 0
        bb = np.concatenate([rfe.b, th])
        model.adapt(bb)
        ll = model.loglik(bb)
        if best is None or ll > best[0]:
            best = (ll, bb, model.mu.copy(), model.C.copy())
    model.mu, model.C = best[2], best[3]
    lab = "log likelihood"
    if not o.get("nolog"):
        out.write("\nFitting fixed-effects model:\n\n", "text")
        print_log(out, rfe, lab)
        out.write("\nRefining starting values:\n\n", "text")
        out.write(f"{'Grid node 0:':<14} {lab} = {_g(best[0], '%10.0g'):>10}\n", "text")
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
    th_names = [f"{t}:_cons" for t in lv.labels()]
    # auxiliary parameters: names in e(b) and the rows of the table
    if fam == "ordinal":
        aux_eq = [f"/cut{j + 1}" for j in range(naux)]
    elif fam == "nbinomial":
        aux_eq = ["/lnalpha" if sp["disp"] == "mean" else "/lndelta"]
    elif fam == "gamma":
        aux_eq = ["/lnphi"]                                     # VERIFY
    elif fam == "gaussian":
        aux_eq = [f"/var(e.{smp.depname})"]
    else:
        aux_eq = []
    ball = np.concatenate([bf, aux, th])
    Vall = np.zeros((ball.size, ball.size))
    Vall[:kk, :kk] = Vf
    Vall[kk:kk + naux, kk:kk + naux] = Va
    Vall[kk:kk + naux, kk + naux:] = V[p:p + naux, p + naux:]
    Vall[kk + naux:, kk:kk + naux] = V[p + naux:, p:p + naux]
    Vall[kk + naux:, kk + naux:] = Vth
    est = Estimates(cmdname, smp.depname, names + ["_cons"] * (naux + th.size), ball, Vall, n, "z", None)
    est.eqnames = [smp.depname] * kk + aux_eq + [f"/{t.split(':')[0]}" for t in th_names]
    est.rows = rows_from_columns(smp.cols, bmap, semap, ds, constant=smp.constant)
    if off_name:
        est.rows.append(CoefRow("fixed", off_name, eq="(exposure)" if o.get("exposure") else "(offset)"))
    if fam in ("ordinal", "nbinomial", "gamma"):
        est.rows.append(CoefRow("sep"))
        for j, nm_ in enumerate(aux_eq):
            se_a = float(np.sqrt(max(Va[j, j], 0)))
            est.rows.append(CoefRow("aux", nm_, float(aux[j]), se_a, name="_cons"))
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = ""
    est.robust = False
    est.offset_name = off_name
    resid = {("bernoulli", "logit"): math.pi ** 2 / 3, ("bernoulli", "probit"): 1.0,
             ("bernoulli", "cloglog"): math.pi ** 2 / 6, ("binomial", "logit"): math.pi ** 2 / 3,
             ("binomial", "probit"): 1.0, ("binomial", "cloglog"): math.pi ** 2 / 6,
             ("ordinal", "logit"): math.pi ** 2 / 3, ("ordinal", "probit"): 1.0}.get((fam, link))
    if fam == "gaussian":
        resid = float(math.exp(aux[0]))
    est.glmm = {"model": model, "theta": th, "Vth": Vth, "beta": beta, "levels": [lv], "aux": aux,
                "Va": Va, "aux_eq": aux_eq, "resid_var": resid, "mask": smp.mask.copy(), "nbin": nbin,
                "nlev": (len(cats) if fam == "ordinal" else 0)}
    cnt = np.bincount(lv.gid)
    est.extra = {"N": float(n), "ll": r.ll, "chi2": wchi2, "p": wp, "df_m": len(slopes), "ll_c": rfe.ll,
                 "chi2_c": chi2_c, "nvc": lv.npar(), "spec": sp, "method": method,
                 "groups": [(lv.var, len(cnt), int(cnt.min()), float(cnt.mean()), int(cnt.max()))]}
    est.extra["npts"] = 1 if method == "laplace" else npts
    est.display = display_me
    est.predict = _predict_me
    est.predict_opts = {"mu": 2, "eta": 3, "reffects": 4, "fixedonly": 5, "pr": 2, "outcome": 3}
    est.predict_multi = {"reffects": _reffects_me}
    if fam == "ordinal":
        est.predict_multi.update({"": _pr_multi, "pr": _pr_multi})
    est.estat = {"icc": (3, estat_icc), "recovariance": (5, estat_recov)}
    est.eform_default = "Odds Ratio" if (link == "logit" and o.get("or")) else \
        ("IRR" if (fam in ("poisson", "nbinomial") and o.get("irr")) else "")
    # e() in the order Stata lists it after melogit (compat 1002); VERIFY for the other commands
    nvar_p = {"identity": 1, "exchangeable": 1}.get(lv.cov, lv.q)
    ncov_p = lv.npar() - nvar_p
    scalars = [("rank", float(ball.size)), ("N", float(n)), ("ic", float(r.iterations)),
               ("converged", 1.0 if r.converged else 0.0), ("rc", 0.0), ("ll", r.ll),
               ("k", float(ball.size)), ("k_eq", float(1 + naux + lv.npar())), ("k_eq_model", 1.0),
               ("ll_c", rfe.ll), ("k_rs", float(nvar_p)), ("k_rc", float(ncov_p)), ("k_dv", 1.0),
               ("k_r", float(lv.q)), ("k_f", float(len(slopes))),   # VERIFY: k_f without _cons
               ("p_c", _p_chibar(chi2_c, lv.npar())), ("df_c", float(lv.npar())), ("chi2_c", chi2_c),
               ("p", wp), ("chi2", wchi2), ("df_m", float(len(slopes)))]
    from ..commands.matrix import Matrix
    from .fvars import base_vars
    from .ml import ilog
    from .results import datasignature
    covs = [c.name for c in smp.cols if c.values is not None]
    sigvars = [lv.var] + [v for v in base_vars(smp.terms) if v != smp.depname] + [smp.depname]
    family = {"bernoulli": "bernoulli", "binomial": "binomial", "poisson": "poisson",
              "nbinomial": "nbinomial", "gaussian": "gaussian", "gamma": "gamma",
              "ordinal": "ordinal"}[fam]
    macros = [("cmd", "meglm"), ("estat_cmd", f"{cmdname}_estat"), ("predict", f"{cmdname}_p"),
              ("cmdline", cmdline(cmdname, args)), ("cmd2", cmdname), ("model", sp["model"]),
              ("title", sp["title"]), ("link", link), ("family", family),
              ("dispersion", sp["disp"]),
              ("n_quad", str(int(est.extra["npts"]))), ("ivars", lv.var), ("chi2type", "Wald"),
              ("datasignaturevars", " ".join(sigvars)),
              ("datasignature", datasignature(s, sigvars, smp.mask)), ("intmethod", method),
              ("depvar", smp.depname), ("marginsnotok", "LATent1(passthru) LATent se(passthru)"),
              ("covariates", " ".join(covs)), ("opt", "moptimize"), ("vce", "oim"),
              ("user", "_gsem_eval_chol()"), ("ml_method", "gf2"), ("technique", "nr"),
              ("which", "max"), ("properties", "b V")]
    g1 = est.extra["groups"][0]
    names_all = list(est.names)
    grad = np.zeros((1, len(names_all)))
    if r.g is not None and len(r.g) == len(names_all):
        grad[0] = r.g
    extra_mats = [(nm, Matrix(np.array([[val]], dtype=float), ["r1"], [lv.var]))
                  for nm, val in (("g_max", g1[4]), ("g_avg", g1[3]), ("g_min", g1[2]), ("N_g", g1[1]))]
    extra_mats += [("gradient", Matrix(grad, ["r1"], names_all)),
                   ("ilog", Matrix(ilog(r), ["r1"], [f"c{i + 1}" for i in range(20)]))]
    post(s, est, scalars, macros, smp.mask, extra_mats=extra_mats)
    display_me(s, est, header=not o.get("noheader"), group=not o.get("nogroup"),
               lrtest=not o.get("nolrtest"))


def display_me(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
               level: float | None = None, group: bool = True, lrtest: bool = True, eform=None, **_kw) -> None:
    out = s.output
    x = est.extra
    sp = x["spec"]
    lev = level or est.level
    if header:
        ll = f"Log likelihood = {_g(x['ll'], '%10.0g'):>10}"
        meth = x["method"]
        mid = []
        if sp["cmd"] == "meglm":
            mid = [f"{'Family:':<15}{_FAM_TITLE[sp['fam']]:>16}", f"{'Link:':<15}{_LINK_TITLE[sp['link']]:>16}"]
        if sp["fam"] == "nbinomial":
            mid.append(f"{'Overdispersion:':<15}{sp['disp']:>16}")      # VERIFY
        group_header(out, sp["title"], x["groups"] if group else [], x["N"],
                     [("", "", ""),
                      (f"Integration method: {meth}", "Integration pts.", f"{int(x['npts'])}"),
                      ("", "", ""),
                      ("", f"Wald chi2({int(x['df_m'])})", _g(x["chi2"], "%10.2f")),
                      (ll, "Prob > chi2", _g(x["p"], "%10.4f"))], gv_right=True, mid=mid)
    if not table:
        return
    z = st.norm.ppf(1 - (1 - lev / 100) / 2)
    ef = est.eform_default if eform is None else eform
    rows = list(est.rows)
    G = est.glmm
    lv = G["levels"][0]
    from .mixed import _vc_rows
    vc = _vc_rows(lv, G["theta"], G["Vth"], False, z)
    rows.append(CoefRow("sep"))
    rows.append(CoefRow("eq", lv.var))
    for lab, v, se, lo, hi in vc:
        if not lab.startswith("cov("):
            rows.append(CoefRow("aux", lab, v, se, eq=f"{float(lo)!r} {float(hi)!r}", noeform=True,
                                tight=True))
    covs = [r_ for r_ in vc if r_[0].startswith("cov(")]
    if covs:
        # covariances: a section of their own, with z tests, and the pair written in reverse
        # order of e(revars) (compat 1002: cov(_cons,x) after var(x) and var(_cons))
        rows.append(CoefRow("sep"))
        rows.append(CoefRow("eq", lv.var))
        for lab, v, se, _lo, _hi in covs:
            a, b_ = lab[4:-1].split(",", 1) if "," in lab else (lab[4:-1], "")
            lab2 = f"cov({b_},{a})" if b_ else lab
            rows.append(CoefRow("coef", lab2, v, se, noeform=True, tight=True))
    if sp["fam"] == "gaussian":
        # residual variance, after the variances of the levels (VERIFY the layout)
        a_, se_a = float(G["aux"][0]), float(np.sqrt(max(G["Va"][0, 0], 0)))
        v_ = math.exp(a_)
        rows.append(CoefRow("sep"))
        rows.append(CoefRow("aux", f"var(e.{est.depvar})", v_, v_ * se_a,
                            eq=f"{math.exp(a_ - z * se_a)!r} {math.exp(a_ + z * se_a)!r}",
                            noeform=True, tight=True))
    coef_table(s, est.depvar, rows, stat="z", df=None, level=lev,
               coef_title=ef or "Coef.", eform=bool(ef))
    if lrtest:
        lr_line(out, sp["vs"], x["chi2_c"], x["nvc"])


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
    rows = np.nonzero(G["mask"])[0]
    lv = G["levels"][0]
    var = ds.get(lv.var).data
    key_of = {}
    for i, g in zip(rows, lv.gid):
        key_of[var[i]] = g
    gidx = np.array([key_of.get(v, -1) for v in var])
    re_all = np.where(gidx[:, None] >= 0, model.mu[np.maximum(gidx, 0)], np.nan)
    Zcols = []
    for nm in lv.names:
        Zcols.append(np.ones(ds.nobs) if nm == "_cons" else ds.get(nm).data.astype(float))
    Z = np.column_stack(Zcols)
    eta = xb + np.where(np.isfinite(re_all).all(axis=1), np.einsum("nq,nq->n", Z, np.nan_to_num(re_all)), 0.0)
    return eta, ok


def _predict_me(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    sp = est.extra["spec"]
    G = est.glmm
    if stat == "xb":
        eta, ok = _posterior_eta(s, est, fixedonly=True)
        return np.where(ok, eta, SYS), ""
    eta, ok = _posterior_eta(s, est, fixedonly=stat == "fixedonly")
    if stat == "eta":
        return np.where(ok, eta, SYS), ""
    mu = _mean_fn(sp)(eta)
    if sp["fam"] == "binomial":
        t = sp["nbin"]
        try:
            mu = mu * float(t)
        except ValueError:
            mu = mu * s.data.get(t).data.astype(float)
    note = "" if stat == "fixedonly" else \
        "(predictions based on fixed effects and posterior means of random effects)"
    if not stat:
        note += "\n(option mu assumed)"
    if stat != "fixedonly":
        note += f"\n(using {int(est.extra['npts'])} quadrature points)"
    return np.where(ok, mu, SYS), note


def _pr_multi(s: "Session", est: Estimates):
    """predict, pr after meologit/meoprobit: the probability of each outcome (or of outcome(#))
    given the fixed effects and the posterior means of the random effects (VERIFY)."""
    G = est.glmm
    sp = est.extra["spec"]
    J = G["nlev"]
    oc = getattr(est, "_pred_opts", {}).get("outcome")
    nreq = 1 if oc else J

    def fn(mask):
        eta, ok = _posterior_eta(s, est)
        cuts = np.concatenate([[-np.inf], G["aux"], [np.inf]])
        F = special.expit if sp["link"] == "logit" else special.ndtr
        ps = [np.where(ok, F(cuts[j + 1] - eta) - F(cuts[j] - eta), SYS) for j in range(J)]
        s.output.write("(predictions based on fixed effects and posterior means of random effects)\n", "text")
        s.output.write(f"(using {int(est.extra['npts'])} quadrature points)\n", "text")
        if oc:
            t = str(oc).strip().lstrip("#")
            k = int(float(t)) - 1 if str(oc).strip().startswith("#") else None
            if k is None:
                cats = [c for c in range(J)]       # outcome(value): VERIFY with non-consecutive values
                k = int(float(t)) - 1 if 0 <= int(float(t)) - 1 < len(cats) else 0
            return [ps[k]]
        return ps
    return nreq, fn, "too few variables specified"


def _reffects_me(s: "Session", est: Estimates):
    """predict newvarlist, reffects: posterior means, one variable per random effect in the
    order of the level's effects (VERIFY: the order for random slopes)."""
    G = est.glmm
    lv = G["levels"][0]
    s.output.write("(calculating posterior means of random effects)\n", "text")

    def fn(mask):
        ds = s.data
        model = G["model"]
        model.adapt(np.concatenate([G["beta"], G["aux"], G["theta"]]), steps=50)
        rows = np.nonzero(G["mask"])[0]
        var = ds.get(lv.var).data
        key_of = {var[i]: g for i, g in zip(rows, lv.gid)}
        gidx = np.array([key_of.get(v, -1) for v in var])
        out = []
        for j in range(lv.q):
            vals = np.where(gidx >= 0, model.mu[np.maximum(gidx, 0), j], SYS)
            out.append(np.where(np.isfinite(vals), vals, SYS))
        s.output.write(f"(using {int(est.extra['npts'])} quadrature points)\n", "text")
        return out
    # Stata's message here has a single line (compat 1002)
    return lv.q, fn, "too few variables specified"


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

"""Mixed-effects generalized linear models: melogit, meprobit and mepoisson
([ME] melogit, [ME] meprobit, [ME] mepoisson, [ME] meglm).

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
    with y broadcast on the first axis)."""
    yc = y[:, None]
    if kind == "logit":
        def f(eta):
            lf = yc * eta - np.logaddexp(0, eta)
            p = special.expit(eta)
            return lf, yc - p, -p * (1 - p)
    elif kind == "probit":
        q = 2 * yc - 1

        def f(eta):
            lc = special.log_ndtr(q * eta)
            lam = q * np.exp(st.norm.logpdf(q * eta) - lc)
            return lc, lam, -lam * (lam + eta)
    else:
        lgy = special.gammaln(yc + 1)

        def f(eta):
            mu = np.exp(np.clip(eta, -700, 700))
            return yc * eta - mu - lgy, yc - mu, -mu
    return f


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
    def __init__(self, kind, y, X, offset, level: Level, npoints: int, method: str):
        self.kind = kind
        self.y, self.X, self.off = y, X, offset
        self.lv = level
        self.q = level.q
        self.gid = level.gid
        self.G = int(self.gid.max()) + 1
        self.Z = level.Z
        self.fam = _family(kind, y)
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
        return b[:p], b[p:]

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
        lf, s1, s2 = self.fam(eta)
        K = self.nodes.shape[0]
        sumlf = np.zeros((self.G, K))
        np.add.at(sumlf, self.gid, lf)
        quad = np.einsum("gki,ij,gkj->gk", U, Sinv, U)
        lprior = -0.5 * (self.q * math.log(2 * math.pi) + ldS + quad)
        _s2, ldC = np.linalg.slogdet(self.C)
        lg = sumlf + lprior - self.lphi[None, :] + ldC[:, None] + self.lw[None, :]
        lL = special.logsumexp(lg, axis=1)
        return lL, lg, s1, s2, U, Sinv, S

    def loglik(self, b) -> float:
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
            _lf, s1, s2 = self.fam(eta[:, None])
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
        lL, lg, s1, s2, U, Sinv, S = r
        pw = np.exp(lg - lL[:, None])
        pwo = pw[self.gid]                                      # N x K
        gb = self.X.T @ np.sum(pwo * s1, axis=1)
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
        return np.concatenate([gb, gth])

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

_TITLE = {"logit": "Mixed-effects logistic regression", "probit": "Mixed-effects probit regression",
          "poisson": "Mixed-effects Poisson regression"}
_VS = {"logit": "logistic model", "probit": "probit model", "poisson": "Poisson model"}
_GLOBAL = {"or": 2, "irr": 3, "intpoints": 4, "intmethod": 7, "nolog": 4, "log": 3, "iterate": 4,
           "tolerance": 3, "ltolerance": 4, "nrtolerance": 4, "level": 1, "noheader": 6, "nogroup": 5,
           "nolrtest": 5, "noconstant": 4, "offset": 3, "exposure": 3, "binomial": 3, "covariance": 3,
           "startgrid": 6, "nofetable": 5, "noretable": 5, "estmetric": 4}


def cmd_me(s: "Session", args: str, kind: str, cmdname: str) -> None:
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
    out = s.output
    ds = s.data
    y = smp.y.astype(np.float64)
    if kind in ("logit", "probit"):
        y = (y != 0).astype(float)
    elif np.any(y < 0):
        raise StataError(459, f"{smp.depname} must be greater than or equal to zero")   # VERIFY
    fvars.drop_collinear(smp.cols, smp.constant)
    est_cols = [col for col in smp.cols if col.values is not None]
    n = smp.n
    X = np.column_stack([col.values for col in est_cols] + ([np.ones(n)] if smp.constant else [])) \
        if (est_cols or smp.constant) else np.zeros((n, 0))
    off = np.zeros(n)
    off_name = ""
    if kind == "poisson" and (o.get("offset") or o.get("exposure")):
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
    model = GLMM(kind, y, X, off, lv, npts, method)
    kw = dict(maxiter=int(o["iterate"]) if o.get("iterate") else 16000,
              tol=float(o["tolerance"]) if o.get("tolerance") else 1e-6,
              ltol=float(o["ltolerance"]) if o.get("ltolerance") else 1e-7,
              nrtol=float(o["nrtolerance"]) if o.get("nrtolerance") else 1e-5)
    # fixed-effects model
    p = X.shape[1]
    fam = _family(kind, y)

    def fe_fun(b, todo):
        eta = (X @ b + off)[:, None]
        lf, s1, s2 = fam(eta)
        ll = float(lf.sum())
        if todo == 0:
            return ll, None, None
        return ll, X.T @ s1[:, 0], (X * s2[:, 0][:, None]).T @ X
    b0 = np.zeros(p)
    if smp.constant:
        ybar = float(y.mean())
        if kind == "logit":
            b0[-1] = math.log(ybar / (1 - ybar))
        elif kind == "probit":
            b0[-1] = float(st.norm.ppf(ybar))
        else:
            b0[-1] = math.log(max(ybar, 1e-10)) - float(off.mean())
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
    Vb = V[:p, :p]
    Vth = V[p:, p:]
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
    ball = np.concatenate([bf, th])
    Vall = np.zeros((ball.size, ball.size))
    Vall[:kk, :kk] = Vf
    Vall[kk:, kk:] = Vth
    est = Estimates(cmdname, smp.depname, names + ["_cons"] * th.size, ball, Vall, n, "z", None)
    est.eqnames = [smp.depname] * kk + [f"/{t.split(':')[0]}" for t in th_names]
    est.rows = rows_from_columns(smp.cols, bmap, semap, ds, constant=smp.constant)
    if off_name:
        est.rows.append(CoefRow("fixed", off_name, eq="(exposure)" if o.get("exposure") else "(offset)"))
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = ""
    est.robust = False
    est.offset_name = off_name
    est.glmm = {"model": model, "theta": th, "Vth": Vth, "beta": beta, "levels": [lv],
                "resid_var": {"logit": math.pi ** 2 / 3, "probit": 1.0}.get(kind), "mask": smp.mask.copy()}
    cnt = np.bincount(lv.gid)
    est.extra = {"N": float(n), "ll": r.ll, "chi2": wchi2, "p": wp, "df_m": len(slopes), "ll_c": rfe.ll,
                 "chi2_c": chi2_c, "nvc": lv.npar(), "kind": kind, "npts": model.nodes.shape[0] ** (1 / lv.q)
                 if lv.q else npts, "method": method,
                 "groups": [(lv.var, len(cnt), int(cnt.min()), float(cnt.mean()), int(cnt.max()))]}
    est.extra["npts"] = 1 if method == "laplace" else npts
    est.display = display_me
    est.predict = _predict_me
    est.predict_opts = {"mu": 2, "eta": 3, "reffects": 4, "fixedonly": 5, "pr": 2}
    est.predict_multi = {"reffects": _reffects_me}
    est.estat = {"icc": (3, estat_icc), "recovariance": (5, estat_recov)}
    est.eform_default = "Odds Ratio" if (kind == "logit" and o.get("or")) else \
        ("IRR" if (kind == "poisson" and o.get("irr")) else "")
    scalars = [("N", float(n)), ("k", float(ball.size)), ("k_f", float(kk)), ("k_r", float(th.size)),
               ("ll", r.ll), ("chi2", wchi2), ("p", wp), ("df_m", float(len(slopes))), ("ll_c", rfe.ll),
               ("chi2_c", chi2_c), ("df_c", float(lv.npar())), ("p_c", _p_chibar(chi2_c, lv.npar())),
               ("n_quad", float(est.extra["npts"])), ("rank", float(ball.size)),
               ("ic", float(r.iterations)), ("converged", 1.0 if r.converged else 0.0)]
    macros = [("cmd", cmdname), ("cmdline", cmdline(cmdname, args)), ("depvar", smp.depname),
              ("ivars", lv.var), ("model", kind), ("title", _TITLE[kind]),
              ("intmethod", method), ("covariances", lv.cov), ("chi2type", "Wald"), ("vce", "oim"),
              ("predict", "meglm_p"), ("estat_cmd", "meglm_estat"), ("properties", "b V")]
    post(s, est, scalars, macros, smp.mask)
    display_me(s, est, header=not o.get("noheader"), group=not o.get("nogroup"),
               lrtest=not o.get("nolrtest"))


def display_me(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
               level: float | None = None, group: bool = True, lrtest: bool = True, eform=None, **_kw) -> None:
    out = s.output
    x = est.extra
    lev = level or est.level
    kind = x["kind"]
    if header:
        ll = f"Log likelihood = {_g(x['ll'], '%10.0g'):>10}"
        meth = x["method"]
        group_header(out, _TITLE[kind], x["groups"] if group else [], x["N"],
                     [("", "", ""),
                      (f"Integration method: {meth}", "Integration pts.", f"{int(x['npts'])}"),
                      ("", "", ""),
                      ("", f"Wald chi2({int(x['df_m'])})", _g(x["chi2"], "%10.2f")),
                      (ll, "Prob > chi2", _g(x["p"], "%10.4f"))])
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
    coef_table(s, est.depvar, rows, stat="z", df=None, level=lev,
               coef_title=ef or "Coef.", eform=bool(ef))
    if lrtest:
        lr_line(out, _VS[kind], x["chi2_c"], x["nvc"])


def _predict_me(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    from .postest import xb_all
    xb, ok = xb_all(s, est)
    ds = s.data
    G = est.glmm
    model = G["model"]
    if est.offset_name:
        nm = est.offset_name
        v = ds.get(nm[3:-1] if nm.startswith("ln(") else nm).data.astype(float)
        xb = xb + np.where(v < SYS, np.log(np.where(v > 0, v, 1)) if nm.startswith("ln(") else v, 0)
    if stat == "xb":
        return np.where(ok, xb, SYS), ""
    # posterior means of the random effects for each group of the estimation sample
    b = np.concatenate([G["beta"], G["theta"]])
    model.adapt(b, steps=50)
    rows = np.nonzero(G["mask"])[0]
    lv = G["levels"][0]
    var = ds.get(lv.var).data
    key_of = {}
    for i, g in zip(rows, lv.gid):
        key_of[var[i]] = g
    gidx = np.array([key_of.get(v, -1) for v in var])
    re_all = np.where(gidx[:, None] >= 0, model.mu[np.maximum(gidx, 0)], np.nan)

    # Z for all observations
    Zcols = []
    for nm in lv.names:
        Zcols.append(np.ones(ds.nobs) if nm == "_cons" else ds.get(nm).data.astype(float))
    Z = np.column_stack(Zcols)
    eta = xb + np.where(np.isfinite(re_all).all(axis=1), np.einsum("nq,nq->n", Z, np.nan_to_num(re_all)), 0.0)
    if stat == "fixedonly":
        eta = xb
    if stat == "eta":
        return np.where(ok, eta, SYS), ""
    kind = est.extra["kind"]
    if kind == "logit":
        mu = special.expit(eta)
    elif kind == "probit":
        mu = st.norm.cdf(eta)
    else:
        mu = np.exp(eta)
    note = "" if stat == "fixedonly" else \
        "(predictions based on fixed effects and posterior means of random effects)"
    if not stat:
        note += "\n(option mu assumed)"
    if stat != "fixedonly":
        note += f"\n(using {int(est.extra['npts'])} quadrature points)"
    return np.where(ok, mu, SYS), note


def _reffects_me(s: "Session", est: Estimates):
    """predict newvarlist, reffects: posterior means, one variable per random effect in the
    order of the level's effects (VERIFY: the order for random slopes)."""
    G = est.glmm
    lv = G["levels"][0]
    s.output.write("(calculating posterior means of random effects)\n", "text")

    def fn(mask):
        ds = s.data
        model = G["model"]
        model.adapt(np.concatenate([G["beta"], G["theta"]]), steps=50)
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
    cmd_me(s, args, "logit", "melogit")


def cmd_meprobit(s, args):
    cmd_me(s, args, "probit", "meprobit")


def cmd_mepoisson(s, args):
    cmd_me(s, args, "poisson", "mepoisson")

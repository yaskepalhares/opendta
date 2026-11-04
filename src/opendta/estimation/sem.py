"""sem ([SEM] sem): linear structural equation models by maximum likelihood.

Path models with observed variables, confirmatory factor analysis and full models with latent
variables. Paths are written as (y <- x1 x2), (x1 x2 -> y) or (F -> a b c); names that start
with a capital letter are latent. Defaults follow [SEM] intro 4:

* each latent variable is anchored by its first indicator (loading fixed at 1) and has mean 0;
* endogenous variables get an intercept and an error variance; error covariances are 0 unless
  given in cov(e.y1*e.y2);
* exogenous variables (observed and latent) have free variances and covariances, and observed
  exogenous variables free means.

The model is written in RAM form, v = A v + u, Var(u) = S, with v = (observed, latent); the
implied mean and covariance of the observed variables are F (I-A)^-1 m and
F (I-A)^-1 S (I-A)^-T F'. The likelihood is the multivariate normal one with means
(method(ml)); standard errors come from the observed information matrix.

Undocumented details are marked VERIFY: starting values, the iteration log, the labels of
long covariance names and the layout (case compat/do/1302_sem.do)."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options, parse_standard
from ..lazy import lazy
from .results import Estimates, g9, post

st = lazy("scipy.stats")
optimize = lazy("scipy.optimize")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
LINE = "-" * 78


# ---------------------------------------------------------------------------
# specification
# ---------------------------------------------------------------------------

@dataclass
class Spec:
    observed: list[str] = field(default_factory=list)
    latent: list[str] = field(default_factory=list)
    paths: list[tuple[str, str]] = field(default_factory=list)       # (from, to)
    fixed: dict[tuple[str, str], float] = field(default_factory=dict)  # path -> value (@#)
    covs: list[tuple[str, str]] = field(default_factory=list)        # extra covariances (e.x or x)
    var_fixed: dict[str, float] = field(default_factory=dict)
    labels: dict[tuple[str, str], str] = field(default_factory=dict)   # path -> symbolic label (@name)

    def all(self) -> list[str]:
        return self.observed + self.latent


def _is_latent(name: str, latent_opt: set[str], nocaps: bool) -> bool:
    if name in latent_opt:
        return True
    return not nocaps and name[:1].isupper()


def _parse_paths(s: "Session", text: str, latent_opt: set[str], nocaps: bool) -> Spec:
    spec = Spec()
    groups = re.findall(r"\(([^()]*)\)", text)
    rest = re.sub(r"\(([^()]*)\)", " ", text).strip()
    if rest:
        raise StataError(198, f"invalid path specification: {rest}")   # VERIFY
    if not groups:
        raise StataError(198, "no paths specified")   # VERIFY

    def resolve(tok: str) -> list[tuple[str, float | str | None]]:
        name, at, val = tok.partition("@")
        fixed: float | str | None = None
        if at:
            try:
                fixed = float(val)
            except ValueError:
                if not re.fullmatch(r"[A-Za-z_]\w*", val):
                    raise StataError(198, f"invalid constraint label {val}")   # VERIFY
                fixed = val
        if _is_latent(name, latent_opt, nocaps):
            if not re.fullmatch(r"[A-Za-z_]\w*", name):
                raise StataError(198, f"{name} invalid name")
            if s.data.has(name):
                raise StataError(198, f"{name} is a variable in the dataset; latent names must "
                                      "not be variable names")   # VERIFY
            if name not in spec.latent:
                spec.latent.append(name)
            return [(name, fixed)]
        names = s.expand_varlist(name)
        for n in names:
            if s.data.get(n).is_string:
                raise StataError(109, "type mismatch")
            if n not in spec.observed:
                spec.observed.append(n)
        return [(n, fixed) for n in names]

    for g in groups:
        parts = re.split(r"\s*(<-|->)\s*", g.strip())
        if len(parts) < 3:
            raise StataError(198, f"invalid path ({g})")   # VERIFY
        for k in range(1, len(parts), 2):
            arrow = parts[k]
            left = [x for tok in parts[k - 1].split() for x in resolve(tok)]
            right = [x for tok in parts[k + 1].split() for x in resolve(tok)]
            src, dst = (right, left) if arrow == "<-" else (left, right)
            for a, fa in src:
                for b, fb in dst:
                    if (a, b) not in spec.paths:
                        spec.paths.append((a, b))
                    fv = fa if fa is not None else fb
                    if isinstance(fv, str):
                        spec.labels[(a, b)] = fv
                    elif fv is not None:
                        spec.fixed[(a, b)] = fv
    return spec


# ---------------------------------------------------------------------------
# model
# ---------------------------------------------------------------------------

class Model:
    def __init__(self, spec: Spec, covs: list[tuple[str, str]], var_fix: dict[str, float],
                 latent_means: tuple[str, ...] = ()):
        self.spec = spec
        self.covs_in = list(covs)
        self.names = spec.all()
        self.p = len(spec.observed)
        self.k = len(self.names)
        idx = {n: i for i, n in enumerate(self.names)}
        self.idx = idx
        self.endog = [n for n in self.names if any(b == n for _, b in spec.paths)]
        self.exog = [n for n in self.names if n not in self.endog]
        # each latent: anchor = first observed indicator it points to (VERIFY: endogenous latent
        # variables anchored the same way)
        self.anchor: dict[str, str] = {}
        for L in spec.latent:
            fixed_one = [b for (a, b), v in spec.fixed.items() if a == L]
            if fixed_one:
                continue
            outs = [b for a, b in spec.paths if a == L and b in spec.observed]
            if not outs:
                outs = [b for a, b in spec.paths if a == L]
            if outs:
                self.anchor[L] = outs[0]
        # free parameters: list of (kind, a, b)
        self.free: list[tuple[str, str, str]] = []
        self.fixed_A: dict[tuple[int, int], float] = {}
        for a, b in spec.paths:
            if (a, b) in spec.fixed:
                self.fixed_A[(idx[b], idx[a])] = spec.fixed[(a, b)]
            elif self.anchor.get(a) == b:
                self.fixed_A[(idx[b], idx[a])] = 1.0
            else:
                self.free.append(("path", a, b))
        # means/intercepts: endogenous observed (intercept), exogenous observed (mean)
        for n in spec.observed:
            self.free.append(("mean", n, n))
        # means (exogenous) or intercepts (endogenous) of latent variables: only in the groups
        # other than the first of a multiple-group model (the first group anchors them at 0)
        for n in latent_means:
            self.free.append(("mean", n, n))
        # variances
        for n in self.names:
            if n in var_fix:
                continue
            self.free.append(("var", n, n))
        self.var_fixed = dict(var_fix)
        # covariances: exogenous among themselves; extra error covariances
        ex = [n for n in self.exog]
        for i in range(len(ex)):
            for j in range(i + 1, len(ex)):
                self.free.append(("cov", ex[i], ex[j]))
        for a, b in covs:
            if a not in idx or b not in idx:
                raise StataError(198, f"cov({a}*{b}): unknown variable")   # VERIFY
            key = ("cov",) + tuple(sorted((a, b), key=lambda n: idx[n]))
            if key not in self.free:
                self.free.append(key)
        self.nfree = len(self.free)
        self.hidden: set = set()
        # paths that share a symbolic label (@name) are constrained equal
        self.ties: dict[tuple[str, str, str], str] = {
            ("path", a, b): lab for (a, b), lab in spec.labels.items() if ("path", a, b) in set(self.free)}

    # implied moments -------------------------------------------------------
    def matrices(self, theta: np.ndarray):
        k = self.k
        A = np.zeros((k, k))
        S = np.zeros((k, k))
        m = np.zeros(k)
        for (i, j), v in self.fixed_A.items():
            A[i, j] = v
        for n, v in self.var_fixed.items():
            S[self.idx[n], self.idx[n]] = v
        for t, (kind, a, b) in zip(theta, self.free):
            if kind == "path":
                A[self.idx[b], self.idx[a]] = t
            elif kind == "mean":
                m[self.idx[a]] = t
            elif kind == "var":
                S[self.idx[a], self.idx[a]] = t
            else:
                i, j = self.idx[a], self.idx[b]
                S[i, j] = S[j, i] = t
        return A, S, m

    def implied(self, theta: np.ndarray):
        A, S, m = self.matrices(theta)
        B = np.linalg.inv(np.eye(self.k) - A)
        Sv = B @ S @ B.T
        mv = B @ m
        p = self.p
        return mv[:p], Sv[:p, :p], (A, S, m, B, Sv, mv)


def _loglik(model: Model, theta: np.ndarray, mean: np.ndarray, cov: np.ndarray, N: float) -> float:
    try:
        mu, Sig, _ = model.implied(theta)
        sign, logdet = np.linalg.slogdet(Sig)
        if sign <= 0:
            return -np.inf
        Si = np.linalg.inv(Sig)
    except np.linalg.LinAlgError:
        return -np.inf
    d = mean - mu
    p = model.p
    return float(-0.5 * N * (p * math.log(2 * math.pi) + logdet + np.trace(cov @ Si) + d @ Si @ d))


def _dmoments(model: Model, theta: np.ndarray):
    """Implied mean and covariance of the observed variables and, for each free parameter, the
    derivative of the mean (or None) and of the covariance (or None) of the observed variables,
    from the exact derivatives of the implied moments (dB = B dA B, dSv = dB S B' + B S dB' +
    B dS B', dmv = dB m + B dm)."""
    mu, Sig, (A, S, m, B, Sv, mv) = model.implied(theta)
    p = model.p
    out = []
    for kind, a, b in model.free:
        ia, ib = model.idx[a], model.idx[b]
        if kind == "path":
            dB = np.outer(B[:, ib], B[ia, :])
            dSv = dB @ S @ B.T
            dSv = dSv + dSv.T
            dmv = dB @ m
        elif kind == "mean":
            dSv = None
            dmv = B[:, ia]
        else:
            E = np.zeros_like(S)
            E[ia, ib] = E[ib, ia] = 1.0
            dSv = B @ E @ B.T
            dmv = None
        out.append((None if dmv is None else dmv[:p], None if dSv is None else dSv[:p, :p]))
    return mu, Sig, out


def _grad_blocks(model: Model, theta: np.ndarray, pats) -> np.ndarray:
    """Analytic gradient of the log likelihood of the pattern summaries `pats`: for each pattern,
    dll/dmu = n Si d and dll/dSigma = -n/2 (Si - Si (S + d d') Si) on the pattern's rows and
    columns, chained through the derivatives of the implied moments."""
    mu, Sig, dm = _dmoments(model, theta)
    p = model.p
    gmu = np.zeros(p)
    G = np.zeros((p, p))
    for q in pats:
        o = q.idx
        Si = np.linalg.inv(Sig[np.ix_(o, o)])
        d = q.mean - mu[o]
        gmu[o] += q.n * Si @ d
        G[np.ix_(o, o)] += -0.5 * q.n * (Si - Si @ (q.cov + np.outer(d, d)) @ Si)
    out = np.zeros(len(theta))
    for k, (dmv, dSv) in enumerate(dm):
        g = 0.0
        if dmv is not None:
            g += float(gmu @ dmv)
        if dSv is not None:
            g += float((G * dSv).sum())
        out[k] = g
    return out


def _loglik_blocks(model: Model, theta: np.ndarray, pats) -> float:
    from ..core import mvnmiss
    try:
        mu, Sig, _ = model.implied(theta)
        return mvnmiss.loglik(pats, mu, Sig)
    except np.linalg.LinAlgError:
        return -np.inf


def _complete_pattern(model: Model, mean: np.ndarray, cov: np.ndarray, N: float):
    from ..core.mvnmiss import Pattern
    return [Pattern(np.arange(model.p), np.zeros(0, dtype=int), N, mean, cov)]


def _grad(model: Model, theta: np.ndarray, mean: np.ndarray, cov: np.ndarray, N: float) -> np.ndarray:
    """Analytic gradient of the log likelihood from the sample mean and covariance."""
    return _grad_blocks(model, theta, _complete_pattern(model, mean, cov, N))


def _start(model: Model, mean: np.ndarray, cov: np.ndarray) -> np.ndarray:
    """Starting values (VERIFY: Stata's are undocumented): sample moments for the observed
    exogenous variables, half the observed variance for errors, latent variances from the
    anchor's variance, paths 0 (loadings 1)."""
    obs = model.spec.observed
    io = {n: i for i, n in enumerate(obs)}
    th = []
    for kind, a, b in model.free:
        if kind == "path":
            th.append(1.0 if a in model.spec.latent and b in model.spec.observed else 0.0)
        elif kind == "mean":
            th.append(float(mean[io[a]]) if a in io else 0.0)
        elif kind == "var":
            if a in io:
                v = float(cov[io[a], io[a]])
                th.append(v if a in model.exog else v / 2)
            else:
                anc = model.anchor.get(a)
                th.append(float(cov[io[anc], io[anc]]) / 2 if anc in io else 1.0)
        else:
            if a in io and b in io and a in model.exog and b in model.exog:
                th.append(float(cov[io[a], io[b]]))
            else:
                th.append(0.0)
    return np.array(th, dtype=float)


def _fit(model: Model, mean, cov, N, s: "Session", log: bool, maxiter: int = 300):
    f = lambda th: _loglik(model, th, mean, cov, N)   # noqa: E731
    gr = lambda th: _grad(model, th, mean, cov, N)   # noqa: E731

    def hess(th):
        # central differences of the analytic gradient
        k = len(th)
        H = np.zeros((k, k))
        for i in range(k):
            e = np.zeros(k)
            hh = 1e-5 * (abs(th[i]) + 1)
            e[i] = hh
            H[:, i] = (gr(th + e) - gr(th - e)) / (2 * hh)
        return (H + H.T) / 2
    th = _start(model, mean, cov)
    out = s.output
    lines = []
    ll = f(th)
    if not np.isfinite(ll):
        # VERIFY: Stata's message when the starting values are not feasible
        raise StataError(430, "initial values not feasible")
    lines.append(ll)
    # quasi-Newton to approach the maximum, Newton steps with the numerical Hessian to finish
    res = optimize.minimize(lambda t: -f(t), th, method="BFGS",
                            jac=lambda t: -gr(t), options={"maxiter": maxiter, "gtol": 1e-7})
    th = res.x
    lines.append(f(th))
    converged = False
    for _ in range(50):
        g = gr(th)
        H = hess(th)
        try:
            step = np.linalg.solve(H, -g)
        except np.linalg.LinAlgError:
            break
        lam = 1.0
        base = f(th)
        while lam > 1e-6 and not (f(th + lam * step) >= base - 1e-12):
            lam /= 2
        th = th + lam * step
        lines.append(f(th))
        if abs(lines[-1] - lines[-2]) < 1e-10 and np.max(np.abs(lam * step)) < 1e-8:
            converged = True
            break
    if log:
        out.write("\nFitting target model:\n\n", "text")
        shown = [lines[0]] + [x for i, x in enumerate(lines[1:], 1) if abs(x - lines[i - 1]) > 1e-9 or i == len(lines) - 1]
        for i, x in enumerate(shown):
            out.write(f"{'Iteration ' + str(i) + ':':<14} log likelihood = {g9(x, '%10.0g'):>10}  \n", "text")
    H = hess(th)
    try:
        V = np.linalg.inv(-H)
    except np.linalg.LinAlgError:
        V = np.linalg.pinv(-H)
    model.ilog = [lines[0]] + [x for i, x in enumerate(lines[1:], 1)
                               if abs(x - lines[i - 1]) > 1e-9 or i == len(lines) - 1]
    return th, f(th), (V + V.T) / 2, converged


# ---------------------------------------------------------------------------
# command
# ---------------------------------------------------------------------------

def _parse_options(args: str):
    c = find_top(args, ",")
    head, opts = (args, "") if c == -1 else (args[:c], args[c + 1:])
    # if/in after the paths
    m = re.search(r"\)\s*((?:if|in)\b.*)$", head, re.S)
    paths_text, tail = (head[:m.start() + 1], m.group(1)) if m else (head, "")
    p = parse_standard(tail + ("," + opts if opts.strip() else ""))
    o = match_options(p.options, {"standardized": 5, "method": 4, "nocapslatent": 9, "latent": 4,
                                  "covariance": 3, "variance": 3, "nolog": 4, "log": 3,
                                  "level": 1, "iterate": 4, "noheader": 6, "nodescribe": 6,
                                  "nofootnote": 7, "coeflegend": 5, "noxconditional": 3,
                                  "allmissing": 3, "group": 5, "ginvariant": 2})
    return paths_text, p, o


def cmd_sem(s: "Session", args: str) -> None:
    from .sem_out import _predict, display_sem
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "sem", args)
        return
    from ..commands._util import touse
    from . import semfit
    paths_text, p, o = _parse_options(args)
    meth = str(o.get("method") or "ml").strip().lower()
    if meth not in ("ml", "mlmv", "adf"):
        raise StataError(198, f"method({meth}) not allowed")   # VERIFY
    latent_opt = set(str(o["latent"]).split()) if isinstance(o.get("latent"), str) else set()
    spec = _parse_paths(s, paths_text, latent_opt, bool(o.get("nocapslatent")))
    covs = []
    for item in str(o.get("covariance") or "").split():
        a, star, b = item.partition("*")
        if not star:
            raise StataError(198, f"invalid covariance() specification: {item}")   # VERIFY
        covs.append((a.replace("e.", ""), b.replace("e.", "")))
    var_fix = {}
    for item in str(o.get("variance") or "").split():
        name, at, val = item.partition("@")
        if at:
            var_fix[name.replace("e.", "")] = float(val)
    ds = s.data
    mask = touse(s, p)
    cols = [ds.get(n).data for n in spec.observed]
    for c in cols:
        if c.dtype.kind in "OU":
            raise StataError(109, "type mismatch")
    if meth == "mlmv":
        # observations with extended missing values (.a to .z) are dropped unless allmissing
        if not o.get("allmissing"):
            for c in cols:
                mask &= ~(c > SYS)
        nonmiss = np.zeros(len(mask), dtype=bool)
        for c in cols:
            nonmiss |= c < SYS
        mask &= nonmiss
    else:
        for c in cols:
            mask &= c < SYS
    gname = None
    levels = [None]
    if o.get("group"):
        gname = str(o["group"]).strip()
        if not ds.has(gname):
            raise StataError(111, f"variable {gname} not found")
        gcol = ds.get(gname)
        if gcol.is_string:
            raise StataError(109, "type mismatch")
        mask &= gcol.data < SYS
        levels = [float(v) for v in np.unique(gcol.data[mask])]
        if len(levels) < 2:
            raise StataError(498, "group() requires at least two groups")   # VERIFY
    elif o.get("ginvariant"):
        raise StataError(198, "option ginvariant() requires option group()")   # VERIFY
    if not mask.any():
        raise StataError(2000, "no observations")
    G = len(levels)
    inv = semfit.parse_ginvariant(o.get("ginvariant") if G > 1 else None)
    # the models: observed exogenous variables are conditioned on (Stata 14's xconditional)
    # unless the moments of one of them are constrained, noxconditional is given, or method(mlmv)
    # has missing values among them; method(adf) always uses all moments (VERIFY)
    X_all = np.column_stack([c.astype(float) for c in cols])
    X_all = np.where(np.column_stack([c >= SYS for c in cols]), np.nan, X_all)
    # e(N_missing)
    n_missing = int(np.isnan(X_all[mask]).sum())      # missing values, not observations (case 1909)
    models, datas = [], []
    for g, lev in enumerate(levels):
        mg = mask if gname is None else (mask & (ds.get(gname).data == lev))
        Xg = X_all[mg]
        lat_means: tuple[str, ...] = ()
        if g > 0 and "mcons" in inv:
            base = Model(spec, covs, var_fix)
            lat_means = tuple(L for L in spec.latent
                              if not (L in base.exog and "meanex" in inv)
                              and not (L in base.endog and "scons" in inv))
        models.append(Model(spec, covs, var_fix, lat_means))
        datas.append(semfit.build_group_data(models[-1], Xg, meth, float(lev) if lev is not None else 0.0))
    model = models[0]
    N = sum(d.N for d in datas)
    pobs = model.p
    prob = semfit.Problem(models, datas, meth, inv)
    df_s = pobs * (pobs + 3) // 2
    df = G * df_s - prob.nglob
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    if not o.get("nodescribe"):
        _describe_vars(s, model)
    if df < 0:
        raise StataError(503, "the model is not identified: more parameters than "
                              "moments")   # VERIFY
    # mlmv fits the saturated and the baseline models, adf the baseline model; Stata shows their
    # logs even with nolog (case 1909; their iteration paths are VERIFY, one line is shown)
    if meth in ("mlmv", "adf"):
        out = s.output
        steps = ([("saturated", "log likelihood", sum(d.ll_s for d in datas))] if meth == "mlmv" else [])
        steps.append(("baseline", "log likelihood" if meth == "mlmv" else "discrepancy",
                      sum(d.ll_b for d in datas) if meth == "mlmv"
                      else sum(d.adf["F_b"] for d in datas)))
        for what, word, v in steps:
            out.write(f"\nFitting {what} model:\n\n", "text")
            out.write(f"Iteration 0:   {word} = {g9(v, '%10.0g'):>10}  \n", "text")
    th, obj, V, conv = semfit.fit_problem(prob, s, not o.get("nolog"))
    ll_s = sum(d.ll_s for d in datas)
    ll_b = sum(d.ll_b for d in datas)
    ex_n = sum(1 for n in spec.observed if n in model.exog) if baseline_ex_correlated(model) else 0
    df_bs = G * df_s - G * (2 * pobs + ex_n * (ex_n - 1) // 2)
    if meth == "adf":
        ll = SYS
        chi2_ms = max(-obj, 0.0)
        chi2_bs = sum(d.N * d.adf["F_b"] for d in datas)
    else:
        ll = obj
        chi2_ms = max(2 * (ll_s - ll), 0.0)
        chi2_bs = max(2 * (ll_s - ll_b), 0.0)
    p_ms = float(st.chi2.sf(chi2_ms, df)) if df > 0 else SYS
    p_bs = float(st.chi2.sf(chi2_bs, df_bs)) if df_bs > 0 else SYS
    mean, cov = datas[0].mean, datas[0].cov
    base_ll = ll_b
    if G > 1:
        from .semgroup import finish_groups
        finish_groups(s, args, prob, spec, th, obj, V, conv, o, dict(
            N=N, ll=ll, df=df, chi2_ms=chi2_ms, p_ms=p_ms, chi2_bs=chi2_bs, df_bs=df_bs, p_bs=p_bs,
            ll_s=ll_s, ll_b=ll_b, level=level, mask=mask, gname=gname, levels=levels, meth=meth))
        return
    names = _param_names(model)
    th_glob, V_glob = th, V
    if prob.nglob != model.nfree:          # paths that share a label: one global entry for several items
        Tm = prob.expand_matrix()
        th, V = Tm @ th, Tm @ V @ Tm.T
    # observed exogenous variables are conditioned on: their means, variances and covariances
    # are not parameters of e(b) and are not shown (compat 1302)
    oxset = {n for n in model.exog if n in spec.observed}
    user_cov = {v for pair in covs for v in pair} | set(var_fix)
    miss_ex = meth == "mlmv" and any(bool(np.isnan(X_all[mask][:, spec.observed.index(n)]).any())
                                     for n in oxset)
    xcond = (bool(oxset) and not (oxset & user_cov) and not o.get("noxconditional")
             and meth == "ml" and not miss_ex)
    keep = [i for i, (kind, a, b) in enumerate(model.free)
            if not (xcond and kind in ("mean", "var", "cov") and a in oxset and b in oxset)]
    model.hidden = set() if not xcond else {model.free[i] for i in range(model.nfree) if i not in keep}
    est = Estimates("sem", " ".join(model.endog), [names[i][1] for i in keep], th[keep].copy(),
                    V[np.ix_(keep, keep)], N, "z", None)
    est.eqnames = [names[i][0] for i in keep]
    est.model = model
    est.level = level
    est.extra = dict(N=N, ll=ll, df=df, chi2_ms=chi2_ms, p_ms=p_ms, chi2_bs=chi2_bs, df_bs=df_bs,
                     p_bs=p_bs, mean=mean, cov=cov, ll_s=ll_s, standardized=bool(o.get("standardized")),
                     th=th.copy(), V=V, xcond=xcond, method=meth, prob=prob, mask=mask,
                     th_glob=th_glob, V_glob=V_glob)
    est.display = display_sem
    est.predict = _predict
    est.predict_opts = {"latent": 3, "xb": 2}
    from .semestat import estat_table
    est.estat = estat_table()
    lat_endog = [n for n in model.endog if n in spec.latent]
    obs_endog = [n for n in model.endog if n in spec.observed]
    lat_exog = [n for n in model.exog if n in spec.latent]
    obs_exog = [n for n in model.exog if n in spec.observed]
    # e() in Stata 14's ereturn list order (compat 1302; VERIFY the matrices' contents)
    ilog = list(prob.ilog) if prob.ilog else [ll]
    scalars = [("converged", 1.0 if conv else 0.0), ("rc", 0.0), ("ic", float(len(ilog) - 1)),
               ("N_groups", 1.0), ("modelmeans", 1.0), ("N", float(N)), ("N_missing", float(n_missing)),
               ("ll", ll), ("critvalue", ll if meth != "adf" else -obj / max(N, 1)),
               ("df_m", float(len(keep))),
               ("rank", float(np.linalg.matrix_rank(est.V))), ("estimates", 1.0),
               ("critvalue_s", ll_s), ("critvalue_b", base_ll),
               ("p_bs", p_bs), ("df_bs", float(df_bs)), ("chi2_bs", chi2_bs), ("p_ms", p_ms),
               ("df_ms", float(df)), ("chi2_ms", chi2_ms), ("df_b", float(df_s - df_bs)), ("df_s", float(df_s))]
    if meth != "mlmv":
        # VERIFY: df_s follows critvalue_b with method(ml) (compat 1302) and comes last with
        # method(mlmv) (compat 1909)
        scalars.insert(14, scalars.pop())
    if meth == "adf":
        scalars = [(k, v) for k, v in scalars if k not in ("ll", "critvalue_s", "critvalue_b")]
    first = obs_endog[0] if obs_endog else ""
    # VERIFY: _NONE when there is no observed exogenous variable (case 1909)
    covariates = " ".join(obs_exog) if xcond else ("_NONE" if not obs_exog else "")
    macros = [("cmdline", f"sem {args.strip()}"), ("covariates", covariates),
              ("cmd", "sem"), ("predict", "sem_p"), ("estat_cmd", "sem_estat"),
              ("marginsdefault", " ".join(f"predict(xb({n}))" for n in obs_endog)),
              ("marginsnotok", "XBLATent1(passthru) XBLATent LATent1(passthru) LATent SCores xb"),
              ("xconditional", "xconditional" if xcond else ""),
              ("vcetype", "OIM"), ("vce", "oim"),
              ("title", "Structural equation model"), ("method", meth), ("data", "raw"),
              ("lxvars", " ".join(lat_exog)), ("oxvars", " ".join(obs_exog)),
              ("lyvars", " ".join(lat_endog)), ("oyvars", " ".join(obs_endog)),
              ("technique", "nr"), ("properties", "b V")]
    from ..commands.matrix import Matrix
    cn = list(est.names)

    def row(v, rn="r1"):
        mt = Matrix(np.asarray(v, float).reshape(1, -1), [rn], cn)
        mt.coleq = list(est.eqnames)
        return mt
    g = (prob.grad(th_glob)[prob.gmap[0]][keep] if meth != "adf" else np.zeros(len(keep)))
    bstd = _standardized(model, th)[keep]
    J = np.zeros((len(keep), len(keep)))
    for jj, i in enumerate(keep):
        e_ = np.zeros(len(th))
        h = 1e-6 * (abs(th[i]) + 1)
        e_[i] = h
        J[:, jj] = (_standardized(model, th + e_)[keep] - _standardized(model, th - e_)[keep]) / (2 * h)
    Vstd = Matrix(J @ est.V @ J.T, cn, cn)
    Vstd.roweq = Vstd.coleq = list(est.eqnames)
    pclass = [{"path": 1, "mean": 2, "var": 3, "cov": 4}[model.free[i][0]] for i in keep]
    mats = [("gradient", row(g)), ("V_std", Vstd), ("b_std", row(bstd)),
            ("admissible", Matrix(np.ones((1, 3)), ["r1"], ["c1", "c2", "c3"])),
            ("b_pclass", row(pclass)), ("nobs", Matrix(np.array([[float(N)]]), ["r1"], ["c1"])),
            ("ilog", Matrix(np.array([ilog + [0.0] * (20 - len(ilog))])[:, :20], ["r1"],
                            [f"c{i}" for i in range(1, 21)]))]
    post(s, est, scalars, macros, mask, extra_mats=mats)
    display_sem(s, est)


def baseline_ex_correlated(model: Model) -> bool:
    """The baseline model leaves the observed exogenous variables correlated only when the model
    has observed endogenous variables; without any, every variable is uncorrelated ([SEM] methods
    and formulas for sem)."""
    return any(n in model.endog for n in model.spec.observed)


def _baseline(model: Model, mean, cov, N):
    """Baseline model: observed endogenous variables uncorrelated with everything; the
    observed exogenous variables keep their covariances (VERIFY)."""
    obs = model.spec.observed
    ex = [i for i, n in enumerate(obs) if n in model.exog] if baseline_ex_correlated(model) else []
    p = len(obs)
    Sig = np.diag(np.diag(cov)).copy()
    for i in ex:
        for j in ex:
            Sig[i, j] = cov[i, j]
    ll = float(-0.5 * N * (p * math.log(2 * math.pi) + np.linalg.slogdet(Sig)[1]
                           + np.trace(cov @ np.linalg.inv(Sig))))
    kfree = p + p + len(ex) * (len(ex) - 1) // 2
    df = p * (p + 3) // 2 - kfree
    return ll, df


def _param_names(model: Model) -> list[tuple[str, str]]:
    out = []
    for kind, a, b in model.free:
        if kind == "path":
            out.append((b, a))
        elif kind == "mean":
            out.append((a, "_cons") if a in model.endog else (f"mean({a})", "_cons"))
        elif kind == "var":
            out.append((f"var(e.{a})" if a in model.endog else f"var({a})", "_cons"))
        else:
            ea = f"e.{a}" if a in model.endog else a
            eb = f"e.{b}" if b in model.endog else b
            out.append((f"cov({ea},{eb})", "_cons"))
    return out


def _describe_vars(s: "Session", model: Model) -> None:
    out = s.output
    spec = model.spec
    obs_endog = [n for n in model.endog if n in spec.observed]
    meas = [n for n in obs_endog if all(a in spec.latent for a, b in spec.paths if b == n)]
    struct_obs = [n for n in obs_endog if n not in meas]
    lat_endog = [n for n in model.endog if n in spec.latent]
    out.write("\nEndogenous variables\n\n", "text")
    if struct_obs:
        out.write("Observed:  " if not meas else "Observed:     ", "text")
        out.write(" ".join(struct_obs) + "\n", "result")
    if meas:
        out.write("Measurement:  ", "text")
        out.write(" ".join(meas) + "\n", "result")
    if lat_endog:
        out.write("Latent:       ", "text")
        out.write(" ".join(lat_endog) + "\n", "result")
    oex = [n for n in model.exog if n in spec.observed]
    lex = [n for n in model.exog if n in spec.latent]
    if oex or lex:
        out.write("\nExogenous variables\n\n", "text")
        if oex:
            out.write("Observed:  " if not lex else "Observed:     ", "text")
            out.write(" ".join(oex) + "\n", "result")
        if lex:
            out.write("Latent:       ", "text")
            out.write(" ".join(lex) + "\n", "result")


# ---------------------------------------------------------------------------
# display
# ---------------------------------------------------------------------------

def _standardized(model: Model, th: np.ndarray) -> np.ndarray:
    """Standardized values of the parameters: paths b * sd(from) / sd(to), variances over the
    total variance, covariances as correlations; means and intercepts unchanged (VERIFY)."""
    *_, (A, S, m, B, Sv, mv) = model.implied(th)
    sd = np.sqrt(np.maximum(np.diag(Sv), 1e-300))
    out = th.copy()
    for k, (kind, a, b) in enumerate(model.free):
        ia, ib = model.idx[a], model.idx[b]
        if kind == "path":
            out[k] = th[k] * sd[ia] / sd[ib]
        elif kind == "var":
            out[k] = th[k] / Sv[ia, ia]
        elif kind == "cov":
            out[k] = th[k] / math.sqrt(S[ia, ia] * S[ib, ib]) if S[ia, ia] > 0 and S[ib, ib] > 0 else 0.0
        else:
            out[k] = th[k] / sd[ia] if a in model.exog else th[k] / sd[ia]
    return out

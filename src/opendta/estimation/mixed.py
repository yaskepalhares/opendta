"""Linear mixed models: mixed ([ME] mixed, [ME] mixed postestimation).

y = Xb + Zu + e, with nested levels of random effects (|| school: || class:),
identity, independent, exchangeable or unstructured covariances, by maximum
likelihood (default) or REML. The likelihood is computed per top-level
group with the Woodbury identity, so the cost depends on the number of
random effects and not on the group sizes; b is profiled out
([ME] mixed, Methods and formulas). The variance components are optimized
in the ln(sd) and atanh(corr) metrics and reported with delta-method
standard errors; intervals use the log scale for variances and the atanh
scale for correlations.

Stata starts with EM iterations that it does not log and then shows the
gradient-based ones; OpenDTA does a silent rough maximization and then
logs the Newton iterations, so the number of iterations may differ
(VERIFY: compat/do/1001_mixed.do).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..core.stats import chi2_tail
from ..lang.syntax import find_top
from ..lazy import lazy
from . import fvars
from .ml import maximize, print_log
from .ml import options as ml_options
from .results import Estimates, coef_table, post, rows_from_columns

# imported on first use (start-up time; see opendta/lazy.py)
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _g(x: float, fmt: str = "%9.0g") -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


# ---------------------------------------------------------------------------
# syntax: fixed || level: [varlist] [, re options] || ... [, options]
# ---------------------------------------------------------------------------

_RE_OPTS = {"covariance": 3, "noconstant": 4, "collinear": 4, "fweight": 2, "pweight": 2}


def split_levels(args: str) -> tuple[str, list[tuple[str, str, dict]], str]:
    """Returns (fixed part, [(level var, re varlist, re options)], global options)."""
    parts = []
    depth, start, i = 0, 0, 0
    quote = False
    while i < len(args):
        c = args[i]
        if c == '"':
            quote = not quote
        elif not quote and c in "([":
            depth += 1
        elif not quote and c in ")]":
            depth -= 1
        elif not quote and depth == 0 and args.startswith("||", i):
            parts.append(args[start:i])
            start = i + 2
            i += 2
            continue
        i += 1
    parts.append(args[start:])
    fixed = parts[0]
    levels = []
    glob = ""
    for k, part in enumerate(parts[1:]):
        c = find_top(part, ",")
        body, opts = (part, "") if c == -1 else (part[:c], part[c + 1:])
        if ":" not in body:
            raise StataError(198, "invalid random-effects specification")   # VERIFY
        lv, revars = body.split(":", 1)
        lv = lv.strip()
        reo: dict = {}
        if opts.strip():
            from ..lang.syntax import parse_options
            rest = []
            for name, arg in parse_options(opts):
                full = next((f for f, m in _RE_OPTS.items() if len(name) >= m and f.startswith(name)), None)
                if full is not None:
                    reo[full] = arg if arg is not None else True
                else:
                    rest.append(name + (f"({arg})" if arg is not None else ""))
            if rest:
                if k != len(parts) - 2:
                    raise StataError(198, f"option {rest[0]} not allowed")
                glob = " ".join(rest)
        levels.append((lv, revars.strip(), reo))
    # global options may also follow the fixed part when there are no levels
    return fixed, levels, glob


def _cov_kind(text: str | bool | None, q: int) -> str:
    if not isinstance(text, str):
        return "identity" if q == 1 else "independent"
    t = text.strip().lower()
    for full, k in (("independent", 3), ("exchangeable", 2), ("identity", 2), ("unstructured", 2)):
        if len(t) >= k and full.startswith(t):
            return full
    raise StataError(198, f"covariance({text}) invalid")   # VERIFY


@dataclass
class Level:
    var: str
    names: list[str]          # random-effect names (variables, then _cons)
    cov: str
    Z: np.ndarray             # n x q over the estimation sample
    gid: np.ndarray           # group index over the sample (nested in the upper levels)

    @property
    def q(self) -> int:
        return len(self.names)

    def npar(self) -> int:
        q = self.q
        return {"identity": 1, "independent": q, "exchangeable": 2,
                "unstructured": q + q * (q - 1) // 2}[self.cov]

    def sigma(self, th: np.ndarray) -> np.ndarray:
        q = self.q
        if self.cov == "identity":
            return np.exp(2 * th[0]) * np.eye(q)
        if self.cov == "independent":
            return np.diag(np.exp(2 * th[:q]))
        if self.cov == "exchangeable":
            s2 = np.exp(2 * th[0])
            r = np.tanh(th[1])
            S = np.full((q, q), r * s2)
            np.fill_diagonal(S, s2)
            return S
        sd = np.exp(th[:q])
        R = np.eye(q)
        k = q
        for i in range(q):
            for j in range(i + 1, q):
                R[i, j] = R[j, i] = np.tanh(th[k])
                k += 1
        return R * np.outer(sd, sd)

    def labels(self) -> list[str]:
        """Names of the parameters in the estimation metric (for e(b))."""
        q = self.q
        if self.cov == "identity":
            return ["lns_1"]
        if self.cov == "independent":
            return [f"lns_{i + 1}" for i in range(q)]
        if self.cov == "exchangeable":
            return ["lns_1", "atr_1_1"]
        return [f"lns_{i + 1}" for i in range(q)] + [f"atr_{i + 1}_{j + 1}" for i in range(q) for j in range(i + 1, q)]


class _Cluster:
    """Precomputed cross products of a top-level group."""

    def __init__(self, rows, X, y, levels: list[Level]):
        self.rows = rows
        Zs = []
        self.blocks = []        # (level index, number of copies)
        for li, lv in enumerate(levels):
            g = lv.gid[rows]
            ug = np.unique(g)
            for u in ug:
                Zb = np.zeros((len(rows), lv.q))
                m = g == u
                Zb[m] = lv.Z[rows][m]
                Zs.append(Zb)
            self.blocks.append((li, len(ug)))
        self.Z = np.column_stack(Zs) if Zs else np.zeros((len(rows), 0))
        Xr, yr = X[rows], y[rows]
        self.n = len(rows)
        self.A = self.Z.T @ self.Z
        self.ZX = self.Z.T @ Xr
        self.Zy = self.Z.T @ yr
        self.XX = Xr.T @ Xr
        self.Xy = Xr.T @ yr
        self.yy = float(yr @ yr)


def _G(cl: _Cluster, sigmas: list[np.ndarray]) -> np.ndarray:
    Q = cl.Z.shape[1]
    G = np.zeros((Q, Q))
    k = 0
    for li, copies in cl.blocks:
        S = sigmas[li]
        q = S.shape[0]
        for _ in range(copies):
            G[k:k + q, k:k + q] = S
            k += q
    return G


class LMM:
    def __init__(self, X, y, levels: list[Level], reml: bool):
        self.X, self.y, self.levels, self.reml = X, y, levels, reml
        self.n, self.p = X.shape
        top = levels[0].gid
        order = np.argsort(top, kind="stable")
        bounds = np.nonzero(np.diff(top[order]))[0] + 1
        self.clusters = [_Cluster(rows, X, y, levels) for rows in np.split(order, bounds)]
        self.sizes = [lv.npar() for lv in levels]

    def split(self, th):
        out, k = [], 0
        for s in self.sizes:
            out.append(th[k:k + s])
            k += s
        return out, th[k]

    def parts(self, th):
        """Aggregates for the profiled likelihood at theta."""
        pieces, lnse = self.split(th)
        sigmas = [lv.sigma(t) for lv, t in zip(self.levels, pieces)]
        s2 = math.exp(2 * lnse)
        p = self.p
        XVX = np.zeros((p, p))
        XVy = np.zeros(p)
        yVy = 0.0
        logdet = 0.0
        cache = []
        for cl in self.clusters:
            G = _G(cl, sigmas)
            Q = G.shape[0]
            K = s2 * np.eye(Q) + G @ cl.A
            sign, ld = np.linalg.slogdet(K)
            if sign <= 0:
                return None
            # V^-1 = (I - Z K^-1 G Z') / s2
            KG = np.linalg.solve(K, G)
            XVX += (cl.XX - cl.ZX.T @ KG @ cl.ZX) / s2
            XVy += (cl.Xy - cl.ZX.T @ KG @ cl.Zy) / s2
            yVy += (cl.yy - cl.Zy @ KG @ cl.Zy) / s2
            logdet += cl.n * math.log(s2) + ld - Q * math.log(s2)
            cache.append((G, KG))
        return XVX, XVy, yVy, logdet, cache, sigmas, s2

    def loglik(self, th) -> float:
        r = self.parts(th)
        if r is None:
            return -np.inf
        XVX, XVy, yVy, logdet, *_ = r
        try:
            b = np.linalg.solve(XVX, XVy)
        except np.linalg.LinAlgError:
            return -np.inf
        q = yVy - 2 * b @ XVy + b @ XVX @ b
        if self.reml:
            s, ldx = np.linalg.slogdet(XVX)
            return -0.5 * ((self.n - self.p) * math.log(2 * math.pi) + logdet + ldx + q)
        return -0.5 * (self.n * math.log(2 * math.pi) + logdet + q)

    def beta(self, th):
        XVX, XVy, *_ = self.parts(th)
        b = np.linalg.solve(XVX, XVy)
        return b, np.linalg.inv(XVX)

    def blups(self, th, b) -> list[np.ndarray]:
        """Empirical Bayes predictions of the random effects of each cluster (stacked
        as in Z) : u = G Z' V^-1 r."""
        _, _, _, _, cache, sigmas, s2 = self.parts(th)
        out = []
        for cl, (G, KG) in zip(self.clusters, cache):
            Zr = cl.Zy - cl.ZX @ b
            out.append((G @ Zr - G @ cl.A @ KG @ Zr) / s2)
        return out


def _numderiv(f, x):
    """Gradient and Hessian by central differences (steps scaled to x)."""
    k = x.size
    h = 1e-4 * (np.abs(x) + 1)
    f0 = f(x)
    g = np.zeros(k)
    H = np.zeros((k, k))
    fp = np.zeros(k)
    fm = np.zeros(k)
    for i in range(k):
        e = np.zeros(k)
        e[i] = h[i]
        fp[i], fm[i] = f(x + e), f(x - e)
        g[i] = (fp[i] - fm[i]) / (2 * h[i])
        H[i, i] = (fp[i] - 2 * f0 + fm[i]) / h[i] ** 2
    for i in range(k):
        for j in range(i + 1, k):
            ei = np.zeros(k)
            ej = np.zeros(k)
            ei[i], ej[j] = h[i], h[j]
            v = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4 * h[i] * h[j])
            H[i, j] = H[j, i] = v
    return f0, g, H


# ---------------------------------------------------------------------------
# the command
# ---------------------------------------------------------------------------

_GLOBAL = {"mle": 3, "reml": 4, "variance": 3, "stddeviations": 3, "nolog": 4, "log": 3,
           "iterate": 4, "tolerance": 3, "ltolerance": 4, "nrtolerance": 4, "level": 1,
           "noheader": 6, "nogroup": 5, "nolrtest": 5, "noretable": 5, "nofetable": 5,
           "noconstant": 4, "emiterate": 3, "emonly": 3, "emlog": 3, "estmetric": 4, "matlog": 4,
           "nostderr": 6, "vce": 3}


def build_levels(s: "Session", mask: np.ndarray, levels_spec, n: int) -> list[Level]:
    from ..core.varlist import resolve_name
    ds = s.data
    out = []
    prev_key = np.zeros(int(mask.sum()), dtype=np.int64)
    for lv, revars, reo in levels_spec:
        if lv == "_all":
            raise StataError(198, "crossed random effects (_all:) are not available in OpenDTA")   # OpenDTA
        var = resolve_name(ds, lv)
        names = s.expand_varlist(revars) if revars else []
        cols = [ds.get(v).data[mask].astype(np.float64) for v in names]
        cons = not reo.get("noconstant")
        if cons:
            cols.append(np.ones(int(mask.sum())))
        if not cols:
            raise StataError(198, "random-effects equation has no variables")   # VERIFY
        Z = np.column_stack(cols)
        raw = ds.get(var).data[mask]
        if ds.get(var).is_string:
            _, codes = np.unique(raw, return_inverse=True)
        else:
            _, codes = np.unique(raw.astype(np.float64), return_inverse=True)
        # nested in the levels above: the group is the combination with the upper key
        key = prev_key * (codes.max() + 1) + codes
        _, gid = np.unique(key, return_inverse=True)
        prev_key = gid.astype(np.int64)
        q = len(names) + (1 if cons else 0)
        out.append(Level(var, names + (["_cons"] if cons else []), _cov_kind(reo.get("covariance"), q), Z, gid))
    return out


def cmd_mixed(s: "Session", args: str, cmdname: str = "mixed") -> None:
    from .postest import replay
    from .regress import Sample, cmdline
    if not args.strip() or args.strip().startswith(","):
        replay(s, "mixed", args)
        return
    fixed, levels_spec, glob = split_levels(args)
    if not levels_spec:
        raise StataError(198, "random-effects equation required")   # VERIFY
    c = find_top(fixed, ",")
    fixed_body, fe_opts = (fixed, "") if c == -1 else (fixed[:c], fixed[c + 1:])
    opts_text = " ".join(x for x in (fe_opts.strip(), glob.strip()) if x)
    extra = " ".join([lv for lv, _r, _o in levels_spec if lv != "_all"] +
                     [r for _l, r, _o in levels_spec if r])
    smp = Sample(s, fixed_body + (f", {opts_text}" if opts_text else ""), opts_spec=_GLOBAL,
                 weights=(), extra_vars=extra)
    o = smp.o
    reml = bool(o.get("reml"))
    out = s.output
    ds = s.data
    fvars.drop_collinear(smp.cols, smp.constant)
    est_cols = [col for col in smp.cols if col.values is not None]
    n = smp.n
    X = np.column_stack([col.values for col in est_cols] + ([np.ones(n)] if smp.constant else [])) \
        if (est_cols or smp.constant) else np.zeros((n, 0))
    y = smp.y.astype(np.float64)
    levels = build_levels(s, smp.mask, levels_spec, n)
    model = LMM(X, y, levels, reml)
    # starting values from OLS
    b_ols, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ b_ols
    s2 = float(res @ res / max(n - X.shape[1], 1))
    th0 = []
    for lv in levels:
        sd0 = []
        for j, nm in enumerate(lv.names):
            zsd = float(np.std(lv.Z[:, j])) if nm != "_cons" else 1.0
            sd0.append(0.5 * math.log(s2 / 4) - (math.log(zsd) if zsd > 0 and nm != "_cons" else 0.0))
        if lv.cov == "identity":
            th0.append(sd0[0])
        elif lv.cov == "independent":
            th0 += sd0
        elif lv.cov == "exchangeable":
            th0 += [sd0[0], 0.0]
        else:
            th0 += sd0 + [0.0] * (lv.q * (lv.q - 1) // 2)
    th0.append(0.5 * math.log(s2 / 2))
    th0 = np.array(th0)

    def fun(th, todo):
        if todo == 0:
            return model.loglik(th), None, None
        return _numderiv(model.loglik, th)

    # silent rough maximization (in place of Stata's EM iterations)
    pre = maximize(fun, th0, maxiter=50, tol=1e-2, ltol=1e-3, nrtol=1e2)
    r = maximize(fun, pre.b, **ml_options(o))
    lab = "log restricted-likelihood" if reml else "log likelihood"
    if not o.get("nolog"):
        out.write("\nPerforming EM optimization: \n\nPerforming gradient-based optimization: \n\n", "text")
        print_log(out, r, lab)
        out.write("\nComputing standard errors:\n", "text")
    th = r.b
    b, Vb = model.beta(th)
    try:
        Vth = np.linalg.inv(-r.H)
        Vth = (Vth + Vth.T) / 2
    except np.linalg.LinAlgError:
        Vth = np.full((th.size, th.size), np.nan)
    # linear model for the LR test (same criterion: ML or REML)
    rss = float(res @ res)
    p = X.shape[1]
    if reml:
        s2r = rss / (n - p)
        _sgn, ldx = np.linalg.slogdet(X.T @ X / s2r)
        ll_lin = -0.5 * ((n - p) * math.log(2 * math.pi) + n * math.log(s2r) + ldx + (n - p))
    else:
        ll_lin = -0.5 * n * (math.log(2 * math.pi) + math.log(rss / n) + 1)
    chi2_c = max(2 * (r.ll - ll_lin), 0.0)
    nvc = sum(lv.npar() for lv in levels)
    # Wald test of the slopes
    slopes = list(range(len(est_cols)))
    try:
        wchi2 = float(b[slopes] @ np.linalg.solve(Vb[np.ix_(slopes, slopes)], b[slopes])) if slopes else SYS
    except np.linalg.LinAlgError:
        wchi2 = SYS
    wp = float(chi2_tail(wchi2, len(slopes))) if slopes and wchi2 < SYS else SYS
    # e(b): fixed part with the full names, then the variance parameters
    names = [col.name if not col.omitted else fvars.omitted_name(col.name) for col in smp.cols]
    if smp.constant:
        names.append("_cons")
    kk = len(names)
    pos = [i for i, col in enumerate(smp.cols) if col.values is not None] + ([kk - 1] if smp.constant else [])
    bf = np.zeros(kk)
    Vf = np.zeros((kk, kk))
    bf[pos] = b
    Vf[np.ix_(pos, pos)] = Vb
    th_names, th_eqs = [], []
    for li, lv in enumerate(levels):
        for lab_ in lv.labels():
            th_eqs.append(f"{lab_.split('_')[0]}{li + 1}_{lab_.split('_', 1)[1]}")
            th_names.append("_cons")
    th_eqs.append("lnsig_e")
    th_names.append("_cons")
    ball = np.concatenate([bf, th])
    Vall = np.zeros((kk + th.size, kk + th.size))
    Vall[:kk, :kk] = Vf
    Vall[kk:, kk:] = Vth
    bmap = {col.name: bf[i] for i, col in enumerate(smp.cols)}
    semap = {col.name: float(np.sqrt(max(Vf[i, i], 0))) for i, col in enumerate(smp.cols)}
    if smp.constant:
        bmap["_cons"], semap["_cons"] = bf[-1], float(np.sqrt(max(Vf[-1, -1], 0)))
    est = Estimates("mixed", smp.depname, names + th_names, ball, Vall, n, "z", None)
    est.eqnames = [smp.depname] * kk + th_eqs
    est.rows = rows_from_columns(smp.cols, bmap, semap, ds, constant=smp.constant)
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = ""
    est.robust = False
    est.lmm = {"model": model, "theta": th, "Vth": Vth, "beta": b, "levels": levels, "reml": reml,
               "kfixed": kk, "mask": smp.mask.copy()}
    est.extra = {"N": float(n), "ll": r.ll, "chi2": wchi2, "p": wp, "df_m": len(slopes), "ll_c": ll_lin,
                 "chi2_c": chi2_c, "nvc": nvc, "reml": reml, "sd": bool(o.get("stddeviations"))}
    est.display = display_mixed
    est.predict = _predict_mixed
    est.predict_opts = {"fitted": 3, "reffects": 4, "residuals": 3, "rstandard": 4}
    est.predict_multi = {"reffects": _reffects}
    est.estat = {"icc": (3, estat_icc), "group": (3, estat_group), "recovariance": (5, estat_recov)}
    est.ml_cmd = "mixed"
    # groups per level
    ginfo = []
    for lv in levels:
        cnt = np.bincount(lv.gid)
        ginfo.append((lv.var, len(cnt), int(cnt.min()), float(cnt.mean()), int(cnt.max())))
    est.extra["groups"] = ginfo
    # e() in the order Stata lists it (compat 1001)
    nvar_p = sum({"identity": 1, "exchangeable": 1}.get(lv.cov, lv.q) for lv in levels)
    ncov_p = nvc - nvar_p
    scalars = [("rank", float(ball.size)), ("k_f", float(kk)), ("k_r", float(th.size)),
               ("k", float(ball.size)), ("k_rs", float(nvar_p + 1)), ("k_rc", float(ncov_p)),
               ("k_res", 0.0), ("converged", 1.0 if r.converged else 0.0), ("ll", r.ll),
               ("ll_c", ll_lin), ("df_c", float(nvc)), ("chi2_c", chi2_c),
               ("p_c", _p_chibar(chi2_c, nvc)), ("N", float(n)), ("nrgroups", float(len(levels))),
               ("rc", 0.0), ("df_m", float(len(slopes))), ("chi2", wchi2), ("p", wp), ("small", 0.0)]
    from .results import datasignature
    sigvars = [smp.depname] + fvars.base_vars(smp.terms)
    for lv in levels:
        for nm in lv.names + [lv.var]:
            if nm != "_cons" and nm not in sigvars and ds.has(nm):
                sigvars.append(nm)
    macros = [("cmdline", cmdline(cmdname, args)), ("datasignaturevars", " ".join(sigvars)),
              ("datasignature", datasignature(s, sigvars, smp.mask)), ("chi2type", "Wald"),
              ("rstructlab", "Independent"), ("rstructure", "independent"),
              ("estat_cmd", "mixed_estat"), ("predict", "mixed_p"),
              ("redim", " ".join(str(lv.q) for lv in levels)),
              ("ivars", " ".join(lv.var for lv in levels)),
              ("vartypes", " ".join({"identity": "Identity", "independent": "Independent",
                                    "exchangeable": "Exchangeable", "unstructured": "Unstructured"}[lv.cov]
                                   for lv in levels)),
              ("revars", " ".join(" ".join(lv.names) for lv in levels)), ("cmd", "mixed"),
              ("title", "Mixed-effects REML regression" if reml else "Mixed-effects ML regression"),
              ("technique", "nr"), ("ml_method", "d0"), ("opt", "moptimize"),
              ("method", "REML" if reml else "ML"), ("optmetric", "matsqrt"), ("depvar", smp.depname),
              ("properties", "b V")]
    from ..commands.matrix import Matrix
    gcols = [lv.var for lv in levels]
    extra_mats = [(nm, Matrix(np.array([[fn(g) for g in ginfo]], dtype=float), ["r1"], gcols))
                  for nm, fn in (("g_max", lambda g: g[4]), ("g_avg", lambda g: g[3]),
                                 ("g_min", lambda g: g[2]), ("N_g", lambda g: g[1]))]
    post(s, est, scalars, macros, smp.mask, extra_mats=extra_mats)
    display_mixed(s, est, header=not o.get("noheader"), table=not o.get("nofetable"),
                  retable=not o.get("noretable"), group=not o.get("nogroup"),
                  lrtest=not o.get("nolrtest"))


def _p_chibar(chi2: float, k: int) -> float:
    if k == 1:
        return 0.5 * float(st.chi2.sf(chi2, 1)) if chi2 > 0 else 1.0
    return float(chi2_tail(chi2, k))


# ---------------------------------------------------------------------------
# display
# ---------------------------------------------------------------------------

def _right(out, left: str, lab: str, val: str) -> None:
    out.write(f"{left:<48}{lab:<18}= ", "text")
    out.write(f"{val:>10}\n", "result")


def group_header(out, title: str, ginfo, N: float, extra_right=(), *, gv_right: bool = False,
                 mid=()) -> None:
    """gv_right: the me commands right-align the group variable ("Group variable:   g",
    compat 1002); mixed writes it after one space."""
    out.write("\n", "text")
    _right(out, title, "Number of obs", f"{int(N):,}")
    for line in mid:                  # Family:/Link: of meglm, Overdispersion: of menbreg
        out.write(line + "\n", "text")
    if len(ginfo) == 1:
        var, ng, mn, av, mx = ginfo[0]
        gv = f"Group variable:{var:>16}" if gv_right else f"Group variable: {var}"
        _right(out, gv, "Number of groups", f"{ng:,}")
        out.write("\n", "text")
        out.write(f"{'':48}Obs per group:\n", "text")
        for lab, val in (("min", f"{mn:,}"), ("avg", f"{av:.1f}"), ("max", f"{mx:,}")):
            out.write(f"{'':48}{lab:>17} = ", "text")
            out.write(f"{val:>10}\n", "result")
    else:
        out.write("\n" + "-" * 61 + "\n", "text")
        out.write("                |     No. of       Observations per Group\n", "text")
        out.write(" Group Variable |     Groups    Minimum    Average    Maximum\n", "text")
        out.write("----------------+--------------------------------------------\n", "text")
        for var, ng, mn, av, mx in ginfo:
            out.write(f"{var[:14]:>15} |", "text")
            out.write(f" {ng:>10,} {mn:>10,} {av:>10.1f} {mx:>10,}\n", "result")
        out.write("-" * 61 + "\n", "text")
    for left, lab, val in extra_right:
        if lab:
            _right(out, left, lab, val)
        else:
            out.write(left + "\n", "text")


def display_mixed(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                  level: float | None = None, retable: bool = True, group: bool = True,
                  lrtest: bool = True, **_kw) -> None:
    out = s.output
    x = est.extra
    lev = level or est.level
    if header:
        title = "Mixed-effects REML regression" if x["reml"] else "Mixed-effects ML regression"
        lab = "Log restricted-likelihood" if x["reml"] else "Log likelihood"
        ll = f"{lab} = {_g(x['ll'], '%10.0g'):>10}"
        group_header(out, title, x["groups"] if group else [], x["N"],
                     [("", "", ""), ("", f"Wald chi2({int(x['df_m'])})", _g(x["chi2"], "%10.2f")),
                      (ll, "Prob > chi2", _g(x["p"], "%10.4f"))])
        out.write("\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=lev)
    if retable:
        re_table(s, est, lev)
    if lrtest:
        lr_line(out, "linear model", x["chi2_c"], x["nvc"])


def lr_line(out, against: str, chi2: float, k: int) -> None:
    if k == 1:
        left = f"LR test vs. {against}: chibar2(01) = {chi2:.2f}"
        right = f"Prob >= chibar2 = {_p_chibar(chi2, 1):.4f}"
    else:
        left = f"LR test vs. {against}: chi2({k}) = {chi2:.2f}"
        right = f"Prob > chi2 = {chi2_tail(chi2, k):.4f}"
    out.write(left + right.rjust(max(78 - len(left), len(right) + 1)) + "\n", "text")
    if k > 1:
        out.write("\nNote: LR test is conservative and provided only for reference.\n", "text")


def _vc_rows(lv: Level, th: np.ndarray, V: np.ndarray, sd: bool, z: float):
    """(label, estimate, se, lo, hi) for the variance components of a level, by the delta method."""
    rows = []
    q = lv.q
    S = lv.sigma(th)

    def num_jac(fn):
        h = 1e-6
        f0 = fn(th)
        J = np.zeros(th.size)
        for i in range(th.size):
            e = np.zeros(th.size)
            e[i] = h
            J[i] = (fn(th + e) - fn(th - e)) / (2 * h)
        return f0, J

    def add(label, fn, kind):
        f0, J = num_jac(fn)
        se = float(np.sqrt(max(J @ V @ J, 0)))
        if kind == "pos":           # variances and sds: interval on the log scale
            lse = se / f0 if f0 > 0 else np.nan
            lo, hi = f0 * math.exp(-z * lse), f0 * math.exp(z * lse)
        elif kind == "corr":
            a = math.atanh(f0) if abs(f0) < 1 else np.nan
            ase = se / (1 - f0 * f0) if abs(f0) < 1 else np.nan
            lo, hi = math.tanh(a - z * ase), math.tanh(a + z * ase)
        else:                       # covariance: VERIFY Stata's interval for cov()
            lo, hi = f0 - z * se, f0 + z * se
        rows.append((label, f0, se, lo, hi))

    pre = "sd" if sd else "var"
    if lv.cov == "identity":
        nm = lv.names[0] if q == 1 else "_cons"
        lab = f"{pre}({nm})" if q == 1 else f"{pre}({lv.names[0]}..{lv.names[-1]})"   # VERIFY
        add(lab, (lambda t: math.exp(t[0])) if sd else (lambda t: math.exp(2 * t[0])), "pos")
        return rows
    if lv.cov == "exchangeable":
        add(f"{pre}({lv.names[0]}..{lv.names[-1]})" if q > 1 else f"{pre}({lv.names[0]})",
            (lambda t: math.exp(t[0])) if sd else (lambda t: math.exp(2 * t[0])), "pos")
        if sd:
            add(f"corr({lv.names[0]}..{lv.names[-1]})", lambda t: math.tanh(t[1]), "corr")
        else:
            add(f"cov({lv.names[0]}..{lv.names[-1]})", lambda t: math.tanh(t[1]) * math.exp(2 * t[0]), "cov")
        return rows
    for i, nm in enumerate(lv.names):
        add(f"{pre}({nm})", (lambda t, i=i: math.exp(t[i])) if sd else (lambda t, i=i: math.exp(2 * t[i])), "pos")
    if lv.cov == "unstructured":
        k = q
        for i in range(q):
            for j in range(i + 1, q):
                if sd:
                    add(f"corr({lv.names[i]},{lv.names[j]})", lambda t, k=k: math.tanh(t[k]), "corr")
                else:
                    add(f"cov({lv.names[i]},{lv.names[j]})",
                        lambda t, i=i, j=j, k=k: math.tanh(t[k]) * math.exp(t[i] + t[j]), "cov")
                k += 1
    del S
    return rows


def _re_row(out, lab, est_, se, lo, hi) -> None:
    out.write(f"{lab[:28]:>28} |", "text")
    out.write(f"{_g(est_):>11}{_g(se):>11}{_g(lo):>14}{_g(hi):>12}\n", "result")


def re_table(s: "Session", est: Estimates, lev: float) -> None:
    out = s.output
    L = est.lmm
    model, th, Vth = L["model"], L["theta"], L["Vth"]
    z = st.norm.ppf(1 - (1 - lev / 100) / 2)
    sd = est.extra["sd"]
    out.write("\n" + "-" * 78 + "\n", "text")
    out.write(f"  Random-effects Parameters  |   Estimate   Std. Err.     [{lev:g}% Conf. Interval]\n", "text")
    out.write("-" * 29 + "+" + "-" * 48 + "\n", "text")
    pieces, lnse = model.split(th)
    k = 0
    for lv, t in zip(L["levels"], pieces):
        npar = lv.npar()
        Vl = Vth[k:k + npar, k:k + npar]
        k += npar
        kind = {"identity": "Identity", "independent": "Independent", "exchangeable": "Exchangeable",
                "unstructured": "Unstructured"}[lv.cov]
        out.write(f"{lv.var + ': ' + kind:<29}|\n", "text")
        for row in _vc_rows(lv, t, Vl, sd, z):
            _re_row(out, *row)
        out.write("-" * 29 + "+" + "-" * 48 + "\n", "text")
    se_l = float(np.sqrt(max(Vth[-1, -1], 0)))
    if sd:
        v = math.exp(lnse)
        se = v * se_l
        _re_row(out, "sd(Residual)", v, se, v * math.exp(-z * se_l), v * math.exp(z * se_l))
    else:
        v = math.exp(2 * lnse)
        se = 2 * v * se_l
        _re_row(out, "var(Residual)", v, se, v * math.exp(-2 * z * se_l), v * math.exp(2 * z * se_l))
    out.write("-" * 78 + "\n", "text")


# ---------------------------------------------------------------------------
# predict and estat
# ---------------------------------------------------------------------------

def _predict_mixed(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    from .postest import xb_all
    xb, ok = xb_all(s, est)
    if stat in ("", "xb"):
        return np.where(ok, xb, SYS), ("" if stat else "(option xb assumed)")   # VERIFY
    L = est.lmm
    model, th = L["model"], L["theta"]
    b = L["beta"]
    ds = s.data
    rows = np.nonzero(L["mask"])[0]
    u = model.blups(th, b)
    fitted = np.full(ds.nobs, SYS)
    for cl, uc in zip(model.clusters, u):
        zu = cl.Z @ uc
        fitted[rows[cl.rows]] = model.X[cl.rows] @ b + zu
    if stat == "fitted":
        return fitted, ""
    yv = ds.get(est.depvar).data.astype(float)
    if stat in ("residuals", "rstandard"):
        r = np.where(fitted < SYS, yv - fitted, SYS)
        if stat == "rstandard":
            _, lnse = model.split(th)
            r = np.where(r < SYS, r / math.exp(lnse), SYS)
        return r, ""
    raise StataError(198, f"option {stat} not allowed")


def _reffects(s: "Session", est: Estimates):
    """predict newvarlist, reffects: one variable per random effect, from the highest level
    down and, within a level, in the order of e(revars) ([ME] mixed postestimation)."""
    L = est.lmm
    model, th, b = L["model"], L["theta"], L["beta"]
    nreq = sum(lv.q for lv in model.levels)

    def fn(mask):
        ds = s.data
        rows = np.nonzero(L["mask"])[0]
        u = model.blups(th, b)
        out = [np.full(ds.nobs, SYS) for _ in range(nreq)]
        offs = np.cumsum([0] + [lv.q for lv in model.levels])
        for cl, uc in zip(model.clusters, u):
            k = 0
            for li, copies in cl.blocks:
                lv = model.levels[li]
                q = lv.q
                gids = np.unique(lv.gid[cl.rows])
                for c in range(copies):
                    m = lv.gid[cl.rows] == gids[c]
                    for j in range(q):
                        out[offs[li] + j][rows[cl.rows[m]]] = uc[k + j]
                    k += q
        return out
    return nreq, fn


def estat_icc(s: "Session", est: Estimates, rest: str) -> None:
    """Intraclass correlations: for level l, the share of the total variance of the random
    intercepts at l and above ([ME] estat icc). Only random-intercept models."""
    out = s.output
    L = getattr(est, "lmm", None)
    lev = est.level
    z = st.norm.ppf(1 - (1 - lev / 100) / 2)
    if L is not None:
        model, th, Vth = L["model"], L["theta"], L["Vth"]
        levels = L["levels"]
        resid = lambda t: math.exp(2 * t[-1])   # noqa: E731
    else:
        model = None
        levels = est.glmm["levels"]
        th, Vth = est.glmm["theta"], est.glmm["Vth"]
        resid = lambda t: est.glmm["resid_var"]   # noqa: E731
    # with covariates in the fixed part the correlation is conditional on them (compat 1001)
    has_cov = any(c.values is not None for c in (est.cols or []))
    title = "Residual intraclass correlation" if has_cov else "Intraclass correlation"
    if L is None and est.glmm["resid_var"] is None:
        raise StataError(321, "estat icc is not available after this model")   # VERIFY
    if any(lv.q != 1 or lv.names[0] != "_cons" for lv in levels):
        raise StataError(498, "estat icc is available only for random-intercept models")   # VERIFY
    out.write(f"\n{title}\n\n", "text")
    out.write("-" * 78 + "\n", "text")
    out.write(f"{'Level':>28} |        ICC   Std. Err.     [{lev:g}% Conf. Interval]\n", "text")
    out.write("-" * 29 + "+" + "-" * 48 + "\n", "text")
    nl = len(levels)
    res = {}
    for k in range(nl):
        def icc(t, k=k):
            v = [math.exp(2 * t[i]) for i in range(nl)]
            return sum(v[:k + 1]) / (sum(v) + resid(t))
        f0 = icc(th)
        J = np.zeros(th.size)
        for i in range(th.size):
            e = np.zeros(th.size)
            e[i] = 1e-6
            J[i] = (icc(th + e) - icc(th - e)) / 2e-6
        se = float(np.sqrt(max(J @ Vth @ J, 0)))
        lg = math.log(f0 / (1 - f0))
        lse = se / (f0 * (1 - f0))
        lo, hi = 1 / (1 + math.exp(-(lg - z * lse))), 1 / (1 + math.exp(-(lg + z * lse)))
        # the innermost level first: class|school (compat 1001)
        lab = "|".join(lv.var for lv in reversed(levels[:k + 1]))
        _re_row(out, lab, f0, se, lo, hi)
        res[k] = (f0, se)
    out.write("-" * 78 + "\n", "text")
    s.r = {"icc2": res[0][0]} if nl == 1 else {f"icc{k + 2}": v[0] for k, v in res.items()}
    del model


def estat_group(s: "Session", est: Estimates, rest: str) -> None:
    out = s.output
    ginfo = est.extra["groups"]
    out.write("\n" + "-" * 61 + "\n", "text")
    out.write("                |     No. of       Observations per Group\n", "text")
    out.write(" Group Variable |     Groups    Minimum    Average    Maximum\n", "text")
    out.write("----------------+--------------------------------------------\n", "text")
    for var, ng, mn, av, mx in ginfo:
        out.write(f"{var[:14]:>15} |", "text")
        out.write(f" {ng:>10,} {mn:>10,} {av:>10.1f} {mx:>10,}\n", "result")
    out.write("-" * 61 + "\n", "text")


def estat_recov(s: "Session", est: Estimates, rest: str) -> None:
    from ..commands.matrix import Matrix
    out = s.output
    L = est.lmm if getattr(est, "lmm", None) is not None else est.glmm
    model_levels = L["levels"]
    th = L["theta"]
    k = 0
    for lv in model_levels:
        npar = lv.npar()
        S = lv.sigma(th[k:k + npar])
        k += npar
        out.write(f"\nRandom-effects covariance matrix for level {lv.var}\n\n", "text")
        m = Matrix(S, list(lv.names), list(lv.names))
        from ..commands.matrix import list_matrix
        list_matrix(s, lv.var, m)
        s.r = {"cov": m}


def cmd_xtmixed(s: "Session", args: str) -> None:
    cmd_mixed(s, args, "xtmixed")

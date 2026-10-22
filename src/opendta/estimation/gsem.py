"""gsem ([SEM] gsem): generalized structural equation models by maximum likelihood.

Paths are written as in sem, with options per path group: (y <- x1 x2 L, family(bernoulli)
link(logit)). Names that start with a capital letter are latent; M[id] is a cluster-level
latent variable (random effect) of the clusters of variable id. Every response has a
conditional distribution from meglm (gaussian, bernoulli, binomial, poisson, nbinomial, gamma,
ordinal) whose linear predictor holds the covariates and the latent variables:

    eta_ri = x_ri b_r + sum_l lambda_rl z_l(unit i),    z ~ N(Xg, Sigma).

The likelihood integrates z out for every unit (an observation for observation-level latent
variables, a cluster for M[id]) by Gauss-Hermite quadrature, adaptive by default
(intmethod(mvaghermite), centered at the mode with the curvature there). Latent variables are
uncorrelated unless covariance() says otherwise; the first path of each latent variable is
fixed at 1. The standard errors come from the numerical observed information.

Undocumented details are marked VERIFY: the iteration log, the layout of ancillary parameters
and the names in e(b) (case compat/do/1909_other_suites.do). Observation-level latent variables
together with cluster-level ones (nested integrals) are not available."""

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
from .results import Estimates, crit, post

st = lazy("scipy.stats")
optimize = lazy("scipy.optimize")
special = lazy("scipy.special")
sparse = lazy("scipy.sparse")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

_FAMS = ("gaussian", "bernoulli", "binomial", "poisson", "nbinomial", "gamma", "ordinal")
_DEFAULT_LINK = {"gaussian": "identity", "bernoulli": "logit", "binomial": "logit",
                 "poisson": "log", "nbinomial": "log", "gamma": "log", "ordinal": "logit"}
_ALLOWED = {"gaussian": ("identity", "log"), "bernoulli": ("logit", "probit", "cloglog"),
            "binomial": ("logit", "probit", "cloglog"), "poisson": ("log",), "nbinomial": ("log",),
            "gamma": ("log",), "ordinal": ("logit", "probit")}
_SHORT = {"regress": ("gaussian", "identity"), "logit": ("bernoulli", "logit"),
          "probit": ("bernoulli", "probit"), "cloglog": ("bernoulli", "cloglog"),
          "poisson": ("poisson", "log"), "nbreg": ("nbinomial", "log"), "ologit": ("ordinal", "logit"),
          "oprobit": ("ordinal", "probit"), "gamma": ("gamma", "log")}
_FAM_TITLE = {"gaussian": "Gaussian", "bernoulli": "Bernoulli", "binomial": "Binomial",
              "gamma": "Gamma", "nbinomial": "Negative binomial", "ordinal": "Ordinal",
              "poisson": "Poisson"}                                   # VERIFY
_LINK_TITLE = {"identity": "Identity", "logit": "Logit", "probit": "Probit", "log": "Log",
               "cloglog": "Complementary log-log"}                    # VERIFY


# ---------------------------------------------------------------------------
# specification
# ---------------------------------------------------------------------------

@dataclass
class Eq:
    name: str
    kind: str                                  # "resp" or "lat"
    terms: list = field(default_factory=list)  # (kind, name, fixed, label): kind x | lat | clat | cons
    family: str | None = None
    link: str | None = None
    nbin: str | None = None
    disp: str = "mean"
    noconst: bool = False
    offset: str | None = None
    exposure: str | None = None


def _split_group(g: str):
    c = find_top(g, ",")
    return (g, "") if c == -1 else (g[:c], g[c + 1:])


def _opt_value(text: str, name: str):
    m = re.search(rf"\b{name}\(([^)]*)\)", text)
    return (m.group(1).strip(), text[:m.start()] + " " + text[m.end():]) if m else (None, text)


def _parse_spec(s: "Session", text: str, latent_opt: set[str], nocaps: bool):
    groups = re.findall(r"\(([^()]*(?:\([^()]*\)[^()]*)*)\)", text)
    rest = re.sub(r"\(([^()]*(?:\([^()]*\)[^()]*)*)\)", " ", text).strip()
    if rest or not groups:
        raise StataError(198, f"invalid path specification: {rest or text.strip()}")   # VERIFY
    eqs: dict[str, Eq] = {}
    latent: list[str] = []
    clat: list[tuple[str, str]] = []           # (latent, cluster variable)
    observed: list[str] = []

    def is_lat(name: str) -> bool:
        if name in latent_opt:
            return True
        return not nocaps and name[:1].isupper() and not s.data.has(name)

    def tok(t: str):
        name, at, val = t.partition("@")
        fixed = label = None
        if at:
            try:
                fixed = float(val)
            except ValueError:
                label = val
        return name, fixed, label

    def get(name: str, kind: str) -> Eq:
        if name not in eqs:
            eqs[name] = Eq(name, kind)
        return eqs[name]

    for g in groups:
        body, opts = _split_group(g)
        fam, v = _opt_value(opts, "family")
        lnk, opts2 = _opt_value(opts, "link")
        fam, opts2 = _opt_value(opts2, "family")
        offs, opts2 = _opt_value(opts2, "offset")
        expo, opts2 = _opt_value(opts2, "exposure")
        words = opts2.split()
        noconst = "noconstant" in words or "nocons" in words or "noconstan" in words
        short = [w for w in words if w in _SHORT]
        parts = re.split(r"\s*(<-|->)\s*", body.strip())
        if len(parts) != 3:
            raise StataError(198, f"invalid path ({body.strip()})")   # VERIFY
        left, arrow, right = parts
        src, dst = (right, left) if arrow == "<-" else (left, right)
        for dtok in dst.split():
            dname, dfx, dlb = tok(dtok)
            if "[" in dname:
                raise StataError(198, "a cluster-level latent variable cannot be an outcome")   # VERIFY
            dlat = is_lat(dname)
            if dlat and dname not in latent:
                latent.append(dname)
            if not dlat:
                for n in s.expand_varlist(dname):
                    if n not in observed:
                        observed.append(n)
            for dn in ([dname] if dlat else s.expand_varlist(dname)):
                e = get(dn, "lat" if dlat else "resp")
                if not dlat:
                    if fam is not None or short or lnk is not None:
                        fname, nb, disp = None, None, "mean"
                        if short:
                            fname, lk = _SHORT[short[0]]
                        else:
                            lk = None
                        if fam is not None:
                            w = fam.split()
                            hit = [f for f in _FAMS if f.startswith(w[0].lower())] if w else []
                            if w and w[0].lower() == "normal":
                                hit = ["gaussian"]
                            if len(hit) != 1:
                                raise StataError(198, f"family({fam}) invalid")   # VERIFY
                            fname = hit[0]
                            if fname == "binomial" and len(w) > 1:
                                nb = w[1]
                            if fname == "nbinomial" and len(w) > 1:
                                disp = "constant" if "constant".startswith(w[1].lower()) else "mean"
                        if lnk is not None:
                            lk = lnk.strip().lower()
                        e.family = fname or e.family or "gaussian"
                        e.link = lk if lk else _DEFAULT_LINK[e.family]
                        e.nbin, e.disp = nb, disp
                    e.noconst = e.noconst or noconst
                    e.offset = offs or e.offset
                    e.exposure = expo or e.exposure
                for stok in src.split():
                    sname, fx, lb = tok(stok)
                    fx = fx if fx is not None else dfx
                    lb = lb if lb is not None else dlb
                    if "[" in sname:
                        m = re.fullmatch(r"([A-Za-z_]\w*)\[(\w+)\]", sname)
                        if not m:
                            raise StataError(198, f"invalid latent variable {sname}")   # VERIFY
                        if dlat:
                            raise StataError(198, "a cluster-level latent variable cannot predict "
                                                  "a latent variable")   # VERIFY
                        if (m.group(1), m.group(2)) not in clat:
                            clat.append((m.group(1), m.group(2)))
                        e.terms.append(("clat", m.group(1), fx, lb))
                    elif is_lat(sname):
                        if dlat:
                            raise StataError(198, "paths between latent variables are not "
                                                  "available in gsem here")   # VERIFY
                        if sname not in latent:
                            latent.append(sname)
                        e.terms.append(("lat", sname, fx, lb))
                    else:
                        for n in s.expand_varlist(sname):
                            if s.data.get(n).is_string:
                                raise StataError(109, "type mismatch")
                            if n not in observed:
                                observed.append(n)
                            if not any(t[0] == "x" and t[1] == n for t in e.terms):
                                e.terms.append(("x", n, fx, lb))
    for e in eqs.values():
        if e.kind == "resp" and e.family is None:
            e.family, e.link = "gaussian", "identity"
        if e.kind == "resp" and e.link not in _ALLOWED[e.family]:
            raise StataError(198, f"link({e.link}) not allowed with family({e.family})")   # VERIFY
    return list(eqs.values()), latent, clat, observed


# ---------------------------------------------------------------------------
# model
# ---------------------------------------------------------------------------

class _Params:
    """Parameter table; a label shared by several entries makes them one parameter."""

    def __init__(self):
        self.eq, self.name, self.init, self.kind = [], [], [], []
        self.labels: dict[str, int] = {}
        self.slots: list[tuple[str, str, int]] = []      # (equation, name, parameter) as in e(b)

    def add(self, eq, name, init, kind="coef", label=None) -> int:
        if label is not None and label in self.labels:
            self.slots.append((eq, name, self.labels[label]))
            return self.labels[label]
        self.slots.append((eq, name, len(self.name)))
        self.eq.append(eq)
        self.name.append(name)
        self.init.append(float(init))
        self.kind.append(kind)
        if label is not None:
            self.labels[label] = len(self.name) - 1
        return len(self.name) - 1


@dataclass
class Resp:
    eq: Eq
    idx: np.ndarray                  # rows (of the sample) with the response observed
    y: np.ndarray
    X: np.ndarray
    off: np.ndarray
    fam: object
    naux: int
    ordered: bool
    xref: list                       # (theta index or -1, fixed value)
    cref: list                       # per latent dimension of the unit: (index or -1, fixed)
    aux: list
    S: object = None                 # unit x rows indicator
    levels: np.ndarray | None = None


class GModel:
    def __init__(self, s, eqs, latent, clat, observed, o, mask):
        from .meglm import _family, family_of, ordinal_family
        ds = s.data
        self.eqs = eqs
        resps = [e for e in eqs if e.kind == "resp"]
        if not resps:
            raise StataError(198, "no outcome variable specified")   # VERIFY
        if latent and clat:
            raise StataError(198, "observation-level and cluster-level latent variables in the "
                                  "same model are not available in OpenDTA yet")   # OpenDTA
        if len({c for _l, c in clat}) > 1:
            raise StataError(198, "nested cluster-level latent variables are not available "
                                  "in OpenDTA yet")   # OpenDTA
        self.cvar = clat[0][1] if clat else None
        self.lat = [c[0] for c in clat] if clat else list(latent)
        lat_regs = [e for e in eqs if e.kind == "lat"]
        q = len(self.lat)
        self.q = q
        n = ds.nobs
        # --- sample
        valid = []
        for e in resps:
            yv = ds.get(e.name).data
            ok = (yv < SYS) if yv.dtype.kind not in "OU" else np.zeros(n, bool)
            for t in e.terms:
                if t[0] == "x":
                    ok &= ds.get(t[1]).data < SYS
            for v in (e.offset, e.exposure):
                if v:
                    ok &= ds.get(v).data < SYS
            valid.append(ok)
        okany = np.any(valid, axis=0)
        if o.get("listwise"):
            okany = np.all(valid, axis=0)
        sample = mask & okany
        for e in lat_regs:
            for t in e.terms:
                if t[0] == "x":
                    sample &= ds.get(t[1]).data < SYS
        if self.cvar:
            if not ds.has(self.cvar):
                raise StataError(111, f"variable {self.cvar} not found")
            cv = ds.get(self.cvar)
            sample &= (cv.data < SYS) if not cv.is_string else (cv.data != "")
        if not sample.any():
            raise StataError(2000, "no observations")
        self.sample = sample
        rows = np.nonzero(sample)[0]
        self.N = len(rows)
        pos = -np.ones(n, dtype=int)
        pos[rows] = np.arange(self.N)
        if self.cvar:
            cvd = ds.get(self.cvar).data[rows]
            _u, gid = np.unique(cvd, return_inverse=True)
            self.unit = gid
            self.U = int(gid.max()) + 1
        else:
            self.unit = np.arange(self.N)
            self.U = self.N
        # --- parameters
        P = _Params()
        self.P = P
        self.var_fixed = o.get("_varfix", {})
        resp_objs: list[Resp] = []
        anchored: set[str] = set()
        first_seen: dict[str, tuple[int, int]] = {}
        for e in resps:
            yfull = ds.get(e.name).data
            sel = valid[resps.index(e)][rows] if False else (valid[resps.index(e)] & sample)
            idx = pos[np.nonzero(sel)[0]]
            y = yfull[sel].astype(float)
            xs = [t for t in e.terms if t[0] == "x"]
            cols = [ds.get(t[1]).data[sel].astype(float) for t in xs]
            cons = (not e.noconst) and e.family != "ordinal"
            X = np.column_stack(cols + ([np.ones(len(y))] if cons else [])) if (cols or cons) \
                else np.zeros((len(y), 0))
            off = np.zeros(len(y))
            if e.offset:
                off = off + ds.get(e.offset).data[sel].astype(float)
            if e.exposure:
                off = off + np.log(ds.get(e.exposure).data[sel].astype(float))
            fam, naux, ordered, levels = self._family(e, y, sel, ds, _family, family_of, ordinal_family)
            if e.family in ("bernoulli",):
                y = (y != 0).astype(float)
            xref = []
            for t in xs:
                if t[2] is not None:
                    xref.append((-1, t[2]))
                else:
                    xref.append((P.add(e.name, t[1], 0.0, "coef", t[3] and f"{t[3]}"), 0.0))
            if cons:
                xref.append((P.add(e.name, "_cons", 0.0), 0.0))
            cref = [(-2, 0.0)] * q          # -2: the latent variable is not in this equation
            for t in e.terms:
                if t[0] not in ("lat", "clat"):
                    continue
                j = self.lat.index(t[1])
                if t[2] is not None:
                    cref[j] = (-1, t[2])
                elif t[3] is not None:
                    cref[j] = (P.add(e.name, t[1], 1.0, "coef", t[3]), 0.0)
                else:
                    if t[1] not in anchored and self._anchor(t[1], o):
                        anchored.add(t[1])
                        cref[j] = (-1, 1.0)
                        first_seen[t[1]] = (len(resp_objs), j)
                    else:
                        cref[j] = (P.add(e.name, t[1], 1.0), 0.0)
            aux = []
            if naux:
                if e.family == "ordinal":
                    aux = [P.add(e.name, f"cut{j + 1}", 0.0, "cut") for j in range(naux)]
                elif e.family == "gaussian":
                    aux = [P.add(e.name, "var", 0.0, "gvar")]
                else:
                    nm = {"nbinomial": "lnalpha" if e.disp == "mean" else "lndelta", "gamma": "logs"}[e.family]
                    aux = [P.add(e.name, nm, 0.0, "aux")]
            r = Resp(e, idx, y, X, off, fam, naux, ordered, xref, cref, aux, levels=levels)
            if q:
                r.S = sparse.csr_matrix((np.ones(len(idx)), (self.unit[idx], np.arange(len(idx)))),
                                        shape=(self.U, len(idx)))
            resp_objs.append(r)
        self.resps = resp_objs
        # latent regressions (observation-level latent variables)
        self.gref: list[list[tuple[int, float]]] = [[] for _ in range(q)]
        self.gcols: list[list[np.ndarray]] = [[] for _ in range(q)]
        self.gnames: list[list[str]] = [[] for _ in range(q)]
        for e in lat_regs:
            j = self.lat.index(e.name) if e.name in self.lat else None
            if j is None or self.cvar:
                raise StataError(198, f"{e.name} is not a latent variable of the model")   # VERIFY
            for t in e.terms:
                if t[0] != "x":
                    continue
                self.gref[j].append((P.add(e.name, t[1], 0.0, "coef", t[3]), 0.0) if t[2] is None
                                    else (-1, t[2]))
                self.gcols[j].append(ds.get(t[1]).data[rows].astype(float))
                self.gnames[j].append(t[1])
        # latent variances and covariances
        self.vref = []
        for j, L in enumerate(self.lat):
            nm = f"{L}[{self.cvar}]" if self.cvar else L
            if L in self.var_fixed:
                self.vref.append((-1, float(self.var_fixed[L])))
            else:
                self.vref.append((P.add(f"var({nm})", "_cons", 0.0, "lvar"), 0.0))
        self.cref_cov = []
        for a, b in o.get("_covs", []):
            if a not in self.lat or b not in self.lat or a == b:
                raise StataError(198, f"invalid covariance({a}*{b})")   # VERIFY
            ia, ib = self.lat.index(a), self.lat.index(b)
            na = f"{a}[{self.cvar}]" if self.cvar else a
            nb = f"{b}[{self.cvar}]" if self.cvar else b
            self.cref_cov.append((ia, ib, P.add(f"cov({na},{nb})", "_cons", 0.0, "cov")))
        self.nint = int(o.get("intpoints") or 7)
        meth = str(o.get("intmethod") or "mvaghermite").lower()
        if meth == "laplace":
            raise StataError(198, "intmethod(laplace) is not available in OpenDTA yet")   # OpenDTA
        self.adaptive = meth in ("mvaghermite", "mcaghermite")
        self.intmethod = meth if meth in ("ghermite", "mvaghermite", "mcaghermite") else "mvaghermite"
        if q:
            from .meglm import _gh_nodes
            self.nodes, self.lw = _gh_nodes(self.nint, q)
            self.lphi = st.norm.logpdf(self.nodes).sum(axis=1)
        self.pm_cols = [np.column_stack(c) if c else None for c in self.gcols]
        self.z_cache = None
        self.anchors = first_seen
        self.n = len(P.name)

    def _anchor(self, L: str, o) -> bool:
        return L not in self.var_fixed

    def _family(self, e, y, sel, ds, _family, family_of, ordinal_family):
        f, lk = e.family, e.link
        if f == "ordinal":
            levels = np.unique(y)
            if len(levels) < 3:
                raise StataError(198, f"{e.name} has fewer than three categories; use a binary family")   # VERIFY
            yi = np.searchsorted(levels, y).astype(int)
            fn, naux = ordinal_family(lk, yi)
            return fn, naux, True, levels
        if f == "poisson":
            return _family("poisson", y), 0, False, None
        if f == "bernoulli":
            yb = (y != 0).astype(float)
            return family_of("binomial", lk, yb, nbin=np.ones(len(y)))[0], 0, False, None
        if f == "binomial":
            if e.nbin is None:
                raise StataError(198, "family(binomial) requires the number of trials")   # VERIFY
            try:
                nb = np.full(len(y), float(e.nbin))
            except ValueError:
                nb = ds.get(e.nbin).data[sel].astype(float)
            return family_of("binomial", lk, y, nbin=nb)[0], 0, False, None
        fn, naux = family_of(f, lk, y, disp=e.disp)
        return fn, naux, False, None

    # --- pieces of the likelihood -----------------------------------------
    @staticmethod
    def _v(th, ref):
        i, v = ref
        return th[i] if i >= 0 else v

    def sigma(self, th):
        q = self.q
        S = np.zeros((q, q))
        for j in range(q):
            S[j, j] = math.exp(th[self.vref[j][0]]) if self.vref[j][0] >= 0 else self.vref[j][1]
        for a, b, i in self.cref_cov:
            S[a, b] = S[b, a] = th[i]
        return S

    def pm(self, th):
        out = np.zeros((self.U, self.q))
        for j in range(self.q):
            if self.pm_cols[j] is not None:
                g = np.array([self._v(th, r) for r in self.gref[j]])
                out[:, j] = self.pm_cols[j] @ g
        return out

    def _eta_parts(self, th, r: Resp):
        beta = np.array([self._v(th, x) for x in r.xref])
        xb = (r.X @ beta if len(beta) else np.zeros(len(r.y))) + r.off
        c = np.array([self._v(th, x) for x in r.cref])
        a = np.array([th[i] for i in r.aux])
        return xb, c, a

    def _bad_cuts(self, th):
        return any(r.ordered and np.any(np.diff([th[i] for i in r.aux]) <= 0) for r in self.resps)

    def _aux_eff(self, r, a):
        if r.eq.family == "gaussian":
            return a                                     # ln variance
        return a

    def lf_points(self, th, Z):
        """Log density of the responses of every unit at the points Z (U x K x q): U x K."""
        K = Z.shape[1]
        lf = np.zeros((self.U, K))
        for r in self.resps:
            xb, c, a = self._eta_parts(th, r)
            eta = xb[:, None] + (Z[self.unit[r.idx]] @ c if self.q else 0.0)
            if eta.ndim == 1 or eta.shape[1] != K:
                eta = np.broadcast_to(eta, (len(xb), K)).copy()
            l, _s1, _s2 = r.fam(eta, a)
            l = np.nan_to_num(l, nan=-np.inf)
            if self.q:
                lf += r.S @ l
            else:
                lf[:, 0] += 0 if False else np.bincount(self.unit[r.idx], l[:, 0], minlength=self.U)
        return lf

    def _derivs(self, th, z):
        """Gradient and Hessian of the log density of the responses with respect to the latent
        variables of every unit, at the points z (U x q)."""
        g = np.zeros((self.U, self.q))
        H = np.zeros((self.U, self.q * self.q))
        for r in self.resps:
            xb, c, a = self._eta_parts(th, r)
            eta = (xb + np.einsum("nq,q->n", z[self.unit[r.idx]], c))[:, None]
            _l, s1, s2 = r.fam(eta, a)
            g += r.S @ (s1[:, 0][:, None] * c[None, :])
            H += r.S @ (s2[:, 0][:, None] * np.outer(c, c).ravel()[None, :])
        return g, H.reshape(self.U, self.q, self.q)

    def _objective(self, th, z, pm, Sinv):
        d = z - pm
        return self.lf_points(th, z[:, None, :])[:, 0] - 0.5 * np.einsum("ui,ij,uj->u", d, Sinv, d)

    def mode(self, th, S, pm, z0):
        Sinv = np.linalg.inv(S)
        z = z0.copy()
        f = self._objective(th, z, pm, Sinv)
        I = np.eye(self.q)
        for _ in range(100):
            g, H = self._derivs(th, z)
            g = g - (z - pm) @ Sinv
            Hn = -H + Sinv
            mn = np.linalg.eigvalsh(Hn)[:, 0]
            ridge = np.maximum(0.0, 1e-6 - mn)
            step = np.linalg.solve(Hn + ridge[:, None, None] * I, g[..., None])[..., 0]
            t = np.ones(self.U)
            for _h in range(30):
                zn = z + t[:, None] * step
                fn = self._objective(th, zn, pm, Sinv)
                bad = ~(fn >= f - 1e-12)
                if not bad.any():
                    break
                t = np.where(bad, t / 2, t)
            else:
                t = np.where(bad, 0.0, t)
                zn = z + t[:, None] * step
                fn = self._objective(th, zn, pm, Sinv)
            z, f = zn, np.where(np.isfinite(fn), fn, f)
            if np.max(np.abs(t[:, None] * step)) < 1e-11:
                break
        g, H = self._derivs(th, z)
        Hn = -H + Sinv
        return z, Hn

    def units(self, th):
        """(log likelihood of every unit, posterior weights, points) at th."""
        U, q = self.U, self.q
        if not q:
            lf = self.lf_points(th, np.zeros((U, 1, 0)))
            return lf[:, 0], None, None
        S = self.sigma(th)
        try:
            Sinv = np.linalg.inv(S)
            L = np.linalg.cholesky(S)
        except np.linalg.LinAlgError:
            return None, None, None
        _s, ldS = np.linalg.slogdet(S)
        pm = self.pm(th)
        if self.adaptive:
            z0 = self.z_cache if (self.z_cache is not None and self.z_cache.shape == pm.shape) else pm
            z, Hn = self.mode(th, S, pm, z0)
            self.z_cache = z
            try:
                C = np.linalg.cholesky(np.linalg.inv(Hn))
            except np.linalg.LinAlgError:
                return None, None, None
        else:
            z = pm
            C = np.broadcast_to(L, (U, q, q))
        Z = z[:, None, :] + np.einsum("uij,kj->uki", C, self.nodes)
        lf = self.lf_points(th, Z)
        d = Z - pm[:, None, :]
        quad = np.einsum("uki,ij,ukj->uk", d, Sinv, d)
        lprior = -0.5 * (q * math.log(2 * math.pi) + ldS + quad)
        ldC = np.linalg.slogdet(C)[1]
        lg = lf + lprior - self.lphi[None, :] + ldC[:, None] + self.lw[None, :]
        lL = special.logsumexp(lg, axis=1)
        return lL, lg, Z

    def loglik(self, th) -> float:
        if self._bad_cuts(th):
            return -np.inf
        try:
            lL, _lg, _Z = self.units(th)
        except (np.linalg.LinAlgError, FloatingPointError):
            return -np.inf
        if lL is None:
            return -np.inf
        v = float(lL.sum())
        return v if np.isfinite(v) else -np.inf

    def posterior_means(self, th):
        lL, lg, Z = self.units(th)
        pw = np.exp(lg - lL[:, None])
        return np.einsum("uk,ukq->uq", pw, Z)


# ---------------------------------------------------------------------------
# fitting
# ---------------------------------------------------------------------------

def _num_grad(f, x, h0=1e-5):
    g = np.zeros(len(x))
    for i in range(len(x)):
        h = h0 * (1 + abs(x[i]))
        xp, xm = x.copy(), x.copy()
        xp[i] += h
        xm[i] -= h
        g[i] = (f(xp) - f(xm)) / (2 * h)
    return g


def _num_hess(gfun, x, h0=1e-4):
    p = len(x)
    H = np.zeros((p, p))
    for i in range(p):
        h = h0 * (1 + abs(x[i]))
        xp, xm = x.copy(), x.copy()
        xp[i] += h
        xm[i] -= h
        H[:, i] = (gfun(xp) - gfun(xm)) / (2 * h)
    return (H + H.T) / 2


def _start(m: GModel, s, log):
    """Starting values: every response fitted alone without latent variables."""
    th = np.array(m.P.init, float)
    for r in m.resps:
        free = [i for i, _v in r.xref if i >= 0]
        nx = len(r.xref)
        aux = r.aux
        beta0 = np.zeros(nx)
        ncut = len(aux) if r.ordered else 0
        if r.ordered:
            p = np.array([np.mean(r.y <= v) for v in r.levels[:-1]])
            p = np.clip(p, 1e-3, 1 - 1e-3)
            cuts = np.log(p / (1 - p)) if r.eq.link == "logit" else st.norm.ppf(p)
            cuts = np.maximum.accumulate(cuts + 1e-3 * np.arange(len(cuts)))
            a0 = np.concatenate([[cuts[0]], np.log(np.maximum(np.diff(cuts), 1e-3))])
        elif r.eq.family == "gaussian":
            a0 = np.array([math.log(max(float(np.var(r.y)), 1e-8))])
        elif r.eq.family in ("nbinomial", "gamma"):
            a0 = np.array([0.0])
        else:
            a0 = np.zeros(0)
        if r.eq.family in ("poisson", "nbinomial", "gamma") or (r.eq.family == "gaussian" and r.eq.link == "log"):
            if r.xref and r.xref[-1][0] >= 0 and r.eq.family and not r.eq.noconst:
                beta0[-1] = math.log(max(float(np.mean(r.y)), 0.1))
        def unpack(v, nx=nx, ncut=ncut, r=r):
            b = v[:nx]
            a = v[nx:]
            if r.ordered:
                a = np.cumsum(np.concatenate([[a[0]], np.exp(a[1:])]))
            return b, a

        def nll(v, r=r, unpack=unpack):
            b, a = unpack(v)
            bv = np.array([th[i] if False else 0.0 for i in range(0)])
            beta = np.array([(b[k] if r.xref[k][0] >= 0 else r.xref[k][1]) for k in range(nx)])
            eta = ((r.X @ beta) if nx else np.zeros(len(r.y))) + r.off
            l = r.fam(eta[:, None], a)[0]
            val = -float(np.sum(l))
            return val if np.isfinite(val) else 1e300
        v0 = np.concatenate([beta0, a0])
        try:
            res = optimize.minimize(nll, v0, method="BFGS", options={"gtol": 1e-6, "maxiter": 200})
            v = res.x
        except Exception:
            v = v0
        b, a = unpack(v)
        for k, (i, _f) in enumerate(r.xref):
            if i >= 0:
                th[i] = b[k]
        for k, i in enumerate(aux):
            th[i] = a[k]
    for i, kd in enumerate(m.P.kind):
        if kd == "lvar":
            th[i] = math.log(0.5)
        elif kd == "cov":
            th[i] = 0.0
    return th


def fit(m: GModel, s, log: bool, maxiter: int):
    th = _start(m, s, log)
    out = s.output
    f = lambda x: -m.loglik(x)                                    # noqa: E731
    gfun = lambda x: _num_grad(f, x)                              # noqa: E731
    if log:
        out.write("\nFitting full model:\n\n", "text")      # VERIFY
    ll0 = m.loglik(th)
    if not np.isfinite(ll0):
        raise StataError(430, "initial values not feasible")
    state = {"k": 0}
    if log:
        out.write(f"Iteration 0:   log likelihood = {ll0:.7f}\n".replace("= -", "= -"), "text")

    def cb(x):
        state["k"] += 1
        if log:
            out.write(f"Iteration {state['k']}:   log likelihood = {m.loglik(x):.7f}\n", "text")
    res = optimize.minimize(f, th, jac=gfun, method="BFGS", callback=cb,
                            options={"gtol": 1e-7, "maxiter": maxiter})
    x = res.x
    # Newton polishing with the numerical Hessian
    cur = -f(x)
    for _ in range(20):
        g = -gfun(x)
        H = _num_hess(lambda y: -gfun(y), x)                      # Hessian of ll
        try:
            step = np.linalg.solve(-H + 1e-10 * np.eye(len(x)), g)
        except np.linalg.LinAlgError:
            break
        t = 1.0
        for _h in range(20):
            new = -f(x + t * step)
            if new >= cur - 1e-12:
                break
            t /= 2
        else:
            break
        xn = x + t * step
        d = new - cur
        x, cur = xn, new
        if abs(d) < 1e-11 * (1 + abs(cur)):
            break
    g = gfun(x)
    Hll = _num_hess(lambda y: -gfun(y), x)
    try:
        V = np.linalg.inv(-Hll)
        V = (V + V.T) / 2
    except np.linalg.LinAlgError:
        V = np.full((len(x), len(x)), np.nan)
    conv = bool(np.max(np.abs(g)) < 1e-4 * (1 + abs(cur)) and np.all(np.isfinite(V)))
    return x, V, cur, conv, state["k"]


# ---------------------------------------------------------------------------
# command
# ---------------------------------------------------------------------------

_OPTS = {"nocapslatent": 9, "latent": 4, "covariance": 3, "variance": 3, "nolog": 4, "log": 3,
         "level": 1, "iterate": 4, "noheader": 6, "nodvheader": 5, "nofootnote": 7, "notable": 5,
         "coeflegend": 5, "listwise": 4, "intmethod": 4, "intpoints": 4, "vce": 3, "from": 4,
         "nodescribe": 6, "noconstant": 4, "nocnsreport": 5, "ntol": 3, "tolerance": 3,
         "ltolerance": 4, "nrtolerance": 4, "technique": 4, "difficult": 4}


def _parse_options(args: str):
    c = find_top(args, ",")
    head, opts = (args, "") if c == -1 else (args[:c], args[c + 1:])
    m = re.search(r"\)\s*((?:if|in)\b.*)$", head, re.S)
    paths_text, tail = (head[:m.start() + 1], m.group(1)) if m else (head, "")
    p = parse_standard(tail + ("," + opts if opts.strip() else ""))
    return paths_text, p, match_options(p.options, _OPTS)


def cmd_gsem(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "gsem", args)
        return
    from ..commands._util import touse
    paths_text, p, o = _parse_options(args)
    vce = str(o.get("vce") or "oim").strip().lower()
    if vce not in ("oim", "oi", "o"):
        raise StataError(198, f"vce({vce}) is not available in OpenDTA yet")   # OpenDTA
    latent_opt = set(str(o["latent"]).split()) if isinstance(o.get("latent"), str) else set()
    eqs, latent, clat, observed = _parse_spec(s, paths_text, latent_opt, bool(o.get("nocapslatent")))
    covs = []
    for item in str(o.get("covariance") or "").split():
        a, star, b = item.partition("*")
        if not star:
            raise StataError(198, f"invalid covariance() specification: {item}")   # VERIFY
        covs.append((a, b))
    var_fix = {}
    for item in str(o.get("variance") or "").split():
        name, at, val = item.partition("@")
        if at:
            var_fix[name] = float(val)
    o["_covs"], o["_varfix"] = covs, var_fix
    mask = touse(s, p)
    m = GModel(s, eqs, latent, clat, observed, o, mask)
    log = not o.get("nolog")
    maxiter = int(o["iterate"]) if o.get("iterate") else 300
    th, V, ll, conv, iters = fit(m, s, log, maxiter)
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", 95))
    _post(s, args, m, th, V, ll, conv, iters, level)


def _report(m: GModel, th, V):
    """Parameters in the reported scale (variances instead of their logs) and their covariance."""
    rep = th.copy()
    D = np.ones(len(th))
    for i, kd in enumerate(m.P.kind):
        if kd in ("lvar", "gvar"):
            rep[i] = math.exp(th[i])
            D[i] = rep[i]
    return rep, D[:, None] * V * D[None, :]


def _post(s, args, m: GModel, th, V, ll, conv, iters, level):
    rep, Vr = _report(m, th, V)
    P = m.P
    names, eqn, ix = [], [], []
    for e, nm, i in P.slots:
        ix.append(i)
        if P.kind[i] == "gvar":
            eqn.append(f"var(e.{e})")
            names.append("_cons")
        elif P.kind[i] in ("cut", "aux"):
            eqn.append(f"/{nm}")
            names.append("_cons")
        else:
            eqn.append(e)
            names.append(nm)
    ix = np.array(ix)
    est = Estimates("gsem", " ".join(e.eq.name for e in m.resps), names, rep[ix].copy(),
                    Vr[np.ix_(ix, ix)], m.N, "z", None)
    est.eqnames = eqn
    est.level = level
    est.model = m
    est.extra = dict(N=m.N, ll=ll, th=th.copy(), V=V, rep=rep, Vr=Vr)
    est.display = display_gsem
    est.predict = _predict
    est.predict_opts = {"mu": 2, "xb": 2, "latent": 3}
    resp_names = [r.eq.name for r in m.resps]
    scalars = [("N", float(m.N)), ("k", float(len(th))), ("k_eq", float(len(m.resps))),
               ("k_dv", float(len(m.resps))), ("converged", 1.0 if conv else 0.0), ("rc", 0.0),
               ("ic", float(iters)), ("ll", ll), ("rank", float(np.linalg.matrix_rank(est.V)))]
    macros = [("cmdline", f"gsem {args.strip()}"), ("cmd", "gsem"), ("predict", "gsem_p"),
              ("title", "Generalized structural equation model"), ("vce", "oim"),
              ("vcetype", "OIM"), ("depvar", " ".join(resp_names)),
              ("intmethod", m.intmethod if m.q else ""), ("properties", "b V")]
    for k, r in enumerate(m.resps, 1):
        macros.append((f"family{k}", r.eq.family))
        macros.append((f"link{k}", r.eq.link))
    if m.q:
        scalars.insert(4, ("n_quad", float(m.nint)))
    post(s, est, scalars, macros, m.sample)
    display_gsem(s, est)


def _g(x, fmt="%9.0g"):
    from ..core.formats import format_value
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def _constraint_lines(m: GModel):
    lines = []
    for r in m.resps:
        nm = r.eq.name
        for t in r.eq.terms:
            pass
        for j, (i, v) in enumerate(r.cref):
            if i == -1:
                L = m.lat[j]
                lab = f"{L}[{m.cvar}]" if m.cvar else L
                lines.append(f"[{nm}]{lab} = {v:g}")
        for k, (i, v) in enumerate(r.xref):
            if i < 0:
                xs = [t for t in r.eq.terms if t[0] == "x"]
                lines.append(f"[{nm}]{xs[k][1]} = {v:g}")
    for j, (i, v) in enumerate(m.vref):
        if i < 0:
            lines.append(f"var({m.lat[j]}) = {v:g}")
    # shared labels
    by_param: dict[int, list[str]] = {}
    for r in m.resps:
        for j, (i, _v) in enumerate(r.cref):
            if i >= 0:
                by_param.setdefault(i, []).append(f"[{r.eq.name}]{m.lat[j]}")
        xs = [t for t in r.eq.terms if t[0] == "x"]
        for k, (i, _v) in enumerate(r.xref):
            if i >= 0 and k < len(xs):
                by_param.setdefault(i, []).append(f"[{r.eq.name}]{xs[k][1]}")
    for i, items in by_param.items():
        for other in items[1:]:
            lines.append(f"{items[0]} - {other} = 0")
    return lines


def display_gsem(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                 level: float | None = None, **_kw) -> None:
    out = s.output
    m: GModel = est.model
    x = est.extra
    lv = level or est.level
    if header:
        out.write("\n" + f"{'Generalized structural equation model':<48}{'Number of obs':<18}= ", "text")
        out.write(f"{x['N']:>10,}\n", "result")
        for r in m.resps:
            e = r.eq
            fam = _FAM_TITLE[e.family]
            out.write(f"\n{'Response':<15}: {e.name}\n{'Family':<15}: {fam}\n"
                      f"{'Link':<15}: {_LINK_TITLE[e.link]}\n", "text")
        out.write("\n", "text")
        out.write("Log likelihood = ", "text")
        out.write(f"{_g(x['ll'], '%10.0g')}\n", "result")
        lines = _constraint_lines(m)
        if lines:
            out.write("\n", "text")
        for k, ln in enumerate(lines, 1):
            out.write(f" ({k:2d})  {ln}\n", "text")
        if not lines:
            out.write("\n", "text")
    if not table:
        return
    rep, Vr = x["rep"], x["Vr"]
    se = np.sqrt(np.maximum(np.diag(Vr), 0))
    c = crit("z", lv, None)
    P = m.P
    labels = []
    for i, kd in enumerate(P.kind):
        if kd in ("lvar",):
            labels.append(f"var({P.eq[i][4:-1]})")
        elif kd == "gvar":
            labels.append(f"var(e.{P.eq[i]})")
        elif kd == "cov":
            labels.append(P.eq[i])
        elif kd in ("cut", "aux"):
            labels.append(f"/{P.name[i]}")
    W = max([12] + [len(lb) - 1 for lb in labels])
    dash = "-" * (W + 1) + "+" + "-" * 64 + "\n"
    out.write("-" * (W + 66) + "\n", "text")
    out.write(f"{'':>{W}} |                 OIM\n", "text")
    out.write(f"{'':>{W}} |      Coef.   Std. Err.      z    P>|z|     [{lv:g}% Conf. Interval]\n", "text")
    out.write(dash, "text")

    def row(label, b, s_, kind="coef"):
        if kind in ("var", "cov", "aux"):
            out.write(f"{label:>{W + 1}}|", "text")
        else:
            out.write(f"{label:>{W}} |", "text")
        if kind == "constrained":
            out.write(f"  {_g(b):>9}  (constrained)\n", "result")
            return
        if s_ > 0 and np.isfinite(s_):
            if kind == "var":
                try:
                    lo, hi = b * math.exp(-c * s_ / b), b * math.exp(c * s_ / b)
                except OverflowError:
                    lo, hi = float("nan"), float("nan")
                out.write(f"  {_g(b):>9}  {_g(s_):>9}" + " " * 17 + f"    {_g(lo):>9}   {_g(hi):>9}\n",
                          "result")
                return
            if kind == "aux":
                out.write(f"  {_g(b):>9}  {_g(s_):>9}" + " " * 17 + f"    {_g(b - c * s_):>9}"
                          f"   {_g(b + c * s_):>9}\n", "result")
                return
            z = b / s_
            pv = 2 * float(st.norm.sf(abs(z)))
            out.write(f"  {_g(b):>9}  {_g(s_):>9}  {z:>7.2f}  {pv:>6.3f}    {_g(b - c * s_):>9}"
                      f"   {_g(b + c * s_):>9}\n", "result")
        elif kind == "var":
            out.write(f"  {_g(b):>9}  {'.':>9}" + " " * 17 + f"    {'.':>9}   {'.':>9}\n", "result")
        else:
            out.write(f"  {_g(b):>9}  {'.':>9}\n", "result")

    first = True
    for e in m.eqs:
        if not first:
            out.write(dash, "text")
        first = False
        out.write(f"{e.name + ' <-':<{W}} |\n", "text")
        if e.kind == "resp":
            r = next(rr for rr in m.resps if rr.eq is e)
            xi = 0
            for t in e.terms:
                if t[0] == "x":
                    i, v = r.xref[xi]
                    xi += 1
                    row(t[1], rep[i], se[i]) if i >= 0 else row(t[1], v, 0, "constrained")
                else:
                    j = m.lat.index(t[1])
                    i, v = r.cref[j]
                    lab = f"{t[1]}[{m.cvar}]" if t[0] == "clat" else t[1]
                    row(lab, rep[i], se[i]) if i >= 0 else row(lab, v, 0, "constrained")
            if len(r.xref) > xi:
                i = r.xref[xi][0]
                row("_cons", rep[i], se[i])
        else:
            j = m.lat.index(e.name)
            for k, nm in enumerate(m.gnames[j]):
                i, v = m.gref[j][k]
                row(nm, rep[i], se[i]) if i >= 0 else row(nm, v, 0, "constrained")
    rows_after = [i for i, kd in enumerate(P.kind) if kd == "lvar"]
    for i in rows_after:
        out.write(dash, "text") if i == rows_after[0] else None
        row(f"var({P.eq[i][4:-1]})", rep[i], se[i], "var")
    for j, (i, v) in enumerate(m.vref):
        if i < 0:
            row(f"var({m.lat[j]})", v, 0, "constrained")
    covs = [i for i, kd in enumerate(P.kind) if kd == "cov"]
    for i in covs:
        out.write(dash, "text") if i == covs[0] else None
        row(P.eq[i], rep[i], se[i], "cov")
    aux = [i for i, kd in enumerate(P.kind) if kd in ("cut", "aux", "gvar")]
    for n_, i in enumerate(aux):
        out.write(dash, "text") if n_ == 0 else None
        if P.kind[i] == "gvar":
            row(f"var(e.{P.eq[i]})", rep[i], se[i], "var")
        else:
            row(f"/{P.name[i]}", rep[i], se[i], "aux")
    out.write("-" * (W + 66) + "\n", "text")


def _predict(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    """predict newvar, mu(y) | xb(y) | latent(L): the mean and the linear predictor of a
    response (latent variables at their empirical Bayes means) and the latent means."""
    m: GModel = est.model
    th = est.extra["th"]
    target = (getattr(est, "_pred_opts", None) or {}).get(stat)
    target = target.strip() if isinstance(target, str) and target.strip() else None
    n = s.data.nobs
    rows = np.nonzero(m.sample)[0]
    out = np.full(n, SYS)
    eb = m.posterior_means(th) if m.q else np.zeros((m.U, 0))
    if stat == "latent":
        L = target or (m.lat[0] if m.lat else None)
        if L not in m.lat:
            raise StataError(198, f"{L} is not a latent variable")   # VERIFY
        out[rows] = eb[m.unit, m.lat.index(L)]
        return out, ""
    r = next((rr for rr in m.resps if rr.eq.name == target), None) if target else m.resps[0]
    if r is None:
        raise StataError(198, f"{target} is not an outcome of the model")   # VERIFY
    ds = s.data
    e = r.eq
    cols = [ds.get(t[1]).data[rows].astype(float) for t in e.terms if t[0] == "x"]
    beta = np.array([m._v(th, x) for x in r.xref])
    Xf = np.column_stack(cols + ([np.ones(m.N)] if len(r.xref) > len(cols) else [])) \
        if (cols or len(r.xref) > len(cols)) else np.zeros((m.N, 0))
    eta = (Xf @ beta if len(beta) else np.zeros(m.N))
    for v in (e.offset,):
        if v:
            eta = eta + ds.get(v).data[rows]
    if e.exposure:
        eta = eta + np.log(ds.get(e.exposure).data[rows])
    if stat in ("mu", "") and m.q:
        c = np.array([m._v(th, x) for x in r.cref])
        eta = eta + eb[m.unit] @ c
    if stat == "xb" or e.link == "identity":
        val = eta
    elif e.link == "log":
        val = np.exp(eta)
    elif e.link == "logit":
        val = special.expit(eta)
    elif e.link == "probit":
        val = special.ndtr(eta)
    else:
        val = -np.expm1(-np.exp(eta))
    if e.family == "ordinal" and stat != "xb":
        raise StataError(198, "predict mu is not available for ordinal outcomes")   # OpenDTA
    out[rows] = val
    return out, ""

"""arch: autoregressive conditional heteroskedasticity models ([TS] arch). The log likelihood and
the variance recursions are plain functions of arrays (ArchSpec describes the model); cmd_arch
parses the syntax, reads the data, maximizes (BHHH then BFGS, as arima) and prints."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from scipy import special

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.syntax import find_top, match_options, parse_standard
from . import fvars
from .ml import maximize_switching, print_log
from .ml import options as ml_options
from .results import CoefRow, Estimates, coef_table, post, rows_from_columns
from .tsops import require_ts

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

# variance terms in the order Stata lists them; (option, number of extra parameters per lag)
TERMS = ["arch", "garch", "saarch", "tarch", "aarch", "narch", "narchk", "abarch", "atarch", "sdgarch",
         "earch", "egarch", "parch", "tparch", "aparch", "nparch", "nparchk", "pgarch"]
POWER = {"parch", "tparch", "aparch", "nparch", "nparchk", "pgarch"}
EXTRA = {"aarch": "aarch_e", "narch": "narch_k", "narchk": "narch_k", "earch": "earch_a",
         "aparch": "aparch_e", "nparch": "nparch_k", "nparchk": "nparch_k"}


@dataclass
class ArchSpec:
    ar: list[int] = field(default_factory=list)
    ma: list[int] = field(default_factory=list)
    terms: list[tuple[str, int]] = field(default_factory=list)       # (kind, lag)
    archm: list[int] = field(default_factory=list)                  # lags of the variance in the mean (0 = current)
    archm_fn: object = None                                          # g(sigma^2) for the mean (default identity)
    dist: str = "gaussian"
    dist_fixed: float | None = None                                  # fixed df or shape
    nhet: int = 0                                                    # columns of the het() matrix (with constant)
    arch0: str = "xb"
    arch0_value: float = 0.0
    condobs: int = 0
    # parameter layout, filled by layout()
    k: int = 0

    @property
    def power(self) -> bool:
        return any(kd in POWER for kd, _ in self.terms)

    @property
    def egarch(self) -> bool:
        return any(kd in ("earch", "egarch") for kd, _ in self.terms)

    def layout(self, k: int) -> dict:
        """Positions of each group of parameters in the vector: beta, psi, ar, ma, het, terms,
        gamma0, power, dist."""
        pos, q = {}, 0
        pos["beta"] = (q, q + k); q += k
        pos["psi"] = (q, q + len(self.archm)); q += len(self.archm)
        pos["ar"] = (q, q + len(self.ar)); q += len(self.ar)
        pos["ma"] = (q, q + len(self.ma)); q += len(self.ma)
        pos["het"] = (q, q + self.nhet); q += self.nhet
        pos["terms"] = []
        for kd, lag in self.terms:
            n_ = 2 if kd in ("aarch", "earch", "aparch") else 1
            if kd in ("narch", "nparch"):
                n_ = 2
            if kd in ("narchk", "nparchk"):
                n_ = 2                                  # alpha and the common kappa
            pos["terms"].append((q, q + n_)); q += n_
        pos["gamma0"] = (q, q + (0 if self.nhet else 1)); q += (0 if self.nhet else 1)
        pos["power"] = (q, q + (1 if self.power else 0)); q += (1 if self.power else 0)
        pos["dist"] = (q, q + (1 if self.dist != "gaussian" and self.dist_fixed is None else 0))
        q = pos["dist"][1]
        pos["n"] = q
        return pos


def _prime(spec: ArchSpec, w: np.ndarray, w0: float) -> float:
    """Presample value of sigma^2 and epsilon^2."""
    m = spec.arch0
    if m == "zero":
        return 0.0
    if m == "value":
        return spec.arch0_value
    if m in ("xbwt", "xb0wt"):
        wt = 0.7 ** np.arange(len(w))
        return float(np.sum(wt * w * w) / np.sum(wt))
    return float(np.mean(w * w))


def archpath(spec: ArchSpec, par: np.ndarray, y: np.ndarray, X: np.ndarray, Z: np.ndarray | None,
             pos: dict, prime0: float | None = None):
    """Conditional variance, innovations and mean for the whole sample. Returns None when the
    parameters are infeasible. Output: dict(s2, e, mu, h)."""
    n = len(y)
    beta = par[pos["beta"][0]:pos["beta"][1]]
    psi = par[pos["psi"][0]:pos["psi"][1]]
    rho = par[pos["ar"][0]:pos["ar"][1]]
    theta = par[pos["ma"][0]:pos["ma"][1]]
    xb = X @ beta if X.shape[1] else np.zeros(n)
    w0 = y - xb
    prime = prime0 if prime0 is not None else _prime(spec, w0, 0.0)
    if spec.arch0 == "zero":
        prime = 0.0
    phi = par[pos["power"][0]] if spec.power else 2.0
    if spec.power and (phi <= 0 or not np.isfinite(phi)):
        return None
    if spec.nhet:
        lam = par[pos["het"][0]:pos["het"][1]]
        g0 = np.exp(Z @ lam)
    else:
        g0 = np.full(n, par[pos["gamma0"][0]])
    maxl = max([lag for _k, lag in spec.terms] + spec.archm + [0])
    # state arrays with presample padding of maxl entries
    P = max(maxl, 1)
    E = np.zeros(n + P)                      # signed innovations (presample 0)
    AE = np.full(n + P, math.sqrt(max(prime, 0.0)))
    S2 = np.full(n + P, prime)
    SD = np.full(n + P, math.sqrt(max(prime, 0.0)))
    LS = np.full(n + P, math.log(prime) if prime > 0 else 0.0)
    ZZ = np.zeros(n + P)
    Wd = np.zeros(n + P)                     # disturbances w_t = y - xb - archm (ARMA part)
    mu = np.zeros(n)
    h = np.zeros(n)
    egarch, power = spec.egarch, spec.power
    sq2pi = math.sqrt(2 / math.pi)
    tp = [par[a:b] for a, b in pos["terms"]]
    gfun = spec.archm_fn
    for t in range(n):
        i = t + P
        v = g0[t]
        for (kd, lag), p_ in zip(spec.terms, tp):
            j = i - lag
            e, ae = E[j], AE[j]
            if kd == "arch":
                v += p_[0] * ae * ae
            elif kd == "garch":
                v += p_[0] * S2[j]
            elif kd == "saarch":
                v += p_[0] * e
            elif kd == "tarch":
                v += p_[0] * ae * ae * (e > 0)
            elif kd == "aarch":
                v += p_[0] * (ae + p_[1] * e) ** 2
            elif kd == "narch":
                v += p_[0] * (e - p_[1]) ** 2
            elif kd == "narchk":
                v += p_[0] * (e - p_[1]) ** 2
            elif kd == "abarch":
                v += p_[0] * ae
            elif kd == "atarch":
                v += p_[0] * ae * (e > 0)
            elif kd == "sdgarch":
                v += p_[0] * SD[j]
            elif kd == "earch":
                zj = ZZ[j]
                v += p_[0] * zj + p_[1] * (abs(zj) - sq2pi)
            elif kd == "egarch":
                v += p_[0] * LS[j]
            elif kd == "parch":
                v += p_[0] * ae ** phi
            elif kd == "tparch":
                v += p_[0] * (e > 0) * ae ** phi
            elif kd == "aparch":
                base = ae + p_[1] * e
                v += p_[0] * (base ** phi if base > 0 else 0.0)
            elif kd in ("nparch", "nparchk"):
                v += p_[0] * abs(e - p_[1]) ** phi
            elif kd == "pgarch":
                v += p_[0] * SD[j] ** phi
        h[t] = v
        if egarch:
            if not np.isfinite(v) or abs(v) > 700:
                return None
            s2 = math.exp(v)
        elif power:
            if not v > 0 or not np.isfinite(v):
                return None
            s2 = v ** (2.0 / phi)
        else:
            if not v > 0 or not np.isfinite(v):
                return None
            s2 = v
        S2[i], SD[i], LS[i] = s2, math.sqrt(s2), math.log(s2)
        m = xb[t]
        for ps, lag in zip(psi, spec.archm):
            sv = S2[i - lag]
            m += ps * (gfun(sv) if gfun is not None else sv)
        w = y[t] - m
        eps = w
        for a_, lag in zip(rho, spec.ar):
            if t - lag >= 0:
                eps -= a_ * Wd[i - lag]
        for b_, lag in zip(theta, spec.ma):
            if t - lag >= 0:
                eps -= b_ * E[i - lag]
        Wd[i] = w
        E[i], AE[i], ZZ[i] = eps, abs(eps), eps / SD[i]
        mu[t] = m
    return {"s2": S2[P:], "e": E[P:], "mu": mu, "h": h}


def archloglik(spec: ArchSpec, par: np.ndarray, y, X, Z, pos, prime0=None) -> np.ndarray | None:
    """Log-likelihood contributions (None when infeasible)."""
    pth = archpath(spec, par, y, X, Z, pos, prime0)
    if pth is None:
        return None
    s2, e = pth["s2"], pth["e"]
    if not np.all(np.isfinite(s2)) or np.any(s2 <= 0):
        return None
    if spec.dist == "gaussian":
        ll = -0.5 * (np.log(2 * math.pi * s2) + e * e / s2)
    elif spec.dist == "t":
        df = spec.dist_fixed if spec.dist_fixed is not None else math.exp(min(par[pos["dist"][0]], 60.0)) + 2
        ll = (special.gammaln((df + 1) / 2) - special.gammaln(df / 2) - 0.5 * np.log(math.pi * (df - 2))
              - 0.5 * np.log(s2) - (df + 1) / 2 * np.log1p(e * e / (s2 * (df - 2))))
    else:
        sh = spec.dist_fixed if spec.dist_fixed is not None else math.exp(min(par[pos["dist"][0]], 60.0))
        lam = math.sqrt(2.0 ** (-2 / sh) * special.gamma(1 / sh) / special.gamma(3 / sh))
        ll = (math.log(sh / lam) - (1 + 1 / sh) * math.log(2) - special.gammaln(1 / sh)
              - 0.5 * np.abs(e / (np.sqrt(s2) * lam)) ** sh - 0.5 * np.log(s2))
    if spec.condobs:
        ll = ll.copy()
        ll[:spec.condobs] = 0.0
    return ll if np.all(np.isfinite(ll)) else None


# ---------------------------------------------------------------------------
# command
# ---------------------------------------------------------------------------

_OPTS = {"noconstant": 4, "archm": 5, "archmlags": 9, "archmexp": 8, "arima": 5, "ar": 2, "ma": 2,
         "distribution": 4, "het": 3, "savespace": 5, "arch0": 5, "arma0": 5, "condobs": 5,
         "vce": 3, "robust": 1, "level": 1, "detail": 3, "nocnsreport": 5, "nolog": 3, "iterate": 4,
         "tolerance": 3, "ltolerance": 4, "nrtolerance": 4, "technique": 4, "gtolerance": 3,
         "difficult": 4, "from": 3, "constraints": 3, "collinear": 4, "trace": 5}
_OPTS.update({"arch": 4, "garch": 5, "aarch": 3, "aparch": 3, "abarch": 3, "atarch": 3, "narch": 5,
              "narchk": 6, "nparch": 5, "nparchk": 6, "pgarch": 3, "parch": 5, "tparch": 3, "sdgarch": 3,
              "saarch": 3, "tarch": 3, "earch": 3, "egarch": 3})


def _numlist(v) -> list[int]:
    from .tsmodels import _lags
    return _lags(v) if isinstance(v, str) else []


def _archm_fn(expr: str):
    safe = {"sqrt": np.sqrt, "ln": np.log, "log": np.log, "exp": np.exp, "abs": np.abs}

    def f(x):
        return float(eval(expr.replace("X", "x"), {"__builtins__": {}}, {**safe, "x": x}))   # noqa: S307
    return f


def cmd_arch(s: "Session", args: str) -> None:
    from .postest import replay
    from .regress import Sample, cmdline
    from .tsmodels import _check_gaps, _hr_start, _numderiv, _right, _g, _sample_text, _scores
    if not args.strip() or args.strip().startswith(","):
        replay(s, "arch", args)
        return
    c = find_top(args, ",")
    body, opts = (args, "") if c == -1 else (args[:c], args[c + 1:])
    o = match_options(opts, _OPTS) if opts.strip() else {}
    for bad in ("constraints", "from", "savespace"):
        if o.get(bad) and bad != "savespace":
            raise StataError(198, f"option {bad}() not implemented in OpenDTA")   # OpenDTA
    spec = ArchSpec()
    spec.ar, spec.ma = _numlist(o.get("ar")), _numlist(o.get("ma"))
    p = parse_standard(body)
    toks = p.varlist.split()
    if isinstance(o.get("arima"), str):
        nums = [int(float(v)) for v in o["arima"].replace(",", " ").split()]
        if len(nums) != 3:
            raise StataError(198, "arima() invalid")   # VERIFY
        pp, d, q = nums
        spec.ar, spec.ma = list(range(1, pp + 1)), list(range(1, q + 1))
        if d:
            toks = [("D." if d == 1 else f"D{d}.") + toks[0]] + toks[1:]
    for kd in TERMS:
        if o.get(kd):
            lags = _numlist(o[kd])
            if not lags:
                raise StataError(198, f"option {kd}() requires a numlist")
            spec.terms += [(kd, lg) for lg in lags]
    spec.terms.sort(key=lambda t: (TERMS.index(t[0]), t[1]))
    if o.get("archm"):
        spec.archm = [0]
    spec.archm += [lg for lg in _numlist(o.get("archmlags")) if lg not in spec.archm]
    spec.archm.sort()
    if isinstance(o.get("archmexp"), str):
        spec.archm_fn = _archm_fn(o["archmexp"])
        if not spec.archm:
            spec.archm = [0]
    if isinstance(o.get("distribution"), str):
        w = o["distribution"].replace(",", " ").split()
        name = w[0].lower()
        if name in ("gaussian", "normal"):
            spec.dist = "gaussian"
        elif name == "t":
            spec.dist = "t"
        elif name == "ged":
            spec.dist = "ged"
        else:
            raise StataError(198, "option distribution() incorrectly specified")
        if len(w) > 1:
            spec.dist_fixed = float(w[1])
            if spec.dist == "t" and spec.dist_fixed <= 2:
                raise StataError(198, "degrees of freedom must be greater than 2")
    if isinstance(o.get("arch0"), str):
        a0 = o["arch0"].strip().lower()
        if a0 in ("xb", "xb0", "xbwt", "xb0wt", "zero"):
            spec.arch0 = a0
        else:
            spec.arch0, spec.arch0_value = "value", float(a0)
    smp_args = " ".join(toks)
    if p.if_:
        smp_args += f" if {p.if_}"
    if p.in_:
        smp_args += f" in {p.in_}"
    extra_o = " ".join(f"{k}({v})" if isinstance(v, str) else k for k, v in o.items()
                       if k in ("noconstant", "level"))
    het_vars = s.expand_varlist(o["het"]) if isinstance(o.get("het"), str) else []
    smp = Sample(s, smp_args + (f", {extra_o}" if extra_o else ""), opts_spec={"noconstant": 4}, weights=(),
                 extra_vars=" ".join(het_vars))
    t0, t1 = _check_gaps(s, smp.mask)
    fvars.drop_collinear(smp.cols, smp.constant)
    est_cols = [col for col in smp.cols if col.values is not None]
    n = smp.n
    k = len(est_cols) + (1 if smp.constant else 0)
    info = require_ts(s.data)
    order = np.argsort(s.data.get(info.tvar).data[smp.mask], kind="stable")
    X = np.column_stack([col.values for col in est_cols] + ([np.ones(n)] if smp.constant else [])) \
        if k else np.zeros((n, 0))
    X, y = X[order], smp.y.astype(float)[order]
    Z = None
    if het_vars:
        hmask = smp.mask
        cols = fvars.build_design(s.data, fvars.expand_fv(s.data, " ".join(het_vars)), hmask, constant=True)
        Z = np.column_stack([c_.values for c_ in cols if c_.values is not None])[order]
        spec.nhet = Z.shape[1]
    spec.condobs = int(float(o["condobs"])) if isinstance(o.get("condobs"), str) else 0
    pos = spec.layout(k)
    npar = pos["n"]
    # starting values: OLS for the mean; small ARCH terms that keep the variance positive
    beta0 = np.linalg.lstsq(X, y, rcond=None)[0] if k else np.zeros(0)
    res0 = y - X @ beta0
    var0 = float(np.var(res0)) or 1.0
    b0 = np.zeros(npar)
    b0[pos["beta"][0]:pos["beta"][1]] = beta0
    if spec.ar or spec.ma:
        a0_, _sg = _hr_start(res0, spec.ar, spec.ma)
        b0[pos["ar"][0]:pos["ma"][1]] = a0_
    nt = len(spec.terms)
    ta = 0.0
    for (kd, lag), (a, bq) in zip(spec.terms, pos["terms"]):
        if kd in ("garch", "sdgarch", "pgarch"):
            b0[a] = 0.5 / max(1, sum(1 for t_ in spec.terms if t_[0] == kd))
        elif kd == "egarch":
            b0[a] = 0.5
        elif kd == "earch":
            b0[a], b0[a + 1] = 0.1, 0.1
        else:
            b0[a] = 0.1 / max(1, sum(1 for t_ in spec.terms if t_[0] == kd))
        ta += b0[a] if kd in ("arch", "garch") else 0.0
    if spec.nhet:
        b0[pos["het"][0]] = math.log(var0)
        if spec.nhet > 1:
            b0[pos["het"][0]:pos["het"][1]] = 0.0
            b0[pos["het"][1] - 1] = math.log(var0)
    else:
        if spec.egarch:
            b0[pos["gamma0"][0]] = math.log(var0) * 0.5
        elif spec.power:
            b0[pos["gamma0"][0]] = var0 ** 1.0 * 0.5
        else:
            b0[pos["gamma0"][0]] = var0 * max(0.05, 1 - ta)
    if spec.power:
        b0[pos["power"][0]] = 2.0 if any(t_[0] in POWER for t_ in spec.terms) else 1.0
    if pos["dist"][1] > pos["dist"][0]:
        b0[pos["dist"][0]] = math.log(6.0) if spec.dist == "t" else 0.0 if spec.dist == "ged" else 0.0
        if spec.dist == "t":
            b0[pos["dist"][0]] = math.log(4.0)
    xb0 = None

    def contrib(b):
        return archloglik(spec, b, y, X, Z, pos, xb0)

    def ll(b):
        c_ = contrib(b)
        return -np.inf if c_ is None else float(c_.sum())

    if not np.isfinite(ll(b0)):
        b0[pos["terms"][0][0]:pos["gamma0"][0]] = 0.0 if spec.terms else 0.0
        if not np.isfinite(ll(b0)):
            raise StataError(430, "initial values not feasible")
    if spec.arch0 in ("xb0", "xb0wt"):
        xb0 = _prime(spec, res0, 0.0)
    r = maximize_switching(ll, lambda bb: _scores(contrib, bb), b0, **ml_options(o))
    _l, _g_, Hm = _numderiv(ll, r.b)
    out = s.output
    if not o.get("nolog"):
        out.write("\n", "text")
        print_log(out, r, "log likelihood")
    b = r.b
    S = _scores(contrib, b)
    try:
        Hinv = np.linalg.inv(-Hm)
    except np.linalg.LinAlgError:
        Hinv = np.linalg.pinv(-Hm)
    vce = str(o.get("vce") or ("robust" if o.get("robust") else "opg")).strip().lower()
    vce = "robust" if vce.startswith("r") else ("oim" if vce == "oim" else "opg")
    if vce == "opg":
        V = np.linalg.inv(S.T @ S)
    elif vce == "oim":
        V = Hinv
    else:
        V = Hinv @ (S.T @ S) @ Hinv
    V = (V + V.T) / 2
    se = np.sqrt(np.maximum(np.diag(V), 0))
    # Wald test of the mean equation (without the constant)
    idx = [i for i in range(len(est_cols))] + list(range(pos["psi"][0], pos["ma"][1]))
    try:
        chi2 = float(b[idx] @ np.linalg.solve(V[np.ix_(idx, idx)], b[idx])) if idx else SYS
    except np.linalg.LinAlgError:
        chi2 = SYS
    pv = float(chi2_tail(chi2, len(idx))) if idx and chi2 < SYS else SYS
    # coefficient table
    eqname = smp.depname.split(".")[-1]
    names, eqs = [], []
    rows: list[CoefRow] = []
    bmap = {col.name: b[i] for i, col in enumerate(est_cols)}
    semap = {col.name: se[i] for i, col in enumerate(est_cols)}
    if smp.constant:
        bmap["_cons"], semap["_cons"] = b[len(est_cols)], se[len(est_cols)]
    if k:
        rows.append(CoefRow("eq", eqname))
        rows += rows_from_columns([c_ for c_ in smp.cols if c_.values is not None], bmap, semap, s.data,
                                  constant=smp.constant)
        names += [c_.name for c_ in est_cols] + (["_cons"] if smp.constant else [])
        eqs += [eqname] * k
    def add_group(title, labels, lo, hi, nm=None):
        nonlocal rows
        if rows:
            rows.append(CoefRow("sep"))
        rows.append(CoefRow("eq", title))
        for j, lab in zip(range(lo, hi), labels):
            rows.append(CoefRow("coef", lab, b[j], se[j], name=f"{title}:{lab}"))
            names.append(lab)
            eqs.append(title)
    if spec.archm:
        add_group("ARCHM", ["sigma2" if lg == 0 else f"L{lg}.sigma2" if lg > 1 else "L.sigma2" for lg in spec.archm],
                  *pos["psi"])
    if spec.ar or spec.ma:
        if rows:
            rows.append(CoefRow("sep"))
        rows.append(CoefRow("eq", "ARMA"))
        j = pos["ar"][0]
        for lab, lags in (("ar", spec.ar), ("ma", spec.ma)):
            if not lags:
                continue
            if rows[-1].kind != "eq":
                rows.append(CoefRow("blank"))
            rows.append(CoefRow("header", lab))
            for lg in lags:
                rows.append(CoefRow("coef", f"L{lg}." if lg > 1 else "L.", b[j], se[j], name=f"ARMA:L{lg}.{lab}"))
                names.append(f"L{lg}.{lab}" if lg > 1 else f"L.{lab}")
                eqs.append("ARMA")
                j += 1
    if spec.nhet:
        hnames = [c_.name for c_ in cols if c_.values is not None]
        hn = [h_ for h_ in hnames]
        add_group("HET", [("_cons" if x_ == "_cons" else x_) for x_ in hn], *pos["het"])
    if spec.terms or not spec.nhet:
        if rows:
            rows.append(CoefRow("sep"))
        rows.append(CoefRow("eq", "ARCH"))
        for kd in dict.fromkeys(kd for kd, _ in spec.terms):
            if rows[-1].kind != "eq":
                rows.append(CoefRow("blank"))
            rows.append(CoefRow("header", kd))
            for (kd2, lg), (a, bq) in zip(spec.terms, pos["terms"]):
                if kd2 != kd:
                    continue
                lab = f"L{lg}." if lg > 1 else "L."
                rows.append(CoefRow("coef", lab, b[a], se[a], name=f"ARCH:{lab}{kd}"))
                names.append(f"{lab}{kd}"); eqs.append("ARCH")
            if kd in EXTRA:
                rows.append(CoefRow("blank"))
                rows.append(CoefRow("header", EXTRA[kd]))
                for (kd2, lg), (a, bq) in zip(spec.terms, pos["terms"]):
                    if kd2 != kd:
                        continue
                    lab = f"L{lg}." if lg > 1 else "L."
                    rows.append(CoefRow("coef", lab, b[a + 1], se[a + 1], name=f"ARCH:{lab}{EXTRA[kd]}"))
                    names.append(f"{lab}{EXTRA[kd]}"); eqs.append("ARCH")
        if not spec.nhet:
            if spec.terms:
                rows.append(CoefRow("blank"))
            g_ = pos["gamma0"][0]
            rows.append(CoefRow("coef", "_cons", b[g_], se[g_], name="ARCH:_cons"))
            names.append("_cons"); eqs.append("ARCH")
    if spec.power:
        add_group("POWER", ["power"], *pos["power"])
    tdf = shape = None
    if pos["dist"][1] > pos["dist"][0]:
        j = pos["dist"][0]
        rows.append(CoefRow("sep"))
        lab = "/lndfm" if spec.dist == "t" else "/lnshape"
        rows.append(CoefRow("coef", lab, b[j], se[j], name=f"{lab[1:]}:_cons"))
        names.append("_cons"); eqs.append(lab[1:])
        val = math.exp(b[j]) + 2 if spec.dist == "t" else math.exp(b[j])
        sev = math.exp(b[j]) * se[j]
        rows.append(CoefRow("aux", "df" if spec.dist == "t" else "shape", val, sev))
        if spec.dist == "t":
            tdf = val
        else:
            shape = val
    elif spec.dist_fixed is not None:
        if spec.dist == "t":
            tdf = spec.dist_fixed
        else:
            shape = spec.dist_fixed
    kk = b.size
    est = Estimates("arch", smp.depname, names, b, V, n, "z", None)
    est.eqnames = eqs
    est.rows = rows
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level = float(o["level"]) if o.get("level") else smp.level
    est.vcetype = {"opg": "OPG", "robust": "Robust", "oim": "OIM"}[vce]
    est.extra = {"N": float(n), "ll": r.ll, "chi2": chi2, "p": pv, "df_m": len(idx), "sample": _sample_text(s, t0, t1),
                 "dist": {"gaussian": "Gaussian", "t": "Student's t", "ged": "GED"}[spec.dist], "ar": bool(spec.ar),
                 "ma": bool(spec.ma)}
    est.arch = {"spec": spec, "pos": pos, "y": y, "X": X, "Z": Z, "order": order, "mask": smp.mask.copy(), "k": k}
    est.display = display_arch
    est.predict = _predict_arch
    est.predict_opts = {"variance": 3, "residuals": 3, "xb": 2, "y": 1}
    scalars = [("N", float(n)), ("N_gaps", 0.0), ("condobs", float(spec.condobs)), ("k", float(kk)),
               ("k_eq", float(len(dict.fromkeys(eqs)))), ("k_dv", 1.0), ("df_m", float(len(idx))),
               ("ll", r.ll), ("chi2", chi2), ("p", pv), ("archi", float(1)), ("archany", float(bool(spec.terms))),
               ("tmin", t0), ("tmax", t1), ("rank", float(kk)), ("ic", float(r.iterations)),
               ("rc", 0.0), ("converged", 1.0 if r.converged else 0.0)]
    if tdf is not None:
        scalars.append(("tdf", tdf))
    if shape is not None:
        scalars.append(("shape", shape))
    macros = [("cmd", "arch"), ("cmdline", cmdline("arch", args)), ("depvar", smp.depname),
              ("eqnames", " ".join(dict.fromkeys(eqs))), ("title", "ARCH family regression"),
              ("dist", {"gaussian": "gaussian", "t": "t", "ged": "ged"}[spec.dist]),
              ("vce", vce), ("vcetype", est.vcetype), ("chi2type", "Wald"), ("tvar", info.tvar),
              ("ar", " ".join(map(str, spec.ar))), ("ma", " ".join(map(str, spec.ma))),
              ("predict", "arch_p"), ("properties", "b V")]
    for kd in dict.fromkeys(kd for kd, _ in spec.terms):
        macros.append((kd, " ".join(str(lg) for k2, lg in spec.terms if k2 == kd)))
    if spec.archm:
        macros.append(("archm", "archm"))
    post(s, est, scalars, macros, smp.mask)
    display_arch(s, est)


def display_arch(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                 level: float | None = None, **_kw) -> None:
    from .tsmodels import _g, _right
    out = s.output
    x = est.extra
    if header:
        suffix = " -- ARMA disturbances" if x["ar"] and x["ma"] else " -- AR disturbances" if x["ar"] \
            else " -- MA disturbances" if x["ma"] else ""
        out.write(f"\nARCH family regression{suffix}\n\n", "text")
        _right(out, x["sample"], "Number of obs", f"{int(x['N']):,}")
        _right(out, f"Distribution: {x['dist']}", f"Wald chi2({int(x['df_m'])})", _g(x["chi2"], "%10.2f"))
        _right(out, f"Log likelihood = {_g(x['ll'], '%9.0g')}", "Prob > chi2", _g(x["p"], "%10.4f"))
        out.write("\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=level or est.level, vcetype=est.vcetype)


def _predict_arch(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    """Over the estimation sample: xb (mean without ARCH-M and ARMA terms), y (conditional mean),
    residuals (epsilon) and variance (conditional variance)."""
    A = est.arch
    pth = archpath(A["spec"], est.b, A["y"], A["X"], A["Z"], A["pos"])
    if pth is None:
        raise StataError(498, "prediction failed")   # VERIFY
    rows = np.nonzero(A["mask"])[0][A["order"]]
    out = np.full(s.data.nobs, SYS)
    beta = est.b[:A["k"]]
    if stat in ("", "xb"):
        out[rows] = A["X"] @ beta
        return out, "" if stat else "(option xb assumed; fitted values)"
    if stat == "y":
        out[rows] = pth["mu"]
    elif stat == "variance":
        out[rows] = pth["s2"]
    elif stat == "residuals":
        out[rows] = pth["e"]
    else:
        raise StataError(198, f"option {stat} not allowed")
    return out, ""

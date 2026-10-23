"""bayesmh and the bayes prefix ([BAYES] bayesmh, [BAYES] bayes), with bayesstats, bayesgraph and
bayestest.

A model is a set of equations, each with a mean expression and a likelihood, and priors on the
parameters {name}, {eq:varlist} (a linear form with a constant) and {eq:}. The posterior is
sampled with the adaptive random-walk Metropolis-Hastings of core/mcmc.py in blocks (the
coefficients of a linear form together, every other parameter alone, or the blocks of block()).
The random numbers come from numpy, so the draws of Stata cannot be reproduced; the summaries
(mean, standard deviation, MCSE, median, credible intervals, efficiency) agree with Stata's up
to Monte Carlo error. The layouts of the header, the model summary and the tables are VERIFY.

Likelihoods: normal(var), lognormal(var), t(var, df), exponential, poisson, logit, probit and
binomial(n) (the last four with the equation's expression as the linear predictor).
Priors: normal(m, v), lognormal(m, v), uniform(a, b), igamma(a, b), gamma(a, b), beta(a, b),
exponential(m), chi2(df), flat and jeffreys (improper: the log marginal likelihood is then
missing)."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from ..core import mcmc
from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options
from ..lazy import lazy
from .results import Estimates, post

special = lazy("scipy.special")
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
_BRACE = re.compile(r"\{([^{}]*)\}")
_LIKELIHOODS = ("normal", "lognormal", "t", "exponential", "poisson", "logit", "probit", "binomial")
_PRIORS = ("normal", "lognormal", "uniform", "igamma", "gamma", "beta", "exponential", "chi2", "flat",
           "jeffreys")
_POSITIVE = ("lognormal", "igamma", "gamma", "beta", "exponential", "chi2", "jeffreys")
_TITLE = {"normal": "normal regression", "lognormal": "lognormal regression", "t": "t regression",
          "exponential": "exponential regression", "poisson": "Poisson regression",
          "logit": "logistic regression", "probit": "probit regression",
          "binomial": "binomial regression"}
_OPTS = {"likelihood": 4, "prior": 5, "block": 5, "initial": 4, "burnin": 6, "mcmcsize": 7,
         "thinning": 5, "rseed": 5, "nomodelsummary": 12, "dots": 4, "saving": 6, "replace": 7,
         "hpd": 3, "clevel": 4, "notable": 5, "noheader": 6, "nomleinitial": 5, "scale": 3,
         "adaptation": 5, "exclude": 4, "noshow": 6, "level": 1, "noisily": 3, "sampling": 4,
         "noconstant": 6, "nodots": 5}


# ---------------------------------------------------------------------------
# parameters and expressions
# ---------------------------------------------------------------------------

@dataclass
class Params:
    names: list[str] = field(default_factory=list)       # "sigma2", "y:x1"
    eqs: list[str] = field(default_factory=list)         # "" for scalars
    init: list[float | None] = field(default_factory=list)
    lcs: dict[str, list[int]] = field(default_factory=dict)       # label -> parameter indices
    lc_vars: dict[str, list[str]] = field(default_factory=dict)   # label -> variables ("" = constant)

    def add(self, name: str, eq: str = "", init=None) -> int:
        if name in self.names:
            i = self.names.index(name)
            if init is not None:
                self.init[i] = init
            return i
        self.names.append(name)
        self.eqs.append(eq)
        self.init.append(init)
        return len(self.names) - 1

    def add_lc(self, label: str, variables: list[str], const: bool) -> list[int]:
        if label in self.lcs:
            return self.lcs[label]
        vs = list(variables) + ([""] if const else [])
        idx = [self.add(f"{label}:{v or '_cons'}", label) for v in vs]
        self.lcs[label] = idx
        self.lc_vars[label] = vs
        return idx


def _parse_brace(body: str):
    """-> ("scalar", name, init) | ("lc", label, vars, noconstant)."""
    body = body.strip()
    if ":" in body:
        label, rest = body.split(":", 1)
        noconst = False
        if "," in rest:
            rest, opt = rest.split(",", 1)
            noconst = opt.strip().startswith("nocons")
        return ("lc", label.strip(), rest.split(), noconst)
    init = None
    name = body
    if "=" in body:
        name, val = body.split("=", 1)
        try:
            init = float(val)
        except ValueError:
            raise StataError(198, f"invalid initial value in {{{body}}}")   # VERIFY
    name = name.strip()
    if not re.fullmatch(r"[A-Za-z_]\w*", name):
        raise StataError(198, f"invalid parameter name {{{body}}}")   # VERIFY
    return ("scalar", name, init)


class Compiled:
    """An expression with parameters {...} compiled once: parameters become scalars of the
    session and the parsed expression is evaluated by the vectorized evaluator."""

    def __init__(self, s: "Session", text: str, P: Params, tag: str):
        from ..lang.expr import parse
        self.s, self.text = s, text
        self.map: list[tuple[str, int]] = []         # (scalar name, parameter index)
        counter = [0]

        def sub(m):
            tk = _parse_brace(m.group(1))
            if tk[0] == "scalar":
                i = P.add(tk[1], "", tk[2])
                nm = f"__bz_{tag}_{counter[0]}"
                counter[0] += 1
                self.map.append((nm, i))
                return nm
            _k, label, vs, noconst = tk
            idx = P.add_lc(label, vs, not noconst)
            terms = []
            for i, v in zip(idx, P.lc_vars[label]):
                nm = f"__bz_{tag}_{counter[0]}"
                counter[0] += 1
                self.map.append((nm, i))
                terms.append(f"{nm}*{v}" if v else nm)
            return "(" + " + ".join(terms) + ")"
        self.rendered = _BRACE.sub(sub, text)
        self.node = parse(self.rendered) if self.rendered.strip() else None
        self.const = not self.map

    def set(self, theta: np.ndarray) -> None:
        sc = self.s.scalars
        for nm, i in self.map:
            sc[nm] = float(theta[i])

    def value(self, theta: np.ndarray):
        from ..lang.vexpr import VectorContext
        from ..commands._util import groups
        from ..lang.vexpr import evaluate_vec
        self.set(theta)
        return evaluate_vec(self.node, VectorContext(self.s, groups(self.s)))

    def scalar(self, theta: np.ndarray) -> float:
        v = self.value(theta)
        return float(np.asarray(v, float).ravel()[0])

    def clean(self) -> None:
        for nm, _ in self.map:
            self.s.scalars.pop(nm, None)


class Arg:
    """A numeric argument of a likelihood or prior: a constant or an expression of parameters."""

    def __init__(self, s, text: str, P: Params, tag: str):
        self.text = text.strip()
        try:
            self.val = float(self.text)
            self.c = None
        except ValueError:
            self.c = Compiled(s, self.text, P, tag)
            self.val = None
            if self.c.const:
                self.val = self.c.scalar(np.zeros(len(P.names)))
                self.c = None

    def __call__(self, theta) -> float:
        return self.val if self.c is None else self.c.scalar(theta)


# ---------------------------------------------------------------------------
# densities
# ---------------------------------------------------------------------------

_LN2PI = math.log(2 * math.pi)


def _lik_logpdf(dist: str, y: np.ndarray, mu: np.ndarray, a: list[float], extra=None) -> float:
    if dist == "normal":
        v = a[0]
        if not v > 0:
            return -np.inf
        return float(-0.5 * len(y) * (_LN2PI + math.log(v)) - np.sum((y - mu) ** 2) / (2 * v))
    if dist == "lognormal":
        v = a[0]
        if not v > 0 or np.any(y <= 0):
            return -np.inf
        z = np.log(y)
        return float(-np.sum(z) - 0.5 * len(y) * (_LN2PI + math.log(v)) - np.sum((z - mu) ** 2) / (2 * v))
    if dist == "t":
        v, df = a[0], a[1]
        if not (v > 0 and df > 0):
            return -np.inf
        r2 = (y - mu) ** 2 / (df * v)
        return float(len(y) * (special.gammaln((df + 1) / 2) - special.gammaln(df / 2)
                               - 0.5 * math.log(math.pi * df * v)) - (df + 1) / 2 * np.sum(np.log1p(r2)))
    if dist == "exponential":
        return float(np.sum(-mu - y * np.exp(-mu)))
    if dist == "poisson":
        return float(np.sum(y * mu - np.exp(np.clip(mu, -700, 700)) - special.gammaln(y + 1)))
    if dist == "logit":
        return float(np.sum(y * mu - np.logaddexp(0, mu)))
    if dist == "probit":
        return float(np.sum(y * special.log_ndtr(mu) + (1 - y) * special.log_ndtr(-mu)))
    if dist == "binomial":
        n = extra
        return float(np.sum(special.gammaln(n + 1) - special.gammaln(y + 1) - special.gammaln(n - y + 1)
                            + y * mu - n * np.logaddexp(0, mu)))
    raise StataError(198, f"likelihood({dist}) is not available")   # VERIFY


def _prior_logpdf(dist: str, x: float, a: list[float]) -> float:
    if dist == "normal":
        m, v = a
        return -0.5 * (_LN2PI + math.log(v)) - (x - m) ** 2 / (2 * v) if v > 0 else -np.inf
    if dist == "lognormal":
        m, v = a
        if x <= 0 or v <= 0:
            return -np.inf
        return -math.log(x) - 0.5 * (_LN2PI + math.log(v)) - (math.log(x) - m) ** 2 / (2 * v)
    if dist == "uniform":
        lo, hi = a
        return -math.log(hi - lo) if lo <= x <= hi and hi > lo else -np.inf
    if dist == "igamma":
        al, be = a
        if x <= 0 or al <= 0 or be <= 0:
            return -np.inf
        return al * math.log(be) - math.lgamma(al) - (al + 1) * math.log(x) - be / x
    if dist == "gamma":
        al, be = a
        if x <= 0 or al <= 0 or be <= 0:
            return -np.inf
        return -math.lgamma(al) - al * math.log(be) + (al - 1) * math.log(x) - x / be
    if dist == "beta":
        al, be = a
        if not (0 < x < 1) or al <= 0 or be <= 0:
            return -np.inf
        return (al - 1) * math.log(x) + (be - 1) * math.log1p(-x) - float(special.betaln(al, be))
    if dist == "exponential":
        m = a[0]
        return -math.log(m) - x / m if x >= 0 and m > 0 else -np.inf
    if dist == "chi2":
        k = a[0] / 2
        return -math.lgamma(k) - k * math.log(2) + (k - 1) * math.log(x) - x / 2 if x > 0 else -np.inf
    if dist == "flat":
        return 0.0
    if dist == "jeffreys":
        return -math.log(x) if x > 0 else -np.inf
    raise StataError(198, f"prior {dist}() is not available")   # VERIFY


# ---------------------------------------------------------------------------
# the model
# ---------------------------------------------------------------------------

@dataclass
class Equation:
    depvar: str
    text: str                     # expression with parameters
    lik: str
    lik_args: list                # Arg
    lik_args_text: list[str]
    mean: Compiled | None = None
    lc: str | None = None         # label when the mean is a single linear form
    rows: np.ndarray | None = None
    y: np.ndarray | None = None
    X: np.ndarray | None = None
    nbin: np.ndarray | None = None


def _split_args(text: str) -> list[str]:
    out, depth, cur = [], 0, ""
    for ch in text:
        if ch in "({":
            depth += 1
        elif ch in ")}":
            depth -= 1
        if ch == "," and depth == 0:
            out.append(cur)
            cur = ""
        else:
            cur += ch
    out.append(cur)
    return [a.strip() for a in out if a.strip() != "" or len(out) > 1]


def _calls(text: str) -> list[tuple[str, str]]:
    """'prior({a}, normal(0,1))' style: returns [(name, inner)] of the top-level name(...) items
    of an option value."""
    return [(m.group(1), m.group(2)) for m in re.finditer(r"(\w+)\(((?:[^()]|\([^()]*\))*)\)", text)]


def _option_values(opts: str, name: str) -> list[str]:
    """All values of a repeatable option name(...) in the option text."""
    out, i = [], 0
    pat = re.compile(rf"\b{name}\(")
    while True:
        m = pat.search(opts, i)
        if not m:
            return out
        depth, j = 1, m.end()
        while j < len(opts) and depth:
            depth += {"(": 1, ")": -1}.get(opts[j], 0)
            j += 1
        out.append(opts[m.end():j - 1])
        i = j


def _strip_option(opts: str, name: str) -> str:
    out, i = "", 0
    pat = re.compile(rf"\b{name}\(")
    while True:
        m = pat.search(opts, i)
        if not m:
            return out + opts[i:]
        depth, j = 1, m.end()
        while j < len(opts) and depth:
            depth += {"(": 1, ")": -1}.get(opts[j], 0)
            j += 1
        out += opts[i:m.start()] + " "
        i = j


class Model:
    def __init__(self, s: "Session", spec: str, o_text: str, mask):
        self.s = s
        self.P = Params()
        self.eqs: list[Equation] = []
        self.priors: list[tuple[int, str, list[Arg], str]] = []     # (param, dist, args, text)
        self.prior_text: list[tuple[list[int], str]] = []
        self.improper = False
        liks = _option_values(o_text, "likelihood")
        self._parse_equations(spec, liks)
        self._parse_priors(_option_values(o_text, "prior"))
        self._data(mask)

    # -- parsing ---------------------------------------------------------
    def _parse_equations(self, spec: str, liks: list[str]) -> None:
        s = self.s
        spec = spec.strip()
        items: list[tuple[str, str | None]] = []
        if spec.startswith("("):
            depth, start = 0, 0
            for i, ch in enumerate(spec):
                if ch == "(":
                    if depth == 0:
                        start = i + 1
                    depth += 1
                elif ch == ")":
                    depth -= 1
                    if depth == 0:
                        body = spec[start:i]
                        own = _option_values(body, "likelihood")
                        body = _strip_option(body, "likelihood")
                        c = find_top(body, ",")
                        head, tail = (body, "") if c == -1 else (body[:c], body[c + 1:])
                        items.append((head.strip() + ("\x00noconstant" if "noconstant" in tail else ""),
                                      own[0] if own else None))
        else:
            items.append((spec, None))
        if not items:
            raise StataError(198, "no model specified")   # VERIFY
        for k, (head, own) in enumerate(items):
            noconst = "\x00noconstant" in head
            head = head.replace("\x00noconstant", "")
            lik_text = own or (liks[k] if k < len(liks) else (liks[0] if liks else None))
            if lik_text is None:
                raise StataError(198, "option likelihood() required")   # VERIFY
            m = re.match(r"\s*(\w+)\s*\((.*)\)\s*$", lik_text, re.S)
            name = lik_text.strip() if not m else m.group(1)
            inner = "" if not m else m.group(2)
            hit = [d for d in _LIKELIHOODS if d == name] or [d for d in _LIKELIHOODS if d.startswith(name)]
            if len(hit) != 1:
                raise StataError(198, f"likelihood({lik_text.strip()}) is not available")   # VERIFY
            dist = hit[0]
            args_t = _split_args(inner) if inner.strip() else []
            if "=" in head and not head.lstrip().startswith("{"):
                lhs, rhs = head.split("=", 1)
                depvar, text = lhs.strip(), rhs.strip()
            else:
                words = head.split()
                depvar = words[0]
                rest = " ".join(words[1:])
                names = s.expand_varlist(depvar)
                depvar = names[0]
                vars_ = s.expand_varlist(rest) if rest.strip() else []
                body = f"{{{depvar}:{' '.join(vars_)}{', noconstant' if noconst else ''}}}"
                text = body
            if not s.data.has(depvar):
                raise StataError(111, f"variable {depvar} not found")
            depvar = s.data.get(depvar).name
            eq = Equation(depvar, text, dist, [], args_t)
            eq.mean = Compiled(s, text, self.P, f"e{k}")
            m2 = re.fullmatch(r"\{([^{}]*:[^{}]*)\}", text.strip())
            if m2:
                eq.lc = _parse_brace(m2.group(1))[1]
            for j, a in enumerate(args_t):
                if dist == "binomial" and j == 0:
                    continue
                eq.lik_args.append(Arg(s, a, self.P, f"l{k}_{j}"))
            if dist == "normal" and not args_t or dist == "lognormal" and not args_t:
                raise StataError(198, f"likelihood({dist}) requires the variance")   # VERIFY
            self.eqs.append(eq)

    def _targets(self, text: str) -> list[int]:
        P = self.P
        out: list[int] = []
        toks = _BRACE.findall(text)
        if not toks:
            raise StataError(198, f"invalid prior() parameter list: {text.strip()}")   # VERIFY
        for body in toks:
            if ":" in body:
                label, rest = body.split(":", 1)
                label = label.strip()
                if label not in P.lcs:
                    raise StataError(198, f"parameter {{{body}}} not found in the model")   # VERIFY
                if not rest.strip():
                    out.extend(P.lcs[label])
                else:
                    for v in rest.replace(",", " ").split():
                        nm = f"{label}:{v}"
                        if nm not in P.names:
                            raise StataError(198, f"parameter {{{nm}}} not found in the model")   # VERIFY
                        out.append(P.names.index(nm))
            else:
                for nm in body.split():
                    nm = nm.split("=")[0]
                    out.append(P.add(nm))
        return out

    def _parse_priors(self, texts: list[str]) -> None:
        for t in texts:
            c = find_top(t, ",")
            if c == -1:
                raise StataError(198, f"invalid prior({t})")   # VERIFY
            head, dist_t = t[:c], t[c + 1:].strip()
            targets = self._targets(head)
            m = re.match(r"(\w+)\s*(?:\((.*)\))?\s*$", dist_t, re.S)
            if not m:
                raise StataError(198, f"invalid prior distribution {dist_t}")   # VERIFY
            name = m.group(1)
            hit = [d for d in _PRIORS if d == name] or [d for d in _PRIORS if d.startswith(name)]
            if len(hit) != 1:
                raise StataError(198, f"prior {name}() is not available")   # VERIFY
            dist = hit[0]
            args = [Arg(self.s, a, self.P, f"p{len(self.priors)}_{j}")
                    for j, a in enumerate(_split_args(m.group(2) or ""))]
            need = {"normal": 2, "lognormal": 2, "uniform": 2, "igamma": 2, "gamma": 2, "beta": 2,
                    "exponential": 1, "chi2": 1, "flat": 0, "jeffreys": 0}[dist]
            if len(args) != need:
                raise StataError(198, f"prior {dist}() requires {need} argument(s)")   # VERIFY
            if dist in ("flat", "jeffreys"):
                self.improper = True
            for i in targets:
                self.priors.append((i, dist, args, dist_t))
            self.prior_text.append((targets, dist_t))
        covered = {i for i, *_ in self.priors}
        miss = [self.P.names[i] for i in range(len(self.P.names)) if i not in covered]
        if miss:
            raise StataError(198, "no prior specified for " + " ".join("{" + n + "}" for n in miss))   # VERIFY

    # -- data ------------------------------------------------------------
    def _data(self, mask) -> None:
        s, ds = self.s, self.s.data
        n = ds.nobs
        ok = mask.copy()
        theta0 = np.array([0.0 if v is None else v for v in self.P.init], float)
        for k, eq in enumerate(self.eqs):
            y = ds.get(eq.depvar).data
            ok &= y < SYS
            if eq.lc is not None:
                for v in self.P.lc_vars[eq.lc]:
                    if v:
                        ok &= ds.get(v).data < SYS
        self.rows = np.nonzero(ok)[0]
        if len(self.rows) == 0:
            raise StataError(2000, "no observations")
        self.N = len(self.rows)
        for eq in self.eqs:
            eq.y = ds.get(eq.depvar).data[self.rows].astype(float)
            if eq.lc is not None:
                cols = [ds.get(v).data[self.rows].astype(float) if v else np.ones(self.N)
                        for v in self.P.lc_vars[eq.lc]]
                eq.X = np.column_stack(cols)
            if eq.lik == "binomial":
                a = eq.lik_args_text[0]
                try:
                    eq.nbin = np.full(self.N, float(a))
                except ValueError:
                    eq.nbin = ds.get(a.strip()).data[self.rows].astype(float)

    # -- posterior -------------------------------------------------------
    def mean(self, eq: Equation, theta: np.ndarray) -> np.ndarray:
        if eq.X is not None:
            return eq.X @ theta[self.P.lcs[eq.lc]]
        v = np.asarray(eq.mean.value(theta), float)
        if v.ndim == 0:
            return np.full(self.N, float(v))
        return v[self.rows]

    def loglik(self, theta: np.ndarray) -> float:
        tot = 0.0
        for eq in self.eqs:
            mu = self.mean(eq, theta)
            if not np.all(np.isfinite(mu)):
                return -np.inf
            a = [g(theta) for g in eq.lik_args]
            v = _lik_logpdf(eq.lik, eq.y, mu, a, eq.nbin)
            if not np.isfinite(v):
                return -np.inf
            tot += v
        return tot

    def logprior(self, theta: np.ndarray) -> float:
        tot = 0.0
        for i, dist, args, _t in self.priors:
            v = _prior_logpdf(dist, float(theta[i]), [g(theta) for g in args])
            if not np.isfinite(v):
                return -np.inf
            tot += v
        return tot

    def logpost(self, theta: np.ndarray) -> float:
        lp = self.logprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        ll = self.loglik(theta)
        return ll + lp if np.isfinite(ll) else -np.inf

    def detach(self) -> None:
        """Drops everything that refers to the session (the estimation results are copied by
        estimates store)."""
        self.clean()
        self.s = None
        for eq in self.eqs:
            eq.mean = None
            eq.lik_args = []
        self.priors = []

    def clean(self) -> None:
        for eq in self.eqs:
            if eq.mean is not None:
                eq.mean.clean()
            for a in eq.lik_args:
                if a.c is not None:
                    a.c.clean()
        for _i, _d, args, _t in self.priors:
            for a in args:
                if a.c is not None:
                    a.c.clean()


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------

def _parse_command(args: str):
    """spec [if] [in], options"""
    c = find_top(args, ",")
    head, opts = (args, "") if c == -1 else (args[:c], args[c + 1:])
    m = re.search(r"(?:^|[\s)])((?:if|in)\b.*)$", head, re.S)
    spec, tail = (head[:m.start(1)], m.group(1)) if m else (head, "")
    return spec.strip(), tail.strip(), opts


def _simple_opts(opts: str) -> dict:
    cleaned = opts
    for name in ("likelihood", "prior", "block", "initial", "saving", "exclude", "adaptation", "scale",
                 "sampling"):
        cleaned = _strip_option(cleaned, name)
    return match_options(cleaned, _OPTS) if cleaned.strip() else {}


def cmd_bayesmh(s: "Session", args: str, *, title: str | None = None, lik_names: dict | None = None) -> None:
    if not args.strip() or args.strip().startswith(","):
        _replay(s, args)
        return
    from ..commands._util import touse
    from ..lang.syntax import parse_standard
    spec, tail, opts = _parse_command(args)
    o = _simple_opts(opts)
    p = parse_standard(tail)
    mask = touse(s, p)
    nburn = int(float(o["burnin"])) if o.get("burnin") else 2500
    nsamp = int(float(o["mcmcsize"])) if o.get("mcmcsize") else 10000
    thin = int(float(o["thinning"])) if o.get("thinning") else 1
    if nburn < 0 or nsamp < 2 or thin < 1:
        raise StataError(198, "invalid burnin(), mcmcsize() or thinning()")   # VERIFY
    seed = int(float(o["rseed"])) if o.get("rseed") else None
    rng = np.random.default_rng(seed)
    saved_scalars = dict(s.scalars)
    model = Model(s, spec, opts, mask)
    try:
        P = model.P
        d = len(P.names)
        theta0 = np.array([1.0 if (v is None and P.eqs[i] == "") else (0.0 if v is None else v)
                           for i, v in enumerate(P.init)], float)
        for item in _option_values(opts, "initial"):
            for body, val in re.findall(r"\{([^{}]*)\}\s*([-+.\w]+)", item):
                tk = _parse_brace(body)
                if tk[0] == "scalar" and tk[1] in P.names:
                    theta0[P.names.index(tk[1])] = float(val)
                elif tk[0] == "lc":
                    for v in tk[2] or [""]:
                        nm = f"{tk[1]}:{v or '_cons'}"
                        if nm in P.names:
                            theta0[P.names.index(nm)] = float(val)
        # blocks
        blocks: list[list[int]] = []
        used: set[int] = set()
        for text in _option_values(opts, "block"):
            c = find_top(text, ",")
            tg = model._targets(text if c == -1 else text[:c])
            blocks.append([i for i in tg if i not in used])
            used |= set(tg)
        for label, idx in P.lcs.items():
            rest = [i for i in idx if i not in used]
            if rest:
                blocks.append(rest)
                used |= set(rest)
        for i in range(d):
            if i not in used:
                blocks.append([i])
        blocks = [b for b in blocks if b]
        # starting point: posterior mode (positive parameters on the log scale)
        pos = {i for i, dist, _a, _t in model.priors if dist in _POSITIVE}

        def to_theta(u):
            t = u.copy()
            for i in pos:
                t[i] = math.exp(min(u[i], 700))
            return t
        u0 = theta0.copy()
        for i in pos:
            u0[i] = math.log(max(theta0[i], 1e-8))
        if not np.isfinite(model.logpost(theta0)):
            raise StataError(498, "initial values not feasible")   # VERIFY
        if o.get("nomleinitial"):
            x0 = theta0
            covs = None
        else:
            u, _ok = mcmc.find_mode(lambda u: model.logpost(to_theta(u)), u0)
            x0 = to_theta(u)
            if not np.isfinite(model.logpost(x0)):
                x0 = theta0
            covs = mcmc.proposal_cov(model.logpost, x0, blocks)
        out = s.output
        out.write("Burn-in ...\n", "text")
        out.write("Simulation ...\n", "text")
        res = mcmc.adaptive_mh(model.logpost, x0, blocks, nburn, nsamp, thin, rng, covs)
        draws = res["draws"]
        lp = res["lp"]
        # log likelihood at the draws for DIC
        ll_draws = np.array([model.loglik(draws[k]) for k in range(0, len(draws), max(len(draws) // 2000, 1))])
        dev_mean = float(-2 * ll_draws.mean())
        dev_at_mean = float(-2 * model.loglik(draws.mean(axis=0)))
        dic = dev_mean + (dev_mean - dev_at_mean)
        lml = float("nan") if model.improper else mcmc.laplace_metropolis(draws, lp)
        level = float(o["clevel"]) if o.get("clevel") else float(s.settings.get("level", 95))
        summ = mcmc.summarize(draws, level, bool(o.get("hpd")))
        eff = [r["ess"] / len(draws) for r in summ]
        info = dict(title=title or _model_title(model), model=model, lik_names=lik_names or {},
                    nomodelsummary=bool(o.get("nomodelsummary")), burn=nburn, size=nsamp, thin=thin,
                    iters=nburn + nsamp * thin, acc=res["acc_rate"], eff=eff, lml=lml, dic=dic,
                    hpd=bool(o.get("hpd")), level=level, draws=draws, lp=lp, N=model.N, seed=seed)
        saving = _option_values(opts, "saving")
        info["saving"] = bool(saving)
        if saving:
            _save_draws(s, saving[0], model, draws, lp, ll_draws)
        _post(s, args, model, draws, info, summ, mask)
    finally:
        model.detach()
        for k in list(s.scalars):
            if k.startswith("__bz_") and k not in saved_scalars:
                s.scalars.pop(k, None)
    display_bayes(s, s._estimates, header=not o.get("noheader"), table=not o.get("notable"),
                  summary=not o.get("nomodelsummary"))


def _replay(s: "Session", args: str) -> None:
    est = _current(s)
    o = match_options(args.lstrip(", "), {"hpd": 3, "clevel": 4, "noheader": 6, "notable": 5,
                                          "level": 1}) if args.strip(" ,") else {}
    lv = float(o["clevel"]) if o.get("clevel") else (float(o["level"]) if o.get("level") else None)
    display_bayes(s, est, header=not o.get("noheader"), table=not o.get("notable"), level=lv,
                  hpd=True if o.get("hpd") else None)


def _model_title(model: Model) -> str:
    if len(model.eqs) == 1:
        return "Bayesian " + _TITLE[model.eqs[0].lik]
    return "Bayesian regression"


def _save_draws(s, text: str, model: Model, draws, lp, ll_draws) -> None:
    from ..core.dataset import Variable
    c = find_top(text, ",")
    path = (text if c == -1 else text[:c]).strip().strip('"')
    replace = c != -1 and "replace" in text[c:]
    if not path.lower().endswith(".dta"):
        path += ".dta"
    n = draws.shape[0]
    s.run("preserve")
    try:
        s.run("quietly clear")
        s.run(f"quietly set obs {n}")
        ds = s.data
        ds.add(Variable("_index", "long", np.arange(1, n + 1, dtype=float)))
        ll_full = np.array([model.loglik(draws[k]) for k in range(n)])
        ds.add(Variable("_loglikelihood", "double", ll_full))
        ds.add(Variable("_logposterior", "double", lp))
        for j, nm in enumerate(model.P.names):
            vn = "eq" + str(j + 1) + "_" + re.sub(r"\W", "_", nm)
            ds.add(Variable(vn[:32], "double", draws[:, j]))
        s.run(f'quietly save "{path}"{", replace" if replace else ""}')
    finally:
        s.run("restore")


def _post(s, args, model: Model, draws, info, summ, mask) -> None:
    P = model.P
    names = [n.split(":", 1)[1] if ":" in n else n for n in P.names]
    eqn = [P.eqs[i] if P.eqs[i] else "_" for i in range(len(P.names))]
    mean = np.array([r["mean"] for r in summ])
    C = np.atleast_2d(np.cov(draws.T))
    est = Estimates("bayesmh", " ".join(e.depvar for e in model.eqs), names, mean.copy(), C, model.N, "z", None)
    est.eqnames = eqn
    est.level = info["level"]
    est.model = model
    est.extra = dict(info)
    est.extra["summ"] = summ
    est.extra["pnames"] = list(P.names)
    est.display = display_bayes
    est.predict = None
    est.estat = {}
    effv = info["eff"]
    sample = np.zeros(s.data.nobs, dtype=bool)
    sample[model.rows] = True
    scalars = [("N", float(model.N)), ("mcmcsize", float(info["size"])), ("burnin", float(info["burn"])),
               ("mcmciter", float(info["iters"])), ("thinning", float(info["thin"])),
               ("arate", info["acc"]), ("corrlag", 0.0), ("eff_min", float(min(effv))),
               ("eff_avg", float(np.mean(effv))), ("eff_max", float(max(effv))),
               ("lml_lm", info["lml"]), ("dic", info["dic"]), ("clevel", info["level"])]
    macros = [("cmdline", f"bayesmh {args.strip()}"), ("cmd", "bayesmh"), ("title", info["title"]),
              ("method", "mh"), ("sampler", "adaptive"), ("depvars", " ".join(e.depvar for e in model.eqs)),
              ("properties", "b V")]
    from ..commands.matrix import Matrix
    cn = list(P.names)

    def row(v):
        return Matrix(np.asarray(v, float).reshape(1, -1), ["r1"], cn)
    mats = [("mean", row(mean)), ("median", row([r["median"] for r in summ])),
            ("sd", row([r["sd"] for r in summ])), ("mcse", row([r["mcse"] for r in summ])),
            ("Cov", Matrix(C, cn, cn))]
    post(s, est, scalars, macros, sample, extra_mats=mats)


def _g(x, fmt="%9.0g"):
    from ..core.formats import format_value
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def _g4(x: float) -> str:
    """Four significant digits without the leading zero (.2889, .05538)."""
    if x is None or not np.isfinite(x):
        return "."
    t = f"{x:.4g}"
    if "e" in t:
        return t
    return t.replace("0.", ".", 1) if t.startswith("0.") else t.replace("-0.", "-.", 1)


def _commas(n: int) -> str:
    return f"{n:,}"


def _model_summary(s, model: Model, info) -> None:
    out = s.output
    line = "-" * 78
    out.write("\nModel summary\n" + line + "\n", "text")
    out.write("Likelihood: \n", "text")
    for eq in model.eqs:
        args = ",".join([f"xb_{eq.lc}" if eq.lc else eq.text] + [a for a in eq.lik_args_text])
        name = info["lik_names"].get(eq.depvar, eq.lik)
        out.write(f"  {eq.depvar} ~ {name}({args})\n", "text")
    out.write("\nPriors: \n", "text")
    notes: list[str] = []
    lines = []
    for targets, dist_t in model.prior_text:
        P = model.P
        labs = {P.eqs[i] for i in targets if P.eqs[i]}
        if targets and all(P.eqs[i] for i in targets) and len(labs) == 1 and \
                set(targets) == set(P.lcs.get(next(iter(labs)), [])):
            lab = "{" + f"{next(iter(labs))}:{' '.join(P.names[i].split(':', 1)[1] for i in targets)}" + "}"
            mark = f"({len(notes) + 1})"
            notes.append(f"({len(notes) + 1}) Parameters are elements of the linear form xb_{next(iter(labs))}.")
        else:
            lab = " ".join("{" + P.names[i] + "}" for i in targets)
            mark = ""
        lines.append((lab, dist_t, mark))
    wl = max([len(x[0]) for x in lines] + [1])
    for lab, dist_t, mark in lines:
        txt = f"  {lab:>{wl}} ~ {dist_t}"
        out.write((f"{txt:<75}{mark}" if mark else txt) + "\n", "text")
    out.write(line + "\n", "text")
    for n in notes:
        out.write(n + "\n", "text")


def display_bayes(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                  summary: bool = False, level: float | None = None, hpd: bool | None = None,
                  **kw) -> None:
    out = s.output
    x = est.extra
    model: Model = est.model
    if summary and not x["nomodelsummary"]:
        _model_summary(s, model, x)
    if header:
        out.write("\n", "text")
        eff = x["eff"]

        def num(v):
            return f"{_g4(v):>10}"
        out.write(f"{x['title']:<49}{'MCMC iterations':<17}= {_commas(x['iters']):>10}\n", "text")
        out.write(f"{'Random-walk Metropolis-Hastings sampling':<49}{'Burn-in':<17}= {_commas(x['burn']):>10}\n", "text")
        out.write(f"{'':<49}{'MCMC sample size':<17}= {_commas(x['size']):>10}\n", "text")
        out.write(f"{'':<49}{'Number of obs':<17}= {_commas(x['N']):>10}\n", "text")
        out.write(f"{'':<49}{'Acceptance rate':<17}= {num(x['acc'])}\n", "text")
        out.write(f"{'':<49}{'Efficiency:  min':<17}= {num(min(eff))}\n", "text")
        out.write(f"{'':<49}{'avg':>16} = {num(float(np.mean(eff)))}\n", "text")
        lml = f"Log marginal likelihood = {_g(x['lml'], '%10.0g'):>10}" if np.isfinite(x["lml"]) else \
            "Log marginal likelihood =          ."
        out.write(f"{lml:<49}{'max':>16} = {num(max(eff))}\n", "text")
        out.write(" \n", "text")      # Stata writes a line with one space (case 1909)
    if not table:
        return
    lv = level or x["level"]
    use_hpd = x["hpd"] if hpd is None else hpd
    summ = x["summ"] if (lv == x["level"] and use_hpd == x["hpd"]) else mcmc.summarize(x["draws"], lv, use_hpd)
    _summary_table(s, est, summ, lv, use_hpd)


def _srow(r) -> str:
    """Mean, Std. Dev., MCSE, Median and the interval as Stata 14 writes them (case 1909)."""
    return (f"{_g(r['mean']):>10}{_g(r['sd']):>11}{_g(r['mcse'], '%8.0g'):>10}{_g(r['median']):>11}"
            f"{_g(r['lo']):>11}{_g(r['hi']):>11}")


def _summary_table(s, est, summ, level, use_hpd, names=None) -> None:
    out = s.output
    P = est.model.P if names is None else None
    names = names or list(est.extra["pnames"])
    eqs = list(P.eqs) if P is not None else [""] * len(names)
    W = 12
    dash = "-" * (W + 1) + "+" + "-" * 64 + "\n"
    out.write("-" * (W + 66) + "\n", "text")
    out.write(f"{'':>{W}} |" + " " * 48 + ("         HPD" if use_hpd else "Equal-tailed") + "\n", "text")
    ci = f"[{level:g}% Cred. Interval]"
    out.write(f"{'':>{W}} |      Mean   Std. Dev.     MCSE     Median  {ci}\n", "text")
    out.write(dash, "text")
    order = []
    seen = []
    for i, e in enumerate(eqs):
        if e and e not in seen:
            seen.append(e)
    for e in seen:
        order.append(("head", e, None))
        for i in range(len(names)):
            if eqs[i] == e:
                order.append(("row", names[i].split(":", 1)[1] if ":" in names[i] else names[i], i))
    scal = [i for i in range(len(names)) if not eqs[i]]
    if scal:
        if order:
            order.append(("sep", "", None))
        for i in scal:
            order.append(("row", names[i], i))
    first = True
    for kind, lab, i in order:
        if kind == "head":
            if not first:
                out.write(dash, "text")
            out.write(f"{lab:<{W}} |\n", "text")
        elif kind == "sep":
            out.write(dash, "text")
        else:
            out.write(f"{lab:>{W}} |" + _srow(summ[i]) + "\n", "result")
        first = False
    out.write("-" * (W + 66) + "\n", "text")


# ---------------------------------------------------------------------------
# bayes: prefix
# ---------------------------------------------------------------------------

_PREFIX = {"regress": ("normal", "Bayesian linear regression", True),
           "logit": ("logit", "Bayesian logistic regression", False),
           "probit": ("probit", "Bayesian probit regression", False),
           "poisson": ("poisson", "Bayesian Poisson regression", False)}


def cmd_bayes(s: "Session", args: str) -> None:
    colon = find_top(args, ":")
    if colon == -1 and (not args.strip() or args.strip().startswith(",")):
        _replay(s, args)
        return
    if colon == -1:
        raise StataError(198, "invalid syntax")
    head, command = args[:colon].strip(), args[colon + 1:].strip()
    if not command:
        raise StataError(198, "invalid syntax")
    word = command.split(None, 1)[0]
    hit = [k for k in _PREFIX if k == word] or [k for k in _PREFIX if k.startswith(word) and len(word) >= 3]
    if len(hit) != 1:
        raise StataError(199, f"bayes: {word} is not supported by the bayes prefix in OpenDTA yet")   # OpenDTA
    name = hit[0]
    lik, title, has_var = _PREFIX[name]
    c = find_top(command, ",")
    chead, copts = (command[len(word):], "") if c == -1 else (command[len(word):c], command[c + 1:])
    m = re.search(r"(?:^|[\s)])((?:if|in)\b.*)$", chead, re.S)
    body, tail = (chead[:m.start(1)], m.group(1)) if m else (chead, "")
    words = body.split()
    if not words:
        raise StataError(100, "varlist required")
    noconst = "noconstant" in copts or "nocons" in copts
    hc = find_top(head, ",") if head else -1
    bopts = head.lstrip(", ").strip()
    o = match_options(_strip_option(_strip_option(_strip_option(bopts, "prior"), "block"), "initial"),
                      {**_OPTS, "normalprior": 7, "igammaprior": 6}) if bopts else {}
    sd = float(o["normalprior"]) if o.get("normalprior") else 100.0
    ig = (o["igammaprior"].split() if o.get("igammaprior") else [".01", ".01"])
    dep, regs = words[0], " ".join(words[1:])
    eq = f"{dep} {regs}".strip()
    user_priors = _option_values(bopts, "prior")
    pri = [f"prior({{{dep}:}}, normal(0,{sd * sd:g}))"]
    if has_var:
        pri.append(f"prior({{sigma2}}, igamma({ig[0]},{ig[1]}))")
    pri = [p for p in pri if not any(u.split(",")[0].strip() in p for u in user_priors)] + \
          [f"prior({u})" for u in user_priors]
    lik_opt = f"likelihood(normal({{sigma2}}))" if has_var else f"likelihood({lik})"
    rest = _strip_option(_strip_option(bopts, "normalprior"), "igammaprior")
    rest = _strip_option(_strip_option(rest, "prior"), "likelihood")
    spec = eq + (", noconstant" if noconst else "")
    if noconst:
        spec = f"({eq}, noconstant)"
    new_args = f"{spec} {tail}, {lik_opt} {' '.join(pri)} {rest}"
    lik_names = {dep: name}
    cmd_bayesmh(s, new_args, title=title, lik_names=lik_names)
    s._estimates.extra["prefix"] = name


# ---------------------------------------------------------------------------
# postestimation
# ---------------------------------------------------------------------------

def _current(s) -> Estimates:
    from .results import current
    est = current(s)
    if est.cmd != "bayesmh":
        raise StataError(301, "last estimates not found")
    return est


def _eval_draws(est, text: str) -> np.ndarray:
    """A function of the parameters, evaluated draw by draw: {name} refer to parameters."""
    names = list(est.extra["pnames"])
    draws = est.extra["draws"]
    py = text
    env = {"exp": np.exp, "ln": np.log, "log": np.log, "log10": np.log10, "sqrt": np.sqrt,
           "abs": np.abs, "invlogit": lambda v: 1 / (1 + np.exp(-v)), "logit": lambda v: np.log(v / (1 - v)),
           "normal": lambda v: special.ndtr(v), "invnormal": lambda v: special.ndtri(v)}

    def sub(m):
        body = m.group(1).strip()
        nm = body if body in names else None
        if nm is None:
            raise StataError(198, f"parameter {{{body}}} not found")   # VERIFY
        env[f"_p{names.index(nm)}"] = draws[:, names.index(nm)]
        return f"_p{names.index(nm)}"
    py = _BRACE.sub(sub, py).replace("^", "**")
    try:
        return np.asarray(eval(py, {"__builtins__": {}}, env), float) * np.ones(draws.shape[0])
    except Exception:
        raise StataError(198, f"invalid expression {text}")   # VERIFY


def _param_items(est, text: str):
    """[(label, draws, index or None)] for a list of {parameters} and (expressions)."""
    names = list(est.extra["pnames"])
    draws = est.extra["draws"]
    text = text.strip()
    if not text:
        return [(n, draws[:, i], i) for i, n in enumerate(names)]
    items = []
    i = 0
    while i < len(text):
        ch = text[i]
        if ch.isspace():
            i += 1
        elif ch == "(":
            depth, j = 1, i + 1
            while j < len(text) and depth:
                depth += {"(": 1, ")": -1}.get(text[j], 0)
                j += 1
            e = text[i + 1:j - 1]
            lm = re.match(r"^\s*(\w+)\s*:(?!\w*\})", e)
            lab, body = (lm.group(1), e[lm.end():]) if lm and "{" not in e[:lm.end()] else (e.strip(), e)
            items.append((lab, _eval_draws(est, body), None))
            i = j
        elif ch == "{":
            j = text.index("}", i)
            body = text[i + 1:j]
            i = j + 1
            if ":" in body:
                label, rest = body.split(":", 1)
                label = label.strip()
                sel = [k for k, n in enumerate(names) if n.startswith(label + ":")]
                if rest.strip():
                    want = set(rest.replace(",", " ").split())
                    sel = [k for k in sel if names[k].split(":", 1)[1] in want]
                for k in sel:
                    items.append((names[k], draws[:, k], k))
            else:
                for nm in body.split():
                    if nm not in names:
                        raise StataError(198, f"parameter {{{nm}}} not found")   # VERIFY
                    k = names.index(nm)
                    items.append((nm, draws[:, k], k))
        else:
            raise StataError(198, f"invalid parameter list: {text}")   # VERIFY
    if not items:
        raise StataError(198, "no parameters selected")   # VERIFY
    return items


def cmd_bayesstats(s: "Session", args: str) -> None:
    est = _current(s)
    m = re.match(r"\s*(\w+)\s*(.*)$", args, re.S)
    if not m:
        raise StataError(198, "bayesstats requires summary, ess or ic")   # VERIFY
    sub = m.group(1).lower()
    rest = m.group(2)
    subs = [w for w in ("summary", "ess", "ic") if w.startswith(sub)]
    if len(subs) != 1:
        raise StataError(198, f"bayesstats {sub}: invalid subcommand")   # VERIFY
    c = find_top(rest, ",")
    body, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    o = match_options(opt, {"hpd": 3, "clevel": 4, "skip": 4, "nolog": 3, "all": 3}) if opt.strip() else {}
    level = float(o["clevel"]) if o.get("clevel") else est.extra["level"]
    if subs[0] == "ic":
        _stats_ic(s, est, body, o)
        return
    items = _param_items(est, body)
    draws = np.column_stack([d for _l, d, _i in items])
    if subs[0] == "summary":
        summ = mcmc.summarize(draws, level, bool(o.get("hpd")))
        s.output.write("\nPosterior summary statistics" + " " * 22 + f"MCMC sample size = {_commas(draws.shape[0]):>9}\n \n", "text")
        _flat_table(s, [l for l, _d, _i in items], summ, level, bool(o.get("hpd")))
    else:
        s.output.write("\nEfficiency summaries    " + f"MCMC sample size = {_commas(draws.shape[0]):>9}\n \n", "text")
        out = s.output
        W = 12
        out.write("-" * (W + 40) + "\n", "text")
        out.write(f"{'':>{W}} |        ESS   Corr. time    Efficiency\n", "text")
        out.write("-" * (W + 1) + "+" + "-" * 38 + "\n", "text")

        def one(lab, d):
            e = mcmc.ess(d)
            out.write(f"{lab:>{W}} |", "text")
            out.write(f"{e:>11.2f}{len(d) / e:>13.2f}{e / len(d):>14.4f}\n", "result")
        _grouped(out, [l for l, _d, _i in items], [d for _l, d, _i in items], one, W, 38)
        out.write("-" * (W + 40) + "\n", "text")


def _grouped(out, labels, payload, one, W, width):
    """Rows of eq:name parameters under an equation header, the others after a separator, as
    Stata's bayesstats tables (case 1909)."""
    seen, done = [], False
    for lab in labels:
        if ":" in lab and lab.split(":", 1)[0] not in seen:
            seen.append(lab.split(":", 1)[0])
    for k, e in enumerate(seen):
        if k:
            out.write("-" * (W + 1) + "+" + "-" * width + "\n", "text")
        out.write(f"{e:<{W}} |\n", "text")
        for lab, pl in zip(labels, payload):
            if ":" in lab and lab.split(":", 1)[0] == e:
                one(lab.split(":", 1)[1], pl)
        done = True
    rest = [(lab, pl) for lab, pl in zip(labels, payload) if ":" not in lab]
    if rest and done:
        out.write("-" * (W + 1) + "+" + "-" * width + "\n", "text")
    for lab, pl in rest:
        one(lab, pl)


def _flat_table(s, labels, summ, level, use_hpd) -> None:
    out = s.output
    W = 12
    out.write("-" * (W + 66) + "\n", "text")
    out.write(f"{'':>{W}} |" + " " * 48 + ("         HPD" if use_hpd else "Equal-tailed") + "\n", "text")
    out.write(f"{'':>{W}} |      Mean   Std. Dev.     MCSE     Median  [{level:g}% Cred. Interval]\n", "text")
    out.write("-" * (W + 1) + "+" + "-" * 64 + "\n", "text")

    def one(lab, r):
        out.write(f"{lab[-W:]:>{W}} |" + _srow(r) + "\n", "result")
    _grouped(out, labels, summ, one, W, 64)
    out.write("-" * (W + 66) + "\n", "text")


def _stored_models(s, names: str):
    from .results import current
    if not names.strip():
        est = _current(s)
        return [("active", est)]
    out = []
    store = getattr(s, "_stored_estimates", {})
    for n in names.split():
        if n == ".":
            out.append(("active", _current(s)))
            continue
        snap = store.get(n)
        if snap is None:
            raise StataError(111, f"estimation result {n} not found")
        est = snap["est"]
        if est is None or est.cmd != "bayesmh":
            raise StataError(198, f"estimation result {n} is not from bayesmh")   # VERIFY
        out.append((n, est))
    return out


def _stats_ic(s, est, body: str, o) -> None:
    models = _stored_models(s, body)
    out = s.output
    out.write("\nBayesian information criteria\n\n", "text")
    W = 12
    out.write("-" * (W + 36) + "\n", "text")
    out.write(f"{'':>{W}} |   Avg. log DIC  log(ML)   log(BF)\n".replace("Avg. log DIC", "       DIC"), "text")
    out.write("-" * (W + 1) + "+" + "-" * 34 + "\n", "text")
    base = models[0][1].extra["lml"]
    for name, e in models:
        x = e.extra
        bf = x["lml"] - base
        out.write(f"{name:>{W}} |  {_g(x['dic']):>9}  {_g(x['lml']):>8}  {_g(bf):>8}\n", "result")
    out.write("-" * (W + 36) + "\n", "text")
    out.write("Note: Marginal likelihood (ML) is computed using Laplace-Metropolis approximation.\n", "text")


def cmd_bayestest(s: "Session", args: str) -> None:
    est = _current(s)
    m = re.match(r"\s*(\w+)\s*(.*)$", args, re.S)
    if not m:
        raise StataError(198, "bayestest requires interval or model")   # VERIFY
    sub = m.group(1).lower()
    rest = m.group(2)
    subs = [w for w in ("interval", "model") if w.startswith(sub)]
    if len(subs) != 1:
        raise StataError(198, f"bayestest {sub}: invalid subcommand")   # VERIFY
    c = find_top(rest, ",")
    body, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    out = s.output
    if subs[0] == "interval":
        o = match_options(opt, {"lower": 3, "upper": 3, "prior": 5}) if opt.strip() else {}
        items = _param_items(est, body)
        n = len(items[0][1])
        out.write("\nInterval tests     " + f"MCMC sample size = {_commas(n):>9}\n\n", "text")
        lo = float(o["lower"]) if o.get("lower") else -np.inf
        hi = float(o["upper"]) if o.get("upper") else np.inf
        res = []
        for k, (lab, d, _i) in enumerate(items, 1):
            cond = (f"{_g(lo)} < {{{lab}}} < {_g(hi)}" if np.isfinite(lo) and np.isfinite(hi) else
                    f"{{{lab}}} > {_g(lo)}" if np.isfinite(lo) else f"{{{lab}}} < {_g(hi)}")
            out.write(f"{'prob' + str(k):>12} : {cond}\n", "text")
            ind = ((d > lo) & (d < hi)).astype(float)
            p = float(ind.mean())
            sd = float(np.std(ind, ddof=1))
            e = mcmc.ess(ind) if sd > 0 else float(len(ind))
            res.append((f"prob{k}", p, sd, sd / math.sqrt(e) if e > 0 else 0.0))
        W = 12
        out.write(" \n" + "-" * (W + 35) + "\n", "text")
        out.write(f"{'':>{W}} |      Mean    Std. Dev.      MCSE\n", "text")
        out.write("-" * (W + 1) + "+" + "-" * 33 + "\n", "text")
        for lab, p, sd, mc in res:
            out.write(f"{lab:>{W}} |", "text")
            out.write(f"{_g(p):>10}{sd:>12.5f}{_g(mc):>11}\n", "result")
        out.write("-" * (W + 35) + "\n", "text")
        return
    models = _stored_models(s, body)
    if len(models) < 2:
        raise StataError(198, "bayestest model requires at least two models")   # VERIFY
    lml = np.array([e.extra["lml"] for _n, e in models])
    if not np.all(np.isfinite(lml)):
        raise StataError(498, "log marginal likelihood is not available for these models")   # VERIFY
    pri = np.full(len(models), 1.0 / len(models))
    post_ = pri * np.exp(lml - lml.max())
    post_ /= post_.sum()
    out.write("\nBayesian model tests\n\n", "text")
    W = 12
    out.write("-" * (W + 36) + "\n", "text")
    out.write(f"{'':>{W}} |    log(ML)    P(M)  P(M|y)\n", "text")
    out.write("-" * (W + 1) + "+" + "-" * 34 + "\n", "text")
    for (n, e), pr, po in zip(models, pri, post_):
        out.write(f"{n:>{W}} |  {_g(e.extra['lml']):>9}  {pr:>6.4f}  {po:>6.4f}\n", "result")
    out.write("-" * (W + 36) + "\n", "text")


def cmd_bayesgraph(s: "Session", args: str) -> None:
    est = _current(s)
    m = re.match(r"\s*(\w+)\s*(.*)$", args, re.S)
    if not m:
        raise StataError(198, "bayesgraph requires a graph type")   # VERIFY
    sub = m.group(1).lower()
    kinds = ("diagnostics", "trace", "histogram", "kdensity", "ac", "cusum")
    hit = [k for k in kinds if k.startswith(sub)]
    if len(hit) != 1:
        raise StataError(198, f"bayesgraph {sub}: invalid graph type")   # VERIFY
    kind = hit[0]
    rest = m.group(2)
    c = find_top(rest, ",")
    body, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    items = _param_items(est, body)
    if len(items) > 1 and kind != "diagnostics":
        pass
    for lab, d, _i in items:
        if kind == "diagnostics":
            names = []
            for k in ("trace", "histogram", "ac", "kdensity"):
                nm = f"__bg_{k}"
                _one_graph(s, k, lab, d, opt, name=nm)
                names.append(nm)
            s.run(f"graph combine {' '.join(names)}")
        else:
            _one_graph(s, kind, lab, d, opt)


def _one_graph(s, kind: str, lab: str, d: np.ndarray, opt: str, name: str | None = None) -> None:
    from ..core.dataset import Variable
    n = len(d)
    ttl = lab
    s.run("preserve")
    try:
        s.run("quietly clear")
        if kind == "trace":
            s.run(f"quietly set obs {n}")
            s.data.add(Variable("iter", "long", np.arange(1, n + 1, dtype=float)))
            s.data.add(Variable("theta", "double", np.asarray(d, float)))
            cmd = f'twoway line theta iter, title("{ttl}") xtitle("Iteration number") ytitle("") legend(off)'
        elif kind == "histogram":
            s.run(f"quietly set obs {n}")
            s.data.add(Variable("theta", "double", np.asarray(d, float)))
            cmd = f'histogram theta, title("{ttl}") xtitle("") legend(off)'
        elif kind == "kdensity":
            s.run(f"quietly set obs {n}")
            s.data.add(Variable("theta", "double", np.asarray(d, float)))
            cmd = f'kdensity theta, title("{ttl}") xtitle("") legend(off) note("")'
        elif kind == "ac":
            ml = min(40, n - 1)
            ac = mcmc.autocorr(np.asarray(d, float), ml)
            s.run(f"quietly set obs {ml}")
            s.data.add(Variable("lag", "long", np.arange(1, ml + 1, dtype=float)))
            s.data.add(Variable("ac", "double", ac[1:]))
            cmd = f'twoway (spike ac lag), title("{ttl}") xtitle("Lag") ytitle("Autocorrelation") legend(off)'
        else:
            s.run(f"quietly set obs {n}")
            cs = np.cumsum(d - np.mean(d))
            s.data.add(Variable("iter", "long", np.arange(1, n + 1, dtype=float)))
            s.data.add(Variable("cusum", "double", cs))
            cmd = f'twoway line cusum iter, title("{ttl}") xtitle("Iteration number") ytitle("Cusum") legend(off)'
        if name:
            cmd += f" name({name}, replace)"
        s.run(cmd)
    finally:
        s.run("restore")

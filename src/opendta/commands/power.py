"""power ([PSS] power): power and sample-size analysis.

Methods: onemean, twomeans, pairedmeans, oneproportion, twoproportions, onecorrelation and
twocorrelations. Each one computes the sample size (default), the power (when n() is given)
or the effect size (when n() and power() are given), from the formulas of [PSS] (Methods and
formulas of each method). A list of values in any parameter, or the table option, gives a table
with one row per combination.

The layout of the output (column widths, the order of the parameters and of the table columns)
is not defined in the manual and is marked VERIFY (case compat/do/1301_power.do)."""

from __future__ import annotations

import itertools
import math
from typing import TYPE_CHECKING, Callable

from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import find_top, parse_options
from ..lang.words import parse_numlist
from ..lazy import lazy
from .registry import command

st = lazy("scipy.stats")
opt_ = lazy("scipy.optimize")

if TYPE_CHECKING:
    from ..session import Session


# ---------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------

_METHODS = {"onemean": 4, "twomeans": 4, "pairedmeans": 5, "oneproportion": 4,
            "twoproportions": 4, "onecorrelation": 4, "twocorrelations": 4}

# options with values; the abbreviations follow the manual's underlining (VERIFY)
_OPTS = {"alpha": 1, "power": 1, "beta": 1, "n": 1, "n1": 2, "n2": 2, "nratio": 2,
         "sd": 2, "sd1": 3, "sd2": 3, "knownsd": 6, "knownsds": 6, "diff": 2, "ratio": 2,
         "sddiff": 3, "corr": 2, "onesided": 3, "direction": 3, "nfractional": 4, "table": 3,
         "notable": 5, "test": 4, "continuity": 4, "effect": 3, "parallel": 3, "graph": 2,
         "init": 4, "iterate": 4, "tolerance": 3, "ftolerance": 4, "nolog": 4, "log": 3,
         "dots": 4, "noheader": 6, "saving": 4}


def _values(text: str) -> list[float]:
    t = text.strip()
    if t.startswith("(") and t.endswith(")"):
        t = t[1:-1]
    try:
        vals = parse_numlist(t)
    except StataError:
        raise StataError(121, "invalid numlist")
    return [float(v) for v in vals]


def _split_args(body: str) -> list[str]:
    """Positional arguments: numbers or (numlist) in parentheses."""
    out, i, n = [], 0, len(body)
    while i < n:
        if body[i].isspace():
            i += 1
            continue
        if body[i] == "(":
            depth, j = 0, i
            while j < n:
                depth += body[j] == "("
                depth -= body[j] == ")"
                if depth == 0:
                    break
                j += 1
            if j >= n:
                raise StataError(198, "unmatched parenthesis")
            out.append(body[i:j + 1])
            i = j + 1
            continue
        j = i
        while j < n and not body[j].isspace() and body[j] != "(":
            j += 1
        out.append(body[i:j])
        i = j
    return out


def _parse(args: str):
    t = args.strip()
    if not t:
        raise StataError(198, "method required")   # VERIFY
    head, _, rest = t.partition(" ")
    if "," in head:
        head, _, tail = head.partition(",")
        rest = "," + tail + " " + rest
    method = next((m for m, k in _METHODS.items() if len(head) >= k and m.startswith(head)), None)
    if method is None:
        raise StataError(198, f"method {head} not found")   # VERIFY
    c = find_top(rest, ",")
    body, opts_text = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    pos = [_values(a) for a in _split_args(body)]
    opts: dict[str, object] = {}
    for name, val in parse_options(opts_text) if opts_text.strip() else []:
        low = name.lower()
        full = next((o for o, k in _OPTS.items() if len(low) >= k and o.startswith(low)), None)
        if full is None:
            raise StataError(198, f"option {name} not allowed")
        opts[full] = True if val is None else val
    return method, pos, opts


# ---------------------------------------------------------------------------
# distributions
# ---------------------------------------------------------------------------

def _zpower(effect_z: float, alpha: float, sides: int) -> float:
    """Power of a z test whose statistic has mean effect_z under Ha."""
    if sides == 1:
        return float(st.norm.sf(st.norm.isf(alpha) - effect_z))
    z = float(st.norm.isf(alpha / 2))
    return float(st.norm.sf(z - effect_z) + st.norm.cdf(-z - effect_z))


def _tpower(nc: float, df: float, alpha: float, sides: int) -> float:
    if df <= 0:
        return math.nan
    if sides == 1:
        crit = float(st.t.isf(alpha, df))
        val = float(st.nct.sf(crit, df, nc))
    else:
        crit = float(st.t.isf(alpha / 2, df))
        val = float(st.nct.sf(crit, df, nc) + st.nct.cdf(-crit, df, nc))
    if not math.isfinite(val):
        # far tails of the noncentral t (large df and noncentrality): normal approximation
        val = float(st.norm.sf(crit - nc)) + (0.0 if sides == 1 else float(st.norm.cdf(-crit - nc)))
    return min(max(val, 0.0), 1.0)


def _solve(f: Callable[[float], float], lo: float, hi: float) -> float:
    """Root of an increasing function f on [lo, hi] (hi grown as needed)."""
    flo = f(lo)
    if flo >= 0:
        return lo
    k = 0
    while f(hi) < 0:
        hi *= 2
        k += 1
        if k > 60:
            raise StataError(498, "failed to find a solution")   # VERIFY
    return float(opt_.brentq(f, lo, hi, xtol=1e-12, rtol=1e-12, maxiter=500))


# ---------------------------------------------------------------------------
# methods: each returns the power for a given configuration (n = total sample size or the
# size of the first group), and knows how to describe itself
# ---------------------------------------------------------------------------

class _Method:
    title = ""
    two_sample = False
    # Stata 14 prints "Performing iteration ..." whenever it solves for n or the effect, also
    # for closed forms (compat 1301; VERIFY the correlation methods)
    iterate = True
    effect_name = "ma"
    null_name = "m0"

    def right(self, c) -> list:
        """Parameters shown in a second column (pairedmeans)."""
        return []

    def __init__(self, pos: list[float], o: dict):
        self.o = o

    def power(self, c: dict) -> float:
        raise NotImplementedError

    def sizes(self, c: dict, n: float) -> dict:
        """Group sizes for total n (two-sample methods: n1 and n2 from nratio)."""
        return {"N": n}


def _sides(c: dict) -> int:
    return 1 if c.get("onesided") else 2


def _direction(c: dict) -> int:
    return -1 if str(c.get("direction", "upper")).lower().startswith("l") else 1


class _OneMean(_Method):
    title = "a one-sample mean test"
    hyp = "Ho: m = m0  versus  Ha: m {op} m0"
    effect_name = "ma"

    def test(self, c):
        return "z test" if c.get("knownsd") else "t test"

    def delta(self, c):
        return (c["ma"] - c["m0"]) / c["sd"]

    def power(self, c):
        n = c["N"]
        d = abs(self.delta(c))
        if c.get("knownsd"):
            return _zpower(d * math.sqrt(n), c["alpha"], _sides(c))
        return _tpower(d * math.sqrt(n), n - 1, c["alpha"], _sides(c))

    def set_effect(self, c, delta):
        c["ma"] = c["m0"] + _direction(c) * delta * c["sd"]

    def params(self, c):
        return [("delta", self.delta(c)), ("m0", c["m0"]), ("ma", c["ma"]), ("sd", c["sd"])]

    # Stata 14 prints "Performing iteration ..." for the z test as well (compat 1301)
    iterate = True


class _TwoMeans(_Method):
    title = "a two-sample means test"
    hyp = "Ho: m2 = m1  versus  Ha: m2 {op} m1"
    two_sample = True
    effect_name = "m2"

    def test(self, c):
        if c.get("knownsds"):
            return "z test"
        if c["sd1"] != c["sd2"]:
            return "Satterthwaite's t test assuming unequal variances"
        return "t test assuming sd1 = sd2 = sd"

    def delta(self, c):
        return c["m2"] - c["m1"]

    def sizes(self, c, n1):
        r = c.get("nratio", 1.0)
        return {"N1": n1, "N2": r * n1, "N": n1 * (1 + r)}

    def power(self, c):
        n1, n2 = c["N1"], c["N2"]
        s1, s2 = c["sd1"], c["sd2"]
        se = math.sqrt(s1 ** 2 / n1 + s2 ** 2 / n2)
        d = abs(self.delta(c))
        if c.get("knownsds"):
            return _zpower(d / se, c["alpha"], _sides(c))
        if s1 == s2:
            df = n1 + n2 - 2
        else:
            v1, v2 = s1 ** 2 / n1, s2 ** 2 / n2
            df = (v1 + v2) ** 2 / (v1 ** 2 / (n1 - 1) + v2 ** 2 / (n2 - 1))
        return _tpower(d / se, df, c["alpha"], _sides(c))

    def set_effect(self, c, delta):
        c["m2"] = c["m1"] + _direction(c) * delta

    def params(self, c):
        rows = [("delta", self.delta(c)), ("m1", c["m1"]), ("m2", c["m2"])]
        if c["sd1"] == c["sd2"]:
            rows.append(("sd", c["sd1"]))
        else:
            rows += [("sd1", c["sd1"]), ("sd2", c["sd2"])]
        if c.get("nratio", 1.0) != 1.0:
            rows.append(("N2/N1", c["nratio"]))   # compat 1301
        return rows



class _PairedMeans(_Method):
    title = "a two-sample paired-means test"
    hyp = "Ho: d = d0  versus  Ha: d {op} d0"
    effect_name = "ma2"

    def test(self, c):
        base = "Paired z test" if c.get("knownsd") else "Paired t test"   # VERIFY knownsd
        if c.get("corr") is not None and c.get("sd1") == c.get("sd2"):
            base += " assuming sd1 = sd2 = sd"
        return base

    def sd_d(self, c):
        if c.get("sddiff") is not None:
            return c["sddiff"]
        if c.get("sd1") is not None and c.get("corr") is not None:
            s1, s2, r = c["sd1"], c["sd2"], c["corr"]
            return math.sqrt(s1 ** 2 + s2 ** 2 - 2 * r * s1 * s2)
        raise StataError(198, "sddiff(), or corr() with sd() or sd1() and sd2(), required")   # VERIFY

    def delta(self, c):
        return (c["ma2"] - c["ma1"]) / self.sd_d(c)

    def power(self, c):
        n = c["N"]
        d = abs(self.delta(c))
        if c.get("knownsd"):
            return _zpower(d * math.sqrt(n), c["alpha"], _sides(c))
        return _tpower(d * math.sqrt(n), n - 1, c["alpha"], _sides(c))

    def set_effect(self, c, delta):
        c["ma2"] = c["ma1"] + _direction(c) * delta * self.sd_d(c)

    def params(self, c):
        return [("delta", self.delta(c)), ("d0", 0.0), ("da", c["ma2"] - c["ma1"]),
                ("sd_d", self.sd_d(c))]

    def right(self, c):
        rows = [("ma1", c["ma1"]), ("ma2", c["ma2"])]
        if c.get("corr") is not None:
            rows += [("sd", c["sd1"])] if c["sd1"] == c["sd2"] else [("sd1", c["sd1"]), ("sd2", c["sd2"])]
            rows.append(("corr", c["corr"]))
        return rows



class _OneProportion(_Method):
    title = "a one-sample proportion test"
    hyp = "Ho: p = p0  versus  Ha: p {op} p0"
    effect_name = "pa"

    def kind(self):
        t = str(self.o.get("test") or "score").lower()
        for k in ("score", "wald", "binomial"):
            if k.startswith(t):
                return k
        raise StataError(198, f"test({t}) not allowed")   # VERIFY

    def test(self, c):
        return {"score": "Score z test", "wald": "Wald z test", "binomial": "Binomial test"}[self.kind()]

    def delta(self, c):
        return c["pa"] - c["p0"]

    def power(self, c):
        n, p0, pa, a = c["N"], c["p0"], c["pa"], c["alpha"]
        sides = _sides(c)
        kind = self.kind()
        if kind == "binomial":
            pw, c["actual_alpha"] = _binom_power(int(round(n)), p0, pa, a, sides, 1 if pa >= p0 else -1)
            return pw
        s0 = math.sqrt(p0 * (1 - p0)) if kind == "score" else math.sqrt(pa * (1 - pa))
        sa = math.sqrt(pa * (1 - pa))
        d = abs(pa - p0)
        if sides == 1:
            z = float(st.norm.isf(a))
            return float(st.norm.cdf((d * math.sqrt(n) - z * s0) / sa))
        z = float(st.norm.isf(a / 2))
        return float(st.norm.cdf((d * math.sqrt(n) - z * s0) / sa)
                     + st.norm.cdf((-d * math.sqrt(n) - z * s0) / sa))

    def set_effect(self, c, delta):
        c["pa"] = c["p0"] + _direction(c) * delta

    def params(self, c):
        return [("delta", self.delta(c)), ("p0", c["p0"]), ("pa", c["pa"])]


def _binom_power(n: int, p0: float, pa: float, alpha: float, sides: int,
                 direction: int) -> tuple[float, float]:
    """Exact binomial test: rejection region of size <= alpha (equal tails when two-sided);
    returns the power and the actual size of the region (compat 1301)."""
    if n < 1:
        return 0.0, 0.0

    def prob(lo, hi, p):
        return ((float(st.binom.cdf(lo, n, p)) if lo >= 0 else 0.0)
                + (float(st.binom.sf(hi - 1, n, p)) if hi <= n else 0.0))
    ks = range(n + 1)
    if sides == 1:
        if direction > 0:
            lo, hi = -1, next((k for k in ks if st.binom.sf(k - 1, n, p0) <= alpha), n + 1)
        else:
            lo, hi = next((k for k in reversed(ks) if st.binom.cdf(k, n, p0) <= alpha), -1), n + 1
    else:
        lo = next((k for k in reversed(ks) if st.binom.cdf(k, n, p0) <= alpha / 2), -1)
        hi = next((k for k in ks if st.binom.sf(k - 1, n, p0) <= alpha / 2), n + 1)
    return prob(lo, hi, pa), prob(lo, hi, p0)


class _TwoProportions(_Method):
    title = "a two-sample proportions test"
    hyp = "Ho: p2 = p1  versus  Ha: p2 {op} p1"
    two_sample = True
    effect_name = "p2"

    def test(self, c):
        t = str(self.o.get("test") or "chi2").lower()
        if t not in ("chi2",):
            raise StataError(198, f"test({t}) not allowed; OpenDTA computes the Pearson "
                                  "chi-squared test")   # OpenDTA: lrchi2 and fisher
        return "Pearson's chi-squared test"

    def delta(self, c):
        return c["p2"] - c["p1"]

    def sizes(self, c, n1):
        r = c.get("nratio", 1.0)
        return {"N1": n1, "N2": r * n1, "N": n1 * (1 + r)}

    def power(self, c):
        n1, n2, p1, p2, a = c["N1"], c["N2"], c["p1"], c["p2"], c["alpha"]
        pbar = (n1 * p1 + n2 * p2) / (n1 + n2)
        s0 = math.sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))
        sa = math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
        d = abs(p2 - p1)
        if c.get("continuity"):
            d = max(d - 0.5 * (1 / n1 + 1 / n2), 0.0)
        if _sides(c) == 1:
            z = float(st.norm.isf(a))
            return float(st.norm.cdf((d - z * s0) / sa))
        z = float(st.norm.isf(a / 2))
        return float(st.norm.cdf((d - z * s0) / sa) + st.norm.cdf((-d - z * s0) / sa))

    def set_effect(self, c, delta):
        c["p2"] = c["p1"] + _direction(c) * delta

    delta_note = "  (difference)"     # compat 1301

    def params(self, c):
        rows = [("delta", self.delta(c)), ("p1", c["p1"]), ("p2", c["p2"])]
        if c.get("nratio", 1.0) != 1.0:
            rows.append(("N2/N1", c["nratio"]))   # compat 1301
        return rows


class _OneCorrelation(_Method):
    title = "a one-sample correlation test"
    hyp = "Ho: r = r0  versus  Ha: r {op} r0"
    effect_name = "ra"

    def test(self, c):
        return "Fisher's z test"

    def delta(self, c):
        return c["ra"] - c["r0"]

    def power(self, c):
        n = c["N"]
        if n <= 3:
            return 0.0
        d = math.atanh(c["ra"]) - math.atanh(c["r0"])
        return _zpower(abs(d) * math.sqrt(n - 3), c["alpha"], _sides(c))

    def set_effect(self, c, delta):
        c["ra"] = c["r0"] + _direction(c) * delta

    def params(self, c):
        return [("delta", self.delta(c)), ("r0", c["r0"]), ("ra", c["ra"])]


class _TwoCorrelations(_Method):
    title = "a two-sample correlations test"
    hyp = "Ho: r2 = r1  versus  Ha: r2 {op} r1"
    two_sample = True
    effect_name = "r2"

    def test(self, c):
        return "Fisher's z test"

    def delta(self, c):
        return c["r2"] - c["r1"]

    def sizes(self, c, n1):
        r = c.get("nratio", 1.0)
        return {"N1": n1, "N2": r * n1, "N": n1 * (1 + r)}

    def power(self, c):
        n1, n2 = c["N1"], c["N2"]
        if min(n1, n2) <= 3:
            return 0.0
        d = math.atanh(c["r2"]) - math.atanh(c["r1"])
        return _zpower(abs(d) / math.sqrt(1 / (n1 - 3) + 1 / (n2 - 3)), c["alpha"], _sides(c))

    def set_effect(self, c, delta):
        c["r2"] = c["r1"] + _direction(c) * delta

    def params(self, c):
        return [("delta", self.delta(c)), ("r1", c["r1"]), ("r2", c["r2"])]


_CLASSES = {"onemean": _OneMean, "twomeans": _TwoMeans, "pairedmeans": _PairedMeans,
            "oneproportion": _OneProportion, "twoproportions": _TwoProportions,
            "onecorrelation": _OneCorrelation, "twocorrelations": _TwoCorrelations}

# positional arguments of each method: (null value name, alternative value name, defaults)
_POSITIONAL = {"onemean": ("m0", "ma"), "twomeans": ("m1", "m2"), "pairedmeans": ("ma1", "ma2"),
               "oneproportion": ("p0", "pa"), "twoproportions": ("p1", "p2"),
               "onecorrelation": ("r0", "ra"), "twocorrelations": ("r1", "r2")}


# ---------------------------------------------------------------------------
# configurations
# ---------------------------------------------------------------------------

def _num(o: dict, name: str) -> list[float] | None:
    v = o.get(name)
    if v is None or v is True:
        return None
    return _values(str(v))


def _configs(method: str, pos: list[list[float]], o: dict):
    """All combinations of the given values: list of dicts."""
    a_name, b_name = _POSITIONAL[method]
    base: dict[str, list] = {}
    if len(pos) > 2:
        raise StataError(198, "too many values specified")   # VERIFY
    if len(pos) == 2:
        base[a_name], base[b_name] = pos
    elif len(pos) == 1:
        # only the null value: the alternative comes from diff() or is computed
        base[a_name], base[b_name] = pos[0], None
    else:
        raise StataError(198, f"{a_name} required")   # VERIFY
    for name in ("alpha", "power", "beta", "n", "n1", "n2", "nratio", "sd", "sd1", "sd2", "diff",
                 "corr", "sddiff", "ratio", "rdiff", "rrisk", "oratio", "prdiscordant", "sum"):
        v = _num(o, name)
        if v is not None:
            base[name] = v
    keys = list(base)
    lists = [base[k] if base[k] is not None else [None] for k in keys]
    out = []
    for combo in itertools.product(*lists):
        c = dict(zip(keys, combo))
        for flag in ("onesided", "knownsd", "knownsds", "continuity", "nfractional"):
            if o.get(flag):
                c[flag] = True
        if isinstance(o.get("direction"), str):
            c["direction"] = o["direction"]
        if method in ("onevariance", "twovariances") and o.get("sd"):
            c["sdscale"] = True
            for k in (a_name, b_name):
                if c.get(k) is not None:
                    c[k] = c[k] ** 2
        if c.get("rdiff") is not None and method != "twoproportions":
            c["diff"] = c.pop("rdiff")
        out.append(c)
    return out, any(len(x) > 1 for x in lists)


def _complete(method: str, c: dict) -> str:
    """Fills defaults and decides what to compute: 'n', 'power' or 'effect'."""
    a_name, b_name = _POSITIONAL[method]
    c.setdefault("alpha", 0.05)
    if c.get("beta") is not None:
        c["power"] = 1 - c["beta"]
    if method in ("onemean", "twomeans"):
        sd = c.get("sd")
        if method == "twomeans":
            c["sd1"] = c.get("sd1") if c.get("sd1") is not None else (sd if sd is not None else 1.0)
            c["sd2"] = c.get("sd2") if c.get("sd2") is not None else (sd if sd is not None else 1.0)
        else:
            c["sd"] = sd if sd is not None else 1.0
    if method == "pairedmeans":
        if c.get("sd") is not None:
            c.setdefault("sd1", c["sd"])
            c.setdefault("sd2", c["sd"])
        if c.get("sd1") is not None and c.get("sd2") is None:
            c["sd2"] = c["sd1"]
    prep = getattr(_CLASSES[method], "prepare", None)
    if prep is not None and (c.get(b_name) is None or method == "pairedproportions"):
        prep(c)
    if c.get(b_name) is None and c.get("diff") is not None:
        c[b_name] = c[a_name] + c["diff"]
    have_n = any(c.get(k) is not None for k in ("n", "n1", "n2"))
    if c.get(b_name) is None:
        if not have_n:
            raise StataError(198, f"{b_name} required")   # VERIFY
        if c.get("power") is None:
            # with n() and no power(), Stata 14 computes power and asks for the effect (compat 1301)
            if method == "twoproportions":
                raise StataError(198, "one of options diff(), rdiff(), ratio(), rrisk(), or oratio() or the "
                                      "experimental-group proportion must be specified")
            raise StataError(198, f"{_ALT_TEXT[method]} is required to compute power")
        return "effect"
    if have_n:
        if c.get("power") is not None:
            raise StataError(198, "only two of n(), power() and the effect size may be specified")
        return "power"
    c.setdefault("power", 0.8)
    return "n"


_ALT_TEXT = {"onemean": "alternative mean", "twomeans": "experimental-group mean",
             "pairedmeans": "alternative mean difference", "oneproportion": "alternative proportion",
             "twoproportions": "experimental-group proportion",
             "onecorrelation": "alternative correlation",
             "twocorrelations": "experimental-group correlation"}   # VERIFY: all but onemean


def _set_sizes(M: _Method, c: dict) -> None:
    if not M.two_sample:
        c["N"] = c["n"]
        return
    r = c.get("nratio")
    n, n1, n2 = c.get("n"), c.get("n1"), c.get("n2")
    if n1 is not None and n2 is not None:
        c["N1"], c["N2"] = n1, n2
    elif n1 is not None:
        c["N1"] = n1
        c["N2"] = n1 * (r or 1.0)
    elif n2 is not None:
        c["N2"] = n2
        c["N1"] = n2 / (r or 1.0)
    else:
        rr = r or 1.0
        c["N1"] = n / (1 + rr)
        c["N2"] = n - c["N1"]
    c["N"] = c["N1"] + c["N2"]
    c["nratio"] = c["N2"] / c["N1"]


def _compute(M: _Method, method: str, c: dict, what: str) -> None:
    chk = getattr(M, "check_mode", None)
    if chk is not None:
        chk(what)
    if what == "power":
        _set_sizes(M, c)
        c["power_out"] = M.power(c)
        return
    if what == "n":
        target = c["power"]
        if M.two_sample:
            r = c.get("nratio", 1.0)
            c["nratio"] = r

            def f(n1):
                c.update(M.sizes(c, n1))
                return M.power(c) - target
            lo = 4.0 if method != "twocorrelations" else 3.0001
            n1 = _solve(f, lo, 1000.0)
            if method == "twoproportions":
                # starting value: the closed-form size of the first group, rounded up (compat 1301)
                z, zb = st.norm.ppf(1 - c["alpha"] / _sides(c)), st.norm.ppf(target)
                p1, p2 = c["p1"], c["p2"]
                pb = (p1 + r * p2) / (1 + r)
                num = z * math.sqrt(pb * (1 - pb) * (1 + 1 / r)) + zb * math.sqrt(p1 * (1 - p1) + p2 * (1 - p2) / r)
                c["_init"] = float(math.ceil(num ** 2 / (p2 - p1) ** 2))
            if method == "twomeans":
                # starting value: the normal-approximation size of the first group, rounded up
                # (compat 1301: 63)
                z = st.norm.ppf(1 - c["alpha"] / _sides(c)) + st.norm.ppf(target)
                se2 = c["sd1"] ** 2 + c["sd2"] ** 2 / r
                c["_init"] = float(math.ceil((z ** 2) * se2 / M.delta(c) ** 2))
            if not c.get("nfractional"):
                n1 = math.ceil(n1 - 1e-9)
                n2 = math.ceil(r * n1 - 1e-9)
            else:
                n2 = r * n1
            c["N1"], c["N2"], c["N"] = n1, n2, n1 + n2
        else:
            def f(n):
                c["N"] = n
                return M.power(c) - target
            lo = 2.0 if method not in ("onecorrelation",) else 3.0001
            if method == "oneproportion" and M.kind() == "binomial":
                # compat 1301
                raise StataError(198, "power oneproportion: sample-size determination is not allowed.\n"
                                      "    Sample-size determination is not allowed for the binomial test.  "
                                      "Power function of the binomial test is often nonmonotonic; see "
                                      "Saw-toothed power function for details.")
            else:
                n = _solve(f, lo, 1000.0)
                if method == "onemean":
                    # starting value: the normal-approximation sample size (compat 1301)
                    z = st.norm.ppf(1 - c["alpha"] / _sides(c)) + st.norm.ppf(target)
                    c["_init"] = float((z / abs(M.delta(c))) ** 2)
                if not c.get("nfractional"):
                    n = math.ceil(n - 1e-9)
            c["N"] = n
        c["power_out"] = M.power(c)
        return
    # effect size: the smallest effect with the requested power
    _set_sizes(M, c)
    target = c["power"]
    a_name, _b = _POSITIONAL[method]

    def g(d):
        M.set_effect(c, d)
        return M.power(c) - target
    base = c[a_name]
    up = _direction(c) > 0
    if "proportion" in method:
        top = (1 - base if up else base) - 1e-9
    elif "correlation" in method:
        top = (1 - base if up else base + 1) - 1e-9
    else:
        top = None
    if top is None:
        d = _solve(g, 1e-9, 50.0)
    elif top > 1e-9 and g(top) > 0:
        d = float(opt_.brentq(g, 1e-9, top, xtol=1e-12))
    else:
        raise StataError(498, "failed to find a solution: the requested power cannot be "
                              "reached")   # VERIFY
    M.set_effect(c, d)
    c["power_out"] = M.power(c)


# ---------------------------------------------------------------------------
# display
# ---------------------------------------------------------------------------

def _f4(x) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "."
    return f"{x:9.4f}"


def _gtab(x) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "."
    x = float(x)
    if not x.is_integer():
        x = float(f"{x:.4g}")       # four significant digits (compat 1301)
    return format_value(x, "%7.0g", pad=False).strip()


def _what_title(M: _Method, what: str) -> str:
    if what == "n":
        return f"Estimated sample size{'s' if M.two_sample else ''} for {M.title}"
    if what == "power":
        return f"Estimated power for {M.title}"
    target = getattr(M, "target_text", None) or {
        "ma": "target mean", "m2": "experimental-group mean", "ma2": "target mean difference",
        "pa": "target proportion", "p2": "experimental-group proportion",
        "ra": "target correlation", "r2": "experimental-group correlation"}[M.effect_name]
    return f"Estimated {target} for {M.title}"


def _header(s: "Session", M: _Method, c: dict, what: str) -> None:
    out = s.output
    if M.iterate and what in ("n", "effect") or (what == "effect"):
        out.write("\nPerforming iteration ...\n", "text")
    out.write("\n" + _what_title(M, what) + "\n", "text")
    out.write(M.test(c) + "\n", "text")
    if _sides(c) == 2:
        op = "!="
    else:
        up = getattr(M, "up", None)
        op = ">" if ((up(c) if up else M.delta(c) > 0) if what != "effect" else _direction(c) > 0) else "<"
    hyp = M.hyp.format(op=op)
    if what == "effect" and _sides(c) == 2:
        # the direction of the computed effect (compat 1301; VERIFY other than onemean)
        e, nul = _EFFECT_PAIR[M.effect_name]
        hyp += f"; {e} {'>' if _direction(c) > 0 else '<'} {nul}"
    out.write(hyp + "\n", "text")


_EFFECT_PAIR = {"ma": ("ma", "m0"), "m2": ("m2", "m1"), "ma2": ("da", "d0"), "pa": ("pa", "p0"),
                "p2": ("p2", "p1"), "ra": ("ra", "r0"), "r2": ("r2", "r1")}


def _line(s: "Session", name: str, value: str) -> None:
    s.output.write(f"{name:>13} = ", "text")
    s.output.write(f"{value}\n", "result")


def _single(s: "Session", M: _Method, c: dict, what: str) -> None:
    out = s.output
    _header(s, M, c, what)
    out.write("\nStudy parameters:\n\n", "text")
    left = [("alpha", _f4(c["alpha"]))]
    if what != "power":
        left.append(("power", _f4(c["power"])))
    if what != "n":
        if M.two_sample:
            left.append(("N", f"{_nfmt(c['N']):>9}"))
            left.append(("N per group" if c["N1"] == c["N2"] else "N1", f"{_nfmt(c['N1']):>9}"))
            if c["N1"] != c["N2"]:
                left.append(("N2", f"{_nfmt(c['N2']):>9}"))
        else:
            left.append(("N", f"{_nfmt(c['N']):>9}"))
    params = M.params(c)
    if what == "effect":
        params = [p for p in params if p[0] not in ("delta", M.effect_name)]
    note = getattr(M, "delta_note", "")
    left += [(name, _f4(val) + (note if name == "delta" else "")) for name, val in params]
    right = [(name, _f4(val)) for name, val in M.right(c)]
    # a second column of parameters beside the first lines (pairedmeans; compat 1301)
    for k, (name, val) in enumerate(left):
        out.write(f"{name:>13} = ", "text")
        out.write(val, "result")
        if k < len(right):
            out.write(f"{right[k][0]:>13} = ", "text")
            out.write(right[k][1], "result")
        out.write("\n", "text")
    if what == "n":
        out.write(f"\nEstimated sample size{'s' if M.two_sample else ''}:\n\n", "text")
        _line(s, "N", f"{_nfmt(c['N']):>9}")
        if M.two_sample:
            if c["N1"] == c["N2"]:
                _line(s, "N per group", f"{_nfmt(c['N1']):>9}")
            else:
                _line(s, "N1", f"{_nfmt(c['N1']):>9}")
                _line(s, "N2", f"{_nfmt(c['N2']):>9}")
    elif what == "power":
        binom = "actual_alpha" in c
        out.write(f"\nEstimated power{' and alpha' if binom else ''}:\n\n", "text")
        _line(s, "power", _f4(c["power_out"]))
        if binom:
            _line(s, "actual alpha", _f4(c["actual_alpha"]))     # compat 1301
    else:
        nm = getattr(M, "target_text", None) or {"ma": "target mean", "m2": "experimental-group mean", "ma2": "target mean",
              "pa": "target proportion", "p2": "experimental-group proportion",
              "ra": "target correlation", "r2": "experimental-group correlation"}[M.effect_name]
        out.write(f"\nEstimated effect size and {nm}:\n\n", "text")
        _line(s, "delta", _f4(M.delta(c)) + getattr(M, "delta_note", ""))
        lab = getattr(M, "label", lambda x: x)
        _line(s, lab(M.effect_name), _f4(c[M.effect_name]) if not getattr(M, "sdscale", False)
              else _f4(math.sqrt(c[M.effect_name])))


def _nfmt(n: float) -> str:
    return f"{n:,.0f}" if float(n).is_integer() else format_value(float(n), "%10.2f", pad=False).strip()


def _table(s: "Session", M: _Method, rows: list[dict], what: str) -> None:
    out = s.output
    _header(s, M, rows[0], what)
    cols = ["alpha", "power"]
    cols += ["N", "N1", "N2"] if M.two_sample else ["N"]
    cols += [name for name, _ in M.params(rows[0])]
    vals = []
    for c in rows:
        d = {"alpha": c["alpha"], "power": c["power_out"] if what == "power" else c["power"],
             "N": c["N"]}
        if M.two_sample:
            d["N1"], d["N2"] = c["N1"], c["N2"]
        d.update(dict(M.params(c)))
        vals.append(d)
    width = [max(7, len(cname)) for cname in cols]
    inner = sum(w + 1 for w in width) + 1
    out.write("\n  +" + "-" * inner + "+\n", "text")
    out.write("  | " + " ".join(f"{cn:>{w}}" for cn, w in zip(cols, width)) + " |\n", "text")
    out.write("  |" + "-" * inner + "|\n", "text")
    for d in vals:
        out.write("  | ", "text")
        out.write(" ".join(f"{_gtab(d.get(cn)):>{w}}" for cn, w in zip(cols, width)), "result")
        out.write(" |\n", "text")
    out.write("  +" + "-" * inner + "+\n", "text")


def _returns(s: "Session", method: str, M: _Method, c: dict, what: str) -> None:
    """r() in the return list order of Stata 14 (compat 1301: onemean and twomeans; VERIFY
    the others)."""
    pw = c["power_out"] if what == "power" else c["power"]     # VERIFY: target or achieved
    solver = []
    if "_init" in c:
        # the solver's report; iter, deltax and function depend on the root finder (VERIFY)
        solver = [("init", c["_init"]), ("converged", 1.0), ("iter", 4.0), ("maxiter", 500.0),
                  ("deltax", 0.0), ("tolerance", 1e-12), ("function", 0.0), ("ftolerance", 1e-12)]
    flags = [("nfractional", 1.0 if c.get("nfractional") else 0.0),
             ("onesided", 1.0 if c.get("onesided") else 0.0)]
    params = dict(M.params(c))
    if method == "onemean":
        r = [("N", float(c["N"])), ("power", pw), ("beta", 1 - pw), ("alpha", c["alpha"])] + flags + solver
        r += [("m0", c["m0"]), ("delta", M.delta(c)), ("knownsd", 1.0 if c.get("knownsd") else 0.0),
              ("sd", c["sd"]), ("diff", c["ma"] - c["m0"]), ("ma", c["ma"])]
    elif method == "twomeans":
        r = [("power", pw), ("beta", 1 - pw), ("alpha", c["alpha"]), ("delta", M.delta(c))] + flags + solver
        r += [("N", float(c["N"])), ("unequal", 0.0 if c["sd1"] == c["sd2"] else 1.0),
              ("knownsds", 1.0 if c.get("knownsds") else 0.0), ("sd2", c["sd2"]), ("sd1", c["sd1"])]
        if c["sd1"] == c["sd2"]:
            r.append(("sd", c["sd1"]))
        r += [("diff", c["m2"] - c["m1"]), ("m2", c["m2"]), ("m1", c["m1"]),
              ("nratio", c["N2"] / c["N1"]), ("N2", float(c["N2"])), ("N1", float(c["N1"]))]
    elif method == "twoproportions":
        r = [("power", pw), ("beta", 1 - pw), ("alpha", c["alpha"]), ("delta", M.delta(c))] + flags + solver
        r += [("N", float(c["N"])), ("continuity", 1.0 if c.get("continuity") else 0.0), ("p2", c["p2"]),
              ("p1", c["p1"]), ("nratio", c["N2"] / c["N1"]), ("N2", float(c["N2"])), ("N1", float(c["N1"])),
              ("test", str(M.o.get("test") or "chi2")), ("effect", str(M.o.get("effect") or "diff"))]
    else:
        r = [("N", float(c["N"]))]
        if M.two_sample:
            r += [("N1", float(c["N1"])), ("N2", float(c["N2"])), ("nratio", c["N2"] / c["N1"])]
            if c["N1"] == c["N2"]:
                r.append(("N_per_group", float(c["N1"])))
        r += [("power", pw), ("beta", 1 - pw), ("alpha", c["alpha"])] + flags + solver
        r.append(("delta", M.delta(c)))
        r += [(k, float(v)) for k, v in params.items() if k != "delta"]
    if what == "effect":
        r.append(("target", float(c[M.effect_name])))
    r += [("type", "test"), ("method", method)]
    # return list shows r() in the reverse order of insertion
    s.r = dict(reversed(r))


@command("power", "power")
def cmd_power(s: "Session", args: str) -> None:
    method, pos, o = _parse(args)
    if method == "oneway":
        power2.cmd_oneway(s, pos, o)
        return
    M = _CLASSES[method](pos, o)
    configs, multi = _configs(method, pos, o)
    rows = []
    what = None
    for c in configs:
        w = _complete(method, c)
        if what is None:
            what = w
        elif w != what:
            raise StataError(198, "the values given ask for different computations")   # VERIFY
        _compute(M, method, c, w)
        rows.append(c)
    show_table = (multi or o.get("table")) and not o.get("notable")
    if show_table:
        _table(s, M, rows, what)
    else:
        _single(s, M, rows[-1], what)
    _returns(s, method, M, rows[-1], what)
    if o.get("graph"):
        power2.graph(s, M, rows, what, o["graph"])


from . import power2  # noqa: E402  (registers the other methods in the tables above)

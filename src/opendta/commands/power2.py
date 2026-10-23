"""power ([PSS] power), more methods: onevariance, twovariances, pairedproportions and oneway, the
likelihood-ratio and Fisher exact tests of twoproportions, the effect options ratio, rrisk and
oratio, and the graph option.

The classes follow the protocol of commands/power.py (power(), params(), set_effect(), ...) and
are registered in its tables when this module is imported. The formulas are those of the
Methods and formulas of [PSS] power onevariance, twovariances, pairedproportions and
twoproportions; the layout of the output is not documented and is marked VERIFY (case
compat/do/1909_other_suites.do)."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core.errors import StataError
from . import power as P
from .power import _Method, _direction, _sides, st

if TYPE_CHECKING:
    from ..session import Session

_SD_LABEL = {"v0": "s0", "va": "sa", "v1": "s1", "v2": "s2"}


# ---------------------------------------------------------------------------
# variances
# ---------------------------------------------------------------------------

class _VarBase(_Method):
    def __init__(self, pos, o):
        super().__init__(pos, o)
        self.sdscale = bool(o.get("sd"))

    def label(self, name: str) -> str:
        return _SD_LABEL.get(name, name) if self.sdscale else name

    def _disp(self, v: float) -> float:
        return math.sqrt(v) if self.sdscale else v

    def test(self, c):
        return "Chi-squared test" if not self.two_sample else "F test"

    @staticmethod
    def _ratio_from_option(c, a: str, b: str) -> None:
        if c.get(b) is None and c.get("ratio") is not None:
            r = c["ratio"] ** 2 if c.get("sdscale") else c["ratio"]
            c[b] = c[a] * r


class _OneVariance(_VarBase):
    title = "a one-sample variance test"
    effect_name = "va"

    def __init__(self, pos, o):
        super().__init__(pos, o)
        if self.sdscale:
            self.title = "a one-sample standard-deviation test"      # case 1909
        v = "s" if self.sdscale else "v"
        self.hyp = f"Ho: {v} = {v}0  versus  Ha: {v} {{op}} {v}0"

    @staticmethod
    def prepare(c):
        _VarBase._ratio_from_option(c, "v0", "va")

    def delta(self, c):
        r = c["va"] / c["v0"]
        return math.sqrt(r) if self.sdscale else r

    def up(self, c):
        return c["va"] > c["v0"]

    def power(self, c):
        n = c["N"]
        if n < 2:
            return 0.0
        df, a = n - 1, c["alpha"]
        r = c["va"] / c["v0"]
        if _sides(c) == 2:
            return float(st.chi2.sf(st.chi2.isf(a / 2, df) / r, df)
                         + st.chi2.cdf(st.chi2.ppf(a / 2, df) / r, df))
        upper = c["va"] > c["v0"] if c.get("va") is not None and not c.get("_effect") else _direction(c) > 0
        if upper:
            return float(st.chi2.sf(st.chi2.isf(a, df) / r, df))
        return float(st.chi2.cdf(st.chi2.ppf(a, df) / r, df))

    def set_effect(self, c, d):
        c["_effect"] = True
        c["va"] = c["v0"] * ((1 + d) if _direction(c) > 0 else 1 / (1 + d)) ** (2 if self.sdscale else 1)

    def params(self, c):
        return [("delta", self.delta(c)), (self.label("v0"), self._disp(c["v0"])),
                (self.label("va"), self._disp(c["va"]))]

    target_text = "target variance"


class _TwoVariances(_VarBase):
    title = "a two-sample variances test"
    two_sample = True
    effect_name = "v2"

    def __init__(self, pos, o):
        super().__init__(pos, o)
        if self.sdscale:
            self.title = "a two-sample standard-deviations test"     # VERIFY
        v = "s" if self.sdscale else "v"
        self.hyp = f"Ho: {v}2 = {v}1  versus  Ha: {v}2 {{op}} {v}1"

    @staticmethod
    def prepare(c):
        _VarBase._ratio_from_option(c, "v1", "v2")

    def delta(self, c):
        r = c["v2"] / c["v1"]
        return math.sqrt(r) if self.sdscale else r

    def up(self, c):
        return c["v2"] > c["v1"]

    def sizes(self, c, n1):
        r = c.get("nratio", 1.0)
        return {"N1": n1, "N2": r * n1, "N": n1 * (1 + r)}

    def power(self, c):
        n1, n2, a = c["N1"], c["N2"], c["alpha"]
        if n1 < 2 or n2 < 2:
            return 0.0
        d1, d2 = n1 - 1, n2 - 1
        r = c["v1"] / c["v2"]            # as in Methods and formulas (matches the manual: nratio(2) gives 225)
        if _sides(c) == 2:
            return float(st.f.sf(r * st.f.isf(a / 2, d1, d2), d1, d2)
                         + st.f.cdf(r * st.f.ppf(a / 2, d1, d2), d1, d2))
        upper = c["v2"] > c["v1"] if not c.get("_effect") else _direction(c) > 0
        if upper:
            return float(st.f.sf(r * st.f.isf(a, d1, d2), d1, d2))
        return float(st.f.cdf(r * st.f.ppf(a, d1, d2), d1, d2))

    def set_effect(self, c, d):
        c["_effect"] = True
        c["v2"] = c["v1"] * ((1 + d) if _direction(c) > 0 else 1 / (1 + d)) ** (2 if self.sdscale else 1)

    def params(self, c):
        rows = [("delta", self.delta(c)), (self.label("v1"), self._disp(c["v1"])),
                (self.label("v2"), self._disp(c["v2"]))]
        if c.get("nratio", 1.0) != 1.0:
            rows.append(("N2/N1", c["nratio"]))
        return rows

    target_text = "experimental-group variance"


# ---------------------------------------------------------------------------
# paired proportions (McNemar)
# ---------------------------------------------------------------------------

class _PairedProportions(_Method):
    title = "a two-sample paired-proportions test"
    hyp = "Ho: p21 = p12  versus  Ha: p21 {op} p12"
    effect_name = "p21"

    def test(self, c):
        return "Large-sample McNemar's test"

    def hypothesis(self, c):
        # with the marginal proportions (corr()), Stata states the test on the margins (case 1909)
        if "pmarg1" in c:
            return "Ho: p+1 = p1+  versus  Ha: p+1 {op} p1+"
        return self.hyp

    @staticmethod
    def prepare(c):
        if c.get("corr") is not None:
            p1, p2 = c["p12"], c["p21"]
            if p2 is None:
                if c.get("diff") is not None:
                    p2 = p1 + c["diff"]
                elif c.get("ratio") is not None:
                    p2 = p1 * c["ratio"]
                elif c.get("oratio") is not None:
                    th = c["oratio"]
                    p2 = th * p1 / (1 - p1 + th * p1)
                else:
                    raise StataError(198, "second marginal proportion or an effect option required")   # VERIFY
            rho = c["corr"]
            p12 = p1 * (1 - p2) - rho * math.sqrt(p1 * (1 - p1) * p2 * (1 - p2))
            c["pmarg1"], c["pmarg2"] = p1, p2
            c["p12"], c["p21"] = p12, p12 + p2 - p1
        else:
            p12 = c["p12"]
            if c.get("p21") is None:
                if c.get("diff") is not None:
                    c["p21"] = p12 + c["diff"]
                elif c.get("ratio") is not None:
                    c["p21"] = p12 * c["ratio"]
                elif c.get("prdiscordant") is not None or c.get("sum") is not None:
                    c["p21"] = (c.get("prdiscordant") if c.get("prdiscordant") is not None else c["sum"]) - p12
        for k in ("diff", "ratio", "oratio", "prdiscordant", "sum"):
            c.pop(k, None)
        if c.get("p21") is not None and not (0 <= c["p12"] and 0 <= c["p21"] and c["p12"] + c["p21"] <= 1):
            raise StataError(198, "invalid discordant proportions")   # VERIFY

    def check_mode(self, what):
        if what == "effect":
            raise StataError(198, "effect-size determination is not implemented for pairedproportions")

    def delta(self, c):
        eff = str(self.o.get("effect") or "diff").lower()
        if eff.startswith("r") and c["p12"] > 0:
            return c["p21"] / c["p12"]
        return c["p21"] - c["p12"]

    def up(self, c):
        return c["p21"] > c["p12"]

    @property
    def delta_note(self):
        eff = str(self.o.get("effect") or "diff").lower()
        return "  (ratio)" if eff.startswith("r") else "  (difference)"     # VERIFY ratio

    def power(self, c):
        n, a = c["N"], c["alpha"]
        pd, pdisc = c["p21"] - c["p12"], c["p21"] + c["p12"]
        v = pdisc - pd * pd
        if n < 1 or v <= 0:
            return 0.0
        d = abs(pd) * math.sqrt(n)
        if _sides(c) == 1:
            return float(st.norm.cdf((d - st.norm.isf(a) * math.sqrt(pdisc)) / math.sqrt(v)))
        z = float(st.norm.isf(a / 2))
        s = math.sqrt(v)
        return float(st.norm.cdf((d - z * math.sqrt(pdisc)) / s) + st.norm.cdf((-d - z * math.sqrt(pdisc)) / s))

    def params(self, c):
        if "pmarg1" in c:
            return [("delta", self.delta(c)), ("p1+", c["pmarg1"]), ("p+1", c["pmarg2"]), ("corr", c["corr"])]
        return [("delta", self.delta(c)), ("p12", c["p12"]), ("p21", c["p21"])]


# ---------------------------------------------------------------------------
# twoproportions: likelihood-ratio and Fisher exact tests
# ---------------------------------------------------------------------------

def _xlogy(a: float, b: float) -> float:
    return 0.0 if a <= 0 else a * math.log(a / b)


def lr_power(n1: float, n2: float, p1: float, p2: float, alpha: float, sides: int) -> float:
    """Power of the likelihood-ratio chi-squared test (normal approximation of the signed
    root statistic, [PSS] power twoproportions)."""
    n = n1 + n2
    w1, w2 = n1 / n, n2 / n
    pb = w1 * p1 + w2 * p2
    kl = (w1 * (_xlogy(p1, pb) + _xlogy(1 - p1, 1 - pb)) + w2 * (_xlogy(p2, pb) + _xlogy(1 - p2, 1 - pb)))
    g = math.sqrt(max(2 * n * kl, 0.0))
    if sides == 1:
        return float(st.norm.cdf(g - st.norm.isf(alpha)))
    z = float(st.norm.isf(alpha / 2))
    return float(st.norm.cdf(g - z) + st.norm.cdf(-g - z))


def fisher_power(n1: int, n2: int, p1: float, p2: float, alpha: float, sides: int,
                 upper: bool) -> tuple[float, float]:
    """Exact power and actual size of Fisher's exact test: the rejection region is built from the
    conditional (hypergeometric) null distribution of the successes of group 2 given the total,
    alpha/2 per tail when two-sided; the power adds the unconditional probabilities of the
    tables inside the region."""
    N = n1 + n2
    f1 = st.binom.pmf(np.arange(n1 + 1), n1, p1)
    f2 = st.binom.pmf(np.arange(n2 + 1), n2, p2)
    null1 = st.binom.pmf(np.arange(n1 + 1), n1, p1)
    null2 = st.binom.pmf(np.arange(n2 + 1), n2, p1)
    pw = sz = 0.0
    tail = alpha / 2 if sides == 2 else alpha
    for m in range(1, N):
        ks = np.arange(max(0, m - n1), min(m, n2) + 1)
        h = st.hypergeom.pmf(ks, N, n2, m)
        sf = np.cumsum(h[::-1])[::-1]           # P(X >= k)
        cdf = np.cumsum(h)                      # P(X <= k)
        rej = np.zeros(len(ks), bool)
        if sides == 2 or upper:
            rej |= sf <= tail + 1e-12
        if sides == 2 or not upper:
            rej |= cdf <= tail + 1e-12
        if not rej.any():
            continue
        k = ks[rej]
        pw += float(np.sum(f1[m - k] * f2[k]))
        sz += float(np.sum(null1[m - k] * null2[k]))
    return min(pw, 1.0), min(sz, 1.0)


def _tp_test(self) -> str:
    t = str(self.o.get("test") or "chi2").lower()
    for k in ("chi2", "lrchi2", "fisher"):
        if k.startswith(t) and t:
            return k
    raise StataError(198, f"test({t}) not allowed")


def _tp_prepare(c):
    if c.get("p2") is None:
        if c.get("rdiff") is not None:
            c["p2"] = c["p1"] + c.pop("rdiff")
        elif c.get("ratio") is not None:
            c["p2"] = c["p1"] * c["ratio"]
        elif c.get("rrisk") is not None:
            c["p2"] = c["p1"] * c["rrisk"]
        elif c.get("oratio") is not None:
            th = c["oratio"]
            c["p2"] = th * c["p1"] / (1 - c["p1"] + th * c["p1"])
        if c.get("p2") is not None and not 0 < c["p2"] < 1:
            raise StataError(198, "invalid experimental-group proportion")   # VERIFY
    for k in ("ratio", "rrisk", "oratio", "rdiff"):
        c.pop(k, None)


def _tp_delta(self, c):
    eff = str(self.o.get("effect") or "diff").lower()
    p1, p2 = c["p1"], c["p2"]
    if eff in ("ratio", "rrisk"):
        return p2 / p1
    if eff == "oratio":
        return (p2 * (1 - p1)) / (p1 * (1 - p2))
    return p2 - p1


_NOTE = {"diff": "  (difference)", "rdiff": "  (risk difference)", "ratio": "  (ratio)",
         "rrisk": "  (relative risk)", "oratio": "  (odds ratio)"}


def _tp_note(self):
    return _NOTE.get(str(self.o.get("effect") or "diff").lower(), "  (difference)")


def _tp_check(self, what):
    if _tp_test(self) == "fisher" and what != "power":
        raise StataError(198, "Fisher's exact method is not allowed with sample-size determination or "
                              "effect-size determination")
    if self.o.get("continuity") and _tp_test(self) != "chi2":
        raise StataError(198, "option continuity is not allowed with test(" + _tp_test(self) + ")")


def _install() -> None:
    TP = P._TwoProportions
    TP.prepare = staticmethod(_tp_prepare)
    TP.delta = _tp_delta
    TP.check_mode = _tp_check
    TP.delta_note = property(_tp_note)
    chi2_power = TP.power

    def power(self, c):
        kind = _tp_test(self)
        if kind == "chi2":
            return chi2_power(self, c)
        n1, n2, p1, p2, a = c["N1"], c["N2"], c["p1"], c["p2"], c["alpha"]
        if kind == "lrchi2":
            return lr_power(n1, n2, p1, p2, a, _sides(c))
        pw, sz = fisher_power(int(round(n1)), int(round(n2)), p1, p2, a, _sides(c), p2 > p1)
        c["actual_alpha"] = sz
        return pw
    TP.power = power

    def test(self, c):
        return {"chi2": "Pearson's chi-squared test", "lrchi2": "Likelihood-ratio test",
                "fisher": "Fisher's exact test"}[_tp_test(self)]
    TP.test = test

    P._METHODS.update({"onevariance": 4, "twovariances": 4, "pairedproportions": 7, "oneway": 4,
                       "pairedmeans": 7})
    P._CLASSES.update({"onevariance": _OneVariance, "twovariances": _TwoVariances,
                       "pairedproportions": _PairedProportions})
    P._POSITIONAL.update({"onevariance": ("v0", "va"), "twovariances": ("v1", "v2"),
                          "pairedproportions": ("p12", "p21")})
    P._ALT_TEXT.update({"onevariance": "alternative variance",
                        "twovariances": "experimental-group variance",
                        "pairedproportions": "second discordant proportion"})
    P._EFFECT_PAIR.update({"va": ("va", "v0"), "v2": ("v2", "v1"), "p21": ("p21", "p12")})
    P._OPTS.update({"rdiff": 2, "rrisk": 2, "oratio": 2, "prdiscordant": 3, "sum": 3,
                    "varerror": 4, "varmeans": 4, "ngroups": 2})


_install()


# ---------------------------------------------------------------------------
# oneway
# ---------------------------------------------------------------------------

def oneway_power(f2: float, g: int, n: float, alpha: float) -> float:
    """Power of the F test of g equal means with n subjects per group; f2 = Var_m / Var_e."""
    N = g * n
    d1, d2 = g - 1, N - g
    if d2 <= 0:
        return 0.0
    crit = float(st.f.isf(alpha, d1, d2))
    return float(st.ncf.sf(crit, d1, d2, N * f2))


def cmd_oneway(s: "Session", pos, o) -> None:
    def one(name, default=None):
        v = P._num(o, name)
        if v is None:
            return default
        if len(v) != 1:
            raise StataError(198, f"{name}() takes a single value for oneway")   # VERIFY
        return v[0]
    means = [p[0] for p in pos if len(p) == 1]
    if len(means) != len(pos):
        raise StataError(198, "group means must be single values for oneway")   # VERIFY
    var_e = one("varerror", 1.0)
    if means:
        g = len(means)
        mb = sum(means) / g
        var_m = sum((m - mb) ** 2 for m in means) / g
    else:
        var_m = one("varmeans")
        g = int(one("ngroups", 0) or 0)
        if var_m is None or g < 2:
            raise StataError(198, "group means, or varmeans() and ngroups(), required")   # VERIFY
    if g < 2:
        raise StataError(198, "at least two groups required")
    alpha = one("alpha", 0.05)
    n_tot, pw = one("n"), one("power")
    if o.get("beta") is not None:
        pw = 1 - one("beta")
    f2 = var_m / var_e
    out = s.output
    if n_tot is not None:
        what = "power"
        n = n_tot / g
        power = oneway_power(f2, g, n, alpha)
    else:
        what = "n"
        pw = 0.8 if pw is None else pw
        n = P._solve(lambda x: oneway_power(f2, g, x, alpha) - pw, 2.0, 100.0)
        if not o.get("nfractional"):
            n = math.ceil(n - 1e-9)
        power = oneway_power(f2, g, n, alpha)
    out.write("\nPerforming iteration ...\n" if what == "n" else "", "text")
    out.write(f"\nEstimated {'sample size' if what == 'n' else 'power'} for one-way ANOVA\n", "text")
    out.write("F test for group effect\n", "text")
    out.write("Ho: delta = 0  versus  Ha: delta != 0\n", "text")
    out.write("\nStudy parameters:\n\n", "text")
    P._line(s, "alpha", P._f4(alpha))
    if what == "n":
        P._line(s, "power", P._f4(pw))
    else:
        P._line(s, "N", f"{P._nfmt(n_tot):>9}")
        P._line(s, "N per group", f"{P._nfmt(n):>9}")
    P._line(s, "delta", P._f4(math.sqrt(f2)))
    P._line(s, "N_g", f"{g:>9}")
    for name, val in [(f"m{i}", m) for i, m in enumerate(means, 1)] + [("Var_m", var_m), ("Var_e", var_e)]:
        P._line(s, name, P._f4(val))
    if what == "n":
        out.write("\nEstimated sample sizes:\n\n", "text")
        P._line(s, "N", f"{P._nfmt(n * g):>9}")
        P._line(s, "N per group", f"{P._nfmt(n):>9}")
    else:
        out.write("\nEstimated power:\n\n", "text")
        P._line(s, "power", P._f4(power))
    s.r = dict(reversed([("N", float(n * g)), ("N_per_group", float(n)), ("power", power),
                         ("beta", 1 - power), ("alpha", alpha), ("delta", math.sqrt(f2)),
                         ("N_g", float(g)), ("Var_m", var_m), ("Var_e", var_e),
                         ("type", "test"), ("method", "oneway")]))


# ---------------------------------------------------------------------------
# graph
# ---------------------------------------------------------------------------

def graph(s: "Session", M: _Method, rows: list[dict], what: str, opt) -> None:
    """The power (or N, or the effect) against the quantity that varies: N when power is
    computed, power otherwise, or the first parameter with several values."""
    from ..core.dataset import Variable
    y_name = {"n": "N", "power": "power"}.get(what, "delta")
    x_name = "N" if what == "power" else "power"
    cols = []
    for c in rows:
        d = {"alpha": c["alpha"], "power": c["power_out"] if what == "power" else c["power"],
             "N": c["N"]}
        d.update(dict(M.params(c)))
        cols.append(d)
    if len({d[x_name] for d in cols}) < 2:
        x_name = next((k for k in cols[0] if k not in (y_name,) and len({d[k] for d in cols}) > 1), None)
    if x_name is None:
        raise StataError(198, "graph requires more than one value in a numlist")   # VERIFY
    order = sorted(range(len(cols)), key=lambda i: cols[i][x_name])
    xs = np.array([cols[i][x_name] for i in order], float)
    ys = np.array([cols[i][y_name] for i in order], float)
    s.run("preserve")
    try:
        s.run("quietly clear")
        s.run(f"quietly set obs {len(xs)}")
        s.data.add(Variable("x", "double", xs))
        s.data.add(Variable("y", "double", ys))
        user = opt if isinstance(opt, str) else ""
        title = "Power" if y_name == "power" else {"N": "Sample size", "delta": "Effect size"}[y_name]
        s.run(f'twoway line y x, ytitle("{y_name}") xtitle("{x_name}") title("{title}") {user}')
    finally:
        s.run("restore")

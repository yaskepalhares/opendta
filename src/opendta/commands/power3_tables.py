"""power cmh, power mcc and power trend ([PSS] power cmh, power mcc, power trend): stratified and
matched 2x2 tables and the trend test. The formulas are in core/power_suites.py."""

from __future__ import annotations

import math
import re

import numpy as np

from ..core import power_suites as C
from ..core.errors import StataError
from . import power as P
from .matrix import get_matrix
from .power3_base import _I, _M, _N, _SOLVER, _fail, _solver_vals, _store


def _probspec(s, body: str):
    """Probabilities per group: numbers or (numlist) per stratum, or a matrix (one row per
    scenario). Returns (list of value lists per column, or list of rows of a matrix)."""
    toks = P._split_args(body)
    if len(toks) == 1 and re.fullmatch(r"[A-Za-z_]\w*", toks[0]):
        m = np.array(get_matrix(s, toks[0]).data, dtype=float)
        return None, [list(r) for r in m]
    return [P._values(t) for t in toks], None


class CMH(_M):
    method = "cmh"
    allowed = {"oratio", "nperstratum", "strweights", "grratios", "continuity", "direction",
               "onesided", "showgrstrsizes"}

    def __init__(self, s, body, o):
        super().__init__(s, body, o)
        self.cols, self.rows = _probspec(s, body)
        self.iterates = not self.flag("onesided")       # Stata 14 prints no iteration message when onesided
        self.K = len(self.cols) if self.cols is not None else len(self.rows[0])
        if self.K < 2:
            _fail("at least two probabilities are required")
        self.sw = self.gr = None
        if o.get("strweights") not in (None, True):
            w = P._values(str(o["strweights"]))
            if len(w) != self.K or min(w) <= 0 or (not self.flag("nfractional") and
                                                      not np.allclose(w, np.round(w))):
                _fail("option strweights() incorrectly specified")
            self.sw = C.standard_weights(w)
        if o.get("grratios") not in (None, True):
            g = P._values(str(o["grratios"]))
            if len(g) != self.K or min(g) <= 0 or max(g) >= 1:
                _fail("option grratios() incorrectly specified")
            self.gr = g
        if self.sw is not None and o.get("nperstratum"):
            _fail("nperstratum() and strweights() may not be combined")
        self.up_dir = not str(o.get("direction") or "upper").lower().startswith("l")

    def base(self):
        b = self.common()
        for k in ("oratio", "nperstratum"):
            v = self.vals(k)
            if v is not None:
                b[k] = v
        if self.cols is not None:
            for j, col in enumerate(self.cols):
                b[f"p1_{j + 1}"] = col
        else:
            b["_row"] = list(range(len(self.rows)))
        return b

    def title(self, what):
        what_ = {"n": "sample size", "power": "power", "effect": "odds ratio"}[what]
        return f"Estimated {what_} for a test of independence in stratified 2x2 tables"

    def hyp_lines(self, c, what):
        op = "!=" if self.sides() == 2 else (">" if (self.up_dir if what == "effect" else c["oratio"] > 1) else "<")
        extra = ", OR " + (">" if self.up_dir else "<") + " 1" if (what == "effect" and self.sides() == 2) else ""
        eq = " = ".join(["OR"] + [f"OR_{j + 1}" for j in range(self.K)])
        cmh = "Cochran-Mantel-Haenszel test" + (" with continuity correction" if self.flag("continuity") else "")
        return [cmh, f"Ho: OR = 1  versus  Ha: OR {op} 1{extra};  {eq}"]

    def p1(self, c):
        if self.cols is not None:
            return np.array([c[f"p1_{j + 1}"] for j in range(self.K)])
        return np.array(self.rows[c["_row"]])

    def compute(self, c):
        self.fix(c)
        K, a, sides, cont = self.K, c["alpha"], self.sides(), self.flag("continuity")
        p1 = self.p1(c)
        if ((p1 <= 0) | (p1 >= 1)).any():
            _fail("probabilities must be between 0 and 1")
        w = self.sw if self.sw is not None else [1.0] * K
        s_ = self.gr if self.gr is not None else [0.5] * K
        sw = sum(w)
        have_n = c.get("n") is not None or c.get("nperstratum") is not None
        theta = c.get("oratio")
        if theta is None:
            if not have_n:
                _fail("oratio() required")                   # VERIFY
            what = "effect"
            c.setdefault("power", 0.8)
        elif have_n:
            if c.get("power") is not None:
                _fail("only two of n(), power() and the effect size may be specified")
            what = "power"
        else:
            c.setdefault("power", 0.8)
            what = "n"

        def power_at(n1, n2, th, up=None):
            return C.cmh_power(n1, n2, p1, C.odds_probability(p1, th), a, sides, cont, up)

        def sizes_total(N):
            nk = np.array(w) * N / sw
            if np.allclose(s_, 0.5):
                return nk / 2, nk / 2
            n2 = nk * np.array(s_)
            return nk - n2, n2
        if what == "n":
            def pw(k):
                n1, n2 = C.cmh_sizes(w, s_, k, self.flag("nfractional"))
                return power_at(n1, n2, theta)
            k = self.size(pw, 1.0 / max(w), c["power"], exact=lambda k: power_at(*C.cmh_sizes(w, s_, k, True), theta))
            if self.exact is not None:
                c["N_a"] = self.exact * sw
            n1, n2 = C.cmh_sizes(w, s_, k, self.flag("nfractional"))
            c["power_out"] = pw(k)
            try:        # r(init) of Stata 14: the root when only the tail in the direction of the effect counts
                c["_init"] = sw * C.solve_up(lambda x: C.cmh_power(*C.cmh_sizes(w, s_, x, True), p1,
                                                                   C.odds_probability(p1, theta), a, sides, cont,
                                                                   tail=True) - c["power"], 1.0 / max(w), 4.0)
            except ValueError:
                c["_init"] = 1.0
        else:
            N = float(c["n"]) if c.get("n") is not None else K * c["nperstratum"]
            n1, n2 = sizes_total(N)
            if what == "power":
                c["power_out"] = power_at(n1, n2, theta)
            else:
                sgn = 1.0 if self.up_dir else -1.0
                try:
                    t = C.solve_up(lambda t: power_at(n1, n2, math.exp(sgn * t), self.up_dir) - c["power"],
                                   1e-9, 1.0)
                except ValueError:
                    raise StataError(498, "failed to find a solution")
                theta = math.exp(sgn * t)
                c["power_out"], c["_init"] = c["power"], 1.0
        c["oratio"], c["n1k"], c["n2k"], c["p1"] = theta, np.array(n1), np.array(n2), p1
        c["N"] = float(np.sum(n1) + np.sum(n2))
        c["delta"] = theta
        nk = c["n1k"] + c["n2k"]
        c["eq_strata"] = bool(np.allclose(nk, nk[0]))
        c["eq_groups"] = bool(abs(c["n1k"].sum() - c["n2k"].sum()) < 1e-9)
        c["balanced"] = c["eq_strata"] and bool(np.allclose(c["n1k"], c["n2k"]))
        return what

    def sizes(self, c):
        items = [("N", _I(c["N"]))]
        nk = c["n1k"] + c["n2k"]
        if c["eq_strata"]:
            items.append(("N per stratum", _I(nk[0])))
        else:
            items += [(f"N{j + 1}", _I(v)) for j, v in enumerate(nk)]
        if c["eq_groups"]:
            items.append(("N per group", _N(c["n1k"].sum())))
        else:
            items += [("N group 1", _I(c["n1k"].sum())), ("N group 2", _I(c["n2k"].sum()))]
        if c["balanced"]:
            items.append(("N per group/stratum", _N(c["n1k"][0])))
        return items

    def show(self, c, what):
        items = [("alpha", P._f4(c["alpha"]))]
        if what != "power":
            items.append(("power", P._f4(c["power"])))
        if what != "n":
            items += self.sizes(c)
        if what != "effect":
            items.append(("delta", P._f4(c["delta"])))
        items.append(("N_s", _I(self.K)))
        items += [(f"p1_{j + 1}", P._f4(v)) for j, v in enumerate(c["p1"])]
        if what != "effect":
            items.append(("oratio", P._f4(c["oratio"])))
        if what == "n":
            tail = ("Estimated sample sizes:", self.sizes(c))
        elif what == "power":
            tail = ("Estimated power:", [("power", P._f4(c["power_out"]))])
        else:
            tail = ("Estimated effect size and odds ratio:",
                    [("delta", P._f4(c["delta"])), ("odds ratio", P._f4(c["oratio"]))])
        self.fit_width(items + tail[1])
        self.put(items)
        self.results(*tail)
        if not c["balanced"] and what == "n":
            self.grstr(c)
        self.returns(c, what)

    def grstr(self, c):
        out = self.s.output
        out.write("\nGroup-per-stratum sample sizes:\n", "text")
        out.write(" " * 15 + "| stratum\n", "text")
        out.write(" " * 15 + "|" + "".join(f"{j + 1:>{9 if j == 0 else 10}}" for j in range(self.K)) + "\n", "text")
        out.write("  " + "-" * 13 + "+" + "-" * (10 * self.K) + "\n", "text")
        out.write("  group".ljust(15) + "|\n", "text")
        for g, arr in ((1, c["n1k"]), (2, c["n2k"])):
            out.write(f"{g:>14} |", "text")
            out.write("".join(f"{P._nfmt(v):>{9 if j == 0 else 10}}" for j, v in enumerate(arr)) + "\n", "result")

    def columns(self, c, what):
        cols = ["alpha", "power", "N", "N_per_stratum" if c["eq_strata"] else "N_avg"]
        if c["balanced"]:
            cols += ["N_per_group", "N_per_grstr"]
        return cols + ["delta", "N_s"] + [f"p1_{j + 1}" for j in range(self.K)] + ["oratio"]

    def row(self, c, what):
        d = {"alpha": c["alpha"], "power": c["power_out"] if what == "power" else c["power"],
             "N": c["N"], "N_per_stratum": c["N"] / self.K, "N_avg": c["N"] / self.K,
             "N_per_group": c["n1k"].sum(), "N_per_grstr": c["n1k"][0],
             "delta": c["delta"], "N_s": self.K, "oratio": c["oratio"]}
        d.update({f"p1_{j + 1}": v for j, v in enumerate(c["p1"])})
        return d

    _ORDER = (["oratio"] + [f"grratio{j}" for j in range(1, 4)] + ["continuity", "N_per_group", "N_per_grstr",
              "balanced", "N_per_stratum"] + _SOLVER + ["onesided", "nfractional", "delta", "alpha", "beta",
              "power", "N"] + [f"p1_{j}" for j in (3, 2, 1)] + ["N_a", "N_s", "N_avg", "type", "method"])

    def returns(self, c, what):
        pw = c["power_out"] if what == "power" else c["power"]
        N = c["N"]
        v = {"N": N, "N_avg": N / self.K, "N_a": c.get("N_a", N), "power": pw, "beta": 1 - pw,
             "alpha": c["alpha"], "delta": c["delta"], "nfractional": float(self.flag("nfractional")),
             "balanced": float(c["balanced"]), "onesided": float(self.sides() == 1),
             "N_per_group": c["n1k"].sum(), "N_s": float(self.K), "oratio": float(c["oratio"]),
             "continuity": float(self.flag("continuity"))}
        if c["eq_strata"]:
            v["N_per_stratum"] = N / self.K
        if c["balanced"]:
            v["N_per_grstr"] = c["n1k"][0]
        s_ = self.gr if self.gr is not None else [0.5] * self.K
        v.update({f"grratio{j + 1}": float(x) for j, x in enumerate(s_)})
        v.update({f"p1_{j + 1}": float(x) for j, x in enumerate(c["p1"])})
        v.update(_solver_vals(c))
        v.update({"type": "test", "method": "cmh"})
        _store(self.s, {k: (float(x) if not isinstance(x, str) else x) for k, x in v.items()}, self._ORDER)


class MCC(_M):
    method = "mcc"
    allowed = {"oratio", "m", "compare", "corr", "direction", "onesided"}

    def __init__(self, s, body, o):
        super().__init__(s, body, o)
        toks = P._split_args(body)
        if len(toks) != 1:
            _fail("p0 required")
        self.p0 = P._values(toks[0])
        self.up_dir = not str(o.get("direction") or "upper").lower().startswith("l")

    def base(self):
        b = self.common()
        b["p0"] = self.p0
        for k in ("oratio", "m", "corr"):
            v = self.vals(k)
            if v is not None:
                b[k] = v
        return b

    def title(self, what):
        w = {"n": "sample size", "power": "power", "effect": "odds ratio"}[what]
        return f"Estimated {w} for a matched case-control study"

    def hyp_lines(self, c, what):
        up = self.up_dir if what == "effect" else c["oratio"] > 1
        return [f"Asymptotic z test, 1:{int(c['M'])} matched design",
                f"Ho: OR = 1  versus  Ha: OR {self.hyp_op(c, up)} 1"]

    def compute(self, c):
        self.fix(c)
        a, sides = c["alpha"], self.sides()
        p0, rho, M = c["p0"], c.get("corr", 0.0), int(c.get("m", 1))
        if not 0 < p0 < 1 or M < 1 or not -1 <= rho <= 1:
            _fail("invalid argument value")                         # VERIFY
        c["M"], c["rho"] = M, rho
        theta = c.get("oratio")
        have_n = c.get("n") is not None
        if theta is None:
            if not have_n:
                _fail("oratio() required")
            what = "effect"
            c.setdefault("power", 0.8)
        elif have_n:
            if c.get("power") is not None:
                _fail("only two of n(), power() and the effect size may be specified")
            what = "power"
        else:
            c.setdefault("power", 0.8)
            what = "n"
        if self.flag("compare") and (what != "n" or M < 2):
            _fail("compare requires sample-size determination and m() of at least 2")   # VERIFY
        try:
            if what == "power":
                c["N"] = float(c["n"])
                c["power_out"] = C.mcc_power(c["N"], p0, theta, rho, M, a, sides)
            elif what == "n":
                c["N"] = self.size(lambda n: C.mcc_power(n, p0, theta, rho, M, a, sides), 1.0, c["power"])
                c["power_out"] = C.mcc_power(c["N"], p0, theta, rho, M, a, sides)
                c["_init"] = C.mcc_size_start(p0, theta, rho, M, a, sides, c["power"])
                if self.flag("compare"):
                    n1 = self.size(lambda n: C.mcc_power(n, p0, theta, rho, 1, a, sides), 1.0, c["power"])
                    c["F_M"] = c["N"] / n1                  # n_M / n_1 with the same corr() (checked, case 1919)
            else:
                sgn = 1.0 if self.up_dir else -1.0
                N = float(c["n"])
                t = C.solve_up(lambda t: C.mcc_power(N, p0, math.exp(sgn * t), rho, M, a, sides, self.up_dir)
                               - c["power"], 1e-9, 1.0)
                theta = math.exp(sgn * t)
                c["N"], c["power_out"], c["_init"] = N, c["power"], 1.0
        except ValueError as e:
            raise StataError(498, f"failed to find a solution: {e}")
        c["oratio"], c["delta"] = theta, theta               # delta is the odds ratio (checked, case 1919)
        c["N"] = float(c["N"])
        return what

    def show(self, c, what):
        items = [("alpha", P._f4(c["alpha"]))]
        if what != "power":
            items.append(("power", P._f4(c["power"])))
        if what != "n":
            items.append(("N cases", _I(c["N"])))
        if what != "effect":
            items += [("delta", P._f4(c["delta"])), ("p0", P._f4(c["p0"])), ("oratio", P._f4(c["oratio"]))]
        else:
            items.append(("p0", P._f4(c["p0"])))
        items += [("corr", P._f4(c["rho"])), ("M", _I(c["M"]))]
        self.put(items)
        if what == "n":
            res = [("N cases", _I(c["N"]))]
            if "F_M" in c:
                res.append(("F_M", P._f4(c["F_M"])))
            self.results("Estimated sample size:", res)
        elif what == "power":
            self.results("Estimated power:", [("power", P._f4(c["power_out"]))])
        else:
            self.results("Estimated effect size and odds ratio:",
                         [("delta", P._f4(c["delta"])), ("odds ratio", P._f4(c["oratio"]))])
        self.returns(c, what)

    def columns(self, c, what):
        return ["alpha", "power", "N", "delta", "M", "p0", "oratio", "corr"]

    def row(self, c, what):
        return {"alpha": c["alpha"], "power": c["power_out"] if what == "power" else c["power"],
                "N": c["N"], "delta": c["delta"], "M": c["M"], "p0": c["p0"], "oratio": c["oratio"],
                "corr": c["rho"]}

    # Stata lists the solver of mcc in the order of the iterations, unlike the other methods
    _ORDER = (["N", "power", "beta", "alpha", "delta", "nfractional", "onesided"] + _SOLVER[::-1]
              + ["p0", "corr", "oratio", "M", "p1", "F_M", "type", "method"])

    def returns(self, c, what):
        pw = c["power_out"] if what == "power" else c["power"]
        v = {"N": c["N"], "power": pw, "beta": 1 - pw, "alpha": c["alpha"], "delta": c["delta"],
             "nfractional": float(self.flag("nfractional")), "onesided": float(self.sides() == 1),
             "p0": c["p0"], "corr": c["rho"], "oratio": c["oratio"], "M": float(c["M"]),
             "p1": C.mcc_p1(c["p0"], c["oratio"], c["rho"])}
        if "F_M" in c:
            v["F_M"] = c["F_M"]
        v.update(_solver_vals(c))
        v.update({"type": "test", "method": "mcc"})
        _store(self.s, {k: (float(x) if not isinstance(x, str) else x) for k, x in v.items()}, self._ORDER)


class Trend(_M):
    method = "trend"
    allowed = {"npergroup", "grweights", "exposure", "continuity", "onesided"}

    def __init__(self, s, body, o):
        super().__init__(s, body, o)
        self.cols, self.rows = _probspec(s, body)
        self.J = len(self.cols) if self.cols is not None else len(self.rows[0])
        if self.J < 2:
            _fail("at least two probabilities are required")
        self.x = None
        if o.get("exposure") not in (None, True):
            x = P._values(str(o["exposure"]))
            if len(x) != self.J or (np.diff(x) <= 0).any():
                _fail("option exposure() incorrectly specified")
            self.x = x
        xs = self.x or list(range(1, self.J + 1))
        self.xs = xs
        if self.flag("continuity") and not np.allclose(np.diff(xs), xs[1] - xs[0]):
            _fail("continuity requires equally spaced exposure levels")
        self.gw = None
        if o.get("grweights") not in (None, True):
            w = P._values(str(o["grweights"]))
            if len(w) != self.J or min(w) <= 0 or (not self.flag("nfractional") and
                                                      not np.allclose(w, np.round(w))):
                _fail("option grweights() incorrectly specified")
            self.gw = C.standard_weights(w)
        if self.gw is not None and o.get("npergroup"):
            _fail("npergroup() and grweights() may not be combined")

    def base(self):
        b = self.common()
        v = self.vals("npergroup")
        if v is not None:
            b["npergroup"] = v
        if self.cols is not None:
            for j, col in enumerate(self.cols):
                b[f"p{j + 1}"] = col
        else:
            b["_row"] = list(range(len(self.rows)))
        return b

    def title(self, what):
        return f"Estimated {'sample size' if what == 'n' else 'power'} for a trend test"

    def hyp_lines(self, c, what):
        up = c["p"][-1] > c["p"][0]
        op = self.hyp_op(c, up)
        cat = "Cochran-Armitage trend test" + (" with continuity correction" if self.flag("continuity") else "")
        return [cat, f"Ho: b = 0  versus  Ha: b {op} 0;  logit(p) = a + b*x"]

    def compute(self, c):
        self.fix(c)
        J, a, sides, cont = self.J, c["alpha"], self.sides(), self.flag("continuity")
        p = (np.array([c[f"p{j + 1}"] for j in range(J)]) if self.cols is not None
             else np.array(self.rows[c["_row"]]))
        if ((p <= 0) | (p >= 1)).any():
            _fail("probabilities must be between 0 and 1")
        c["p"] = p
        c["mono"] = C.is_monotonic(p)
        w = self.gw if self.gw is not None else [1.0] * J
        sw = sum(w)
        have_n = c.get("n") is not None or c.get("npergroup") is not None
        if have_n:
            if c.get("power") is not None:
                _fail("only two of n(), power() and the effect size may be specified")
            N = float(c["n"]) if c.get("n") is not None else float(J * c["npergroup"])
            nj = np.array(w) * N / sw
            c["power_out"] = C.trend_power(nj, p, self.xs, a, sides, cont)
            what = "power"
        else:
            c.setdefault("power", 0.8)
            k = self.size(lambda k: C.trend_power(np.array(w) * k, p, self.xs, a, sides, cont),
                          1.0 / max(w), c["power"])
            nj = np.array(w) * k
            if self.exact is not None:
                c["N_a"] = self.exact * sw
            N = float(nj.sum())
            c["power_out"] = C.trend_power(nj, p, self.xs, a, sides, cont)
            try:        # r(init) of Stata 14: the root when only the tail in the direction of the effect counts
                c["_init"] = sw * C.solve_up(lambda x: C.trend_power(np.array(w) * x, p, self.xs, a, sides, cont,
                                                                     tail=True) - c["power"], 1.0 / max(w), 4.0)
            except ValueError:
                c["_init"] = 1.0
            what = "n"
        c["N"], c["nj"] = N, nj
        c["balanced"] = bool(np.allclose(nj, nj[0]))
        return what

    def sizes(self, c):
        if c["balanced"]:
            return [("N", _I(c["N"])), ("N per group", _I(c["nj"][0]))]
        return [("N", _I(c["N"])), ("Average N", P._f4(c["N"] / self.J))] + [
            (f"N{j + 1}", _I(v)) for j, v in enumerate(c["nj"])]

    def show(self, c, what):
        out = self.s.output
        items = [("alpha", P._f4(c["alpha"]))]
        if what == "n":
            items.append(("power", P._f4(c["power"])))
        else:
            items += self.sizes(c)
        items.append(("N_g", _I(self.J)))
        out.write("\nStudy parameters:\n\n", "text")
        for name, val in items:
            self.line(name, val)
        # the exposure beside each probability when it is given
        for j, v in enumerate(c["p"]):
            out.write(f"{'p' + str(j + 1):>{self.w}} = ", "text")
            out.write(P._f4(v), "result")
            if self.x is not None:
                out.write(f"{'x' + str(j + 1):>{self.w}} = ", "text")
                out.write(P._f4(self.x[j]), "result")
            out.write("\n", "text")
        if what == "n":
            self.results("Estimated sample sizes:", self.sizes(c))
        else:
            self.results("Estimated power:", [("power", P._f4(c["power_out"]))])
        if not c["mono"]:
            out.write("\nWarning: Alternative probabilities are not monotonic.\n", "text")
        self.returns(c, what)

    def header(self, c, what, lines):
        if self.x is None:
            self.s.output.write("\nnote: exposure levels are assumed to be equally spaced\n", "text")
        super().header(c, what, lines)

    def columns(self, c, what):
        cols = ["alpha", "power", "N", "N_per_group" if c["balanced"] else "N_avg", "N_g"]
        return cols + [f"p{j + 1}" for j in range(self.J)]

    def row(self, c, what):
        d = {"alpha": c["alpha"], "power": c["power_out"] if what == "power" else c["power"], "N": c["N"],
             "N_per_group": c["N"] / self.J,
             "N_avg": c["N"] / self.J, "N_g": self.J}
        d.update({f"p{j + 1}": v for j, v in enumerate(c["p"])})
        return d

    def returns(self, c, what):
        pw = c["power_out"] if what == "power" else c["power"]
        J, N = self.J, c["N"]
        v = {"N": N, "N_per_group": N / J, "N_avg": N / J, "N_a": c.get("N_a", N), "power": pw,
             "beta": 1 - pw, "alpha": c["alpha"], "nfractional": float(self.flag("nfractional")),
             "balanced": float(c["balanced"]), "onesided": float(self.sides() == 1), "N_g": float(J),
             "continuity": float(self.flag("continuity"))}
        v.update({f"p{j + 1}": float(x) for j, x in enumerate(c["p"])})
        if self.x is not None:
            v.update({f"x{j + 1}": float(x) for j, x in enumerate(self.x)})      # VERIFY: place of x_j
        v.update(_solver_vals(c))
        v.update({"type": "test", "method": "trend"})
        order = (["continuity", "balanced", "N_per_group", "N_avg", "N_g", "N_a"]
                 + [f"p{j + 1}" for j in range(J)] + [f"x{j + 1}" for j in range(J)]
                 + ["N", "power", "beta", "alpha"] + _SOLVER + ["onesided", "nfractional", "type", "method"])
        _store(self.s, {k: (float(x) if not isinstance(x, str) else x) for k, x in v.items()}, order)

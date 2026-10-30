"""power ([PSS] power), part 3: twoway, repeated, cmh, mcc, trend, cox, exponential and logrank.

The formulas are in core/power_suites.py (the Methods and formulas of each [PSS] page); this
module reads the arguments and options, chooses what to compute (sample size, power or effect
size), writes the output and stores r(). The commands are added to the method tables of
commands/power.py when this module is imported.

The layout of the output, the order of r() and the rounding of the sizes are only partly given
by the manual and are marked VERIFY (case compat/do/1919_other_suites2.do). Options that are
not implemented raise "not available yet": graph(), saving(), show, st1(), n#()."""

from __future__ import annotations

import itertools
import math
import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import power_suites as C
from ..core.errors import StataError
from . import power as P
from .matrix import Matrix, get_matrix

if TYPE_CHECKING:
    from ..session import Session

METHODS = {"twoway": 4, "repeated": 3, "cmh": 3, "mcc": 3, "trend": 3, "cox": 3,
           "exponential": 3, "logrank": 3}          # minimal abbreviations: VERIFY

_NEW_OPTS = {"npercell": 5, "cellweights": 4, "nrows": 5, "ncols": 5, "factor": 3,
             "vareffect": 4, "varrow": 6, "varcolumn": 9, "varrowcolumn": 12, "varerror": 6,
             "varbetween": 9, "varwithin": 7, "varbwithin": 8, "showcellsizes": 5,
             "showmatrices": 8, "showmeans": 7, "npergroup": 5, "grweights": 3, "ngroups": 4,
             "nrepeated": 4, "covmatrix": 3, "nperstratum": 5, "strweights": 3, "grratios": 3,
             "showgrstrsizes": 9, "showasmatrix": 9, "m": 1, "compare": 4, "exposure": 3,
             "hratio": 2, "lnhratio": 4, "hdifference": 5, "difference": 4, "loghazard": 4,
             "unconditional": 4, "studytime": 4, "fperiod": 3, "aperiod": 3, "aprob": 3,
             "aptime": 3, "atime": 3, "ashape": 3, "lossprob": 8, "lossprob1": 9, "lossprob2": 9,
             "losstime": 5, "losshaz": 7, "losshaz1": 8, "losshaz2": 8, "time": 4, "show": 4,
             "schoenfeld": 3, "simpson": 3, "st1": 3, "wdprob": 3, "r2": 2, "eventprob": 3,
             "failprob": 3, "nodots": 5, "notitle": 6, "hdiff": 5}

_COMMON = {"alpha", "power", "beta", "n", "nfractional", "parallel", "table", "notable", "nolog",
           "log", "dots", "nodots", "init", "iterate", "tolerance", "ftolerance", "notitle"}
_NOT_YET = {"graph", "saving", "show", "st1", "showmatrices", "showmeans", "showasmatrix",
            "n1", "n2", "effect_exp", "gridlines"}

OP = {"!=": "!=", "up": ">", "lo": "<"}


def _I(x) -> str:
    return f"{P._nfmt(x):>9}"


def _N(x) -> str:
    """A size that is fractional when N is odd and the groups are equal (N per group)."""
    return _I(x) if abs(x - round(x)) < 1e-9 else P._f4(x)


def _T(txt: str) -> str:
    return f"{txt:>9}"


def _fail(what: str):
    raise StataError(198, what)


def _grid(text: str) -> list[list[float]]:
    text = text.strip()
    if text.startswith("(") and text.endswith(")"):
        text = text[1:-1]
    rows = []
    for r in text.split("\\"):
        try:
            rows.append([float(x) for x in r.replace(",", " ").split()])
        except ValueError:
            _fail("invalid matrix specification")
    if len({len(r) for r in rows}) != 1 or not rows[0]:
        _fail("invalid matrix specification")
    return rows


def _combos(base: dict, parallel: bool):
    keys = list(base)
    lists = [base[k] for k in keys]
    if parallel:
        lens = {len(v) for v in lists if len(v) > 1}
        if len(lens) > 1:
            _fail("parallel requires lists of the same length")        # VERIFY
        L = max(lens) if lens else 1
        rows = zip(*[v * L if len(v) == 1 else v for v in lists])
    else:
        rows = itertools.product(*lists)
    return [dict(zip(keys, r)) for r in rows], any(len(v) > 1 for v in lists)


# The solver's items in the order Stata lists them (VERIFY: iter, deltax and function depend on
# Stata's root finder; init is its undocumented starting value)
_SOLVER = ["ftolerance", "function", "tolerance", "deltax", "maxiter", "iter", "converged", "init"]


def _MAT(txt: str = "<matrix>") -> str:
    return f"{txt:>10}"


def _solver_vals(c: dict) -> dict:
    if "_init" not in c:
        return {}
    return {"init": float(c["_init"]), "converged": 1.0, "iter": 3.0, "maxiter": 500.0,
            "deltax": 0.0, "tolerance": 1e-12, "function": 0.0, "ftolerance": 1e-12}


def _store(s, vals: dict, order: list) -> None:
    """r() in the order of Stata's `return list` (observed in case 1919: it is not the logical
    order and differs by method); names missing from `order` go last. VERIFY: the order of
    variants of a command that the case does not list."""
    shown = {k: vals[k] for k in order if k in vals}
    shown.update({k: v for k, v in vals.items() if k not in shown})
    s.r = dict(reversed(list(shown.items())))


class _M:
    method = ""
    allowed: set = set()
    iterates = True           # prints "Performing iteration ..." for sample size and effect size
    effect_title = ""

    def __init__(self, s: "Session", body: str, o: dict):
        self.s, self.o, self.body = s, o, body
        for k in o:
            if k in _NOT_YET:
                _fail(f"option {k} is not available yet in OpenDTA")
            if k not in _COMMON and k not in self.allowed:
                _fail(f"option {k} not allowed")
        if o.get("onesided") and "onesided" not in self.allowed:
            _fail("option onesided not allowed")

    w = 13                     # width of the labels; fit_width() widens it to the longest label plus two

    def line(self, name, val) -> None:
        out = self.s.output
        out.write(f"{name:>{self.w}} = ", "text")
        out.write(f"{val}\n", "result")

    def fit_width(self, items: list) -> None:
        """Stata widens the labels to the longest one plus two blanks (case 1919)."""
        self.w = max(13, max(len(n) for n, _ in items) + 2)

    def vals(self, name):
        return P._num(self.o, name)

    def one(self, name, default=None):
        v = self.vals(name)
        if v is None:
            return default
        if len(v) != 1:
            _fail(f"{name}() takes a single value for {self.method}")     # VERIFY
        return v[0]

    def flag(self, name) -> bool:
        return bool(self.o.get(name))

    def common(self) -> dict:
        b = {"alpha": self.vals("alpha") or [0.05]}
        for k in ("power", "beta", "n"):
            v = self.vals(k)
            if v is not None:
                b[k] = v
        return b

    def sides(self) -> int:
        return 1 if self.flag("onesided") else 2

    def fix(self, c: dict) -> None:
        if c.get("beta") is not None:
            c["power"] = 1 - c["beta"]

    def size(self, pw, lo, target, start=None, fractional=None, exact=None):
        """Smallest size (integer unless nfractional) whose power reaches the target. The size
        before rounding (r(N_a)) is left in self.exact, in the same units."""
        frac = self.flag("nfractional") if fractional is None else fractional
        f = exact or pw
        self.exact = None
        try:
            self.exact = float(C.solve_up(lambda x: f(x) - target, lo, max(lo * 2, 4.0)))
        except ValueError:
            pass
        try:
            if frac:
                if self.exact is None:
                    raise ValueError
                return self.exact
            return float(C.smallest_int(pw, target, int(math.ceil(lo))))
        except ValueError:
            raise StataError(498, "failed to find a solution")           # VERIFY

    # -- output ---------------------------------------------------------------------------
    def header(self, c: dict, what: str, lines: list[str]) -> None:
        out = self.s.output
        if self.iterates and what in ("n", "effect") and not self.flag("nolog"):
            out.write("\nPerforming iteration ...\n", "text")
        out.write("\n" + self.title(what) + "\n", "text")
        for ln in lines:
            out.write(ln + "\n", "text")

    def put(self, items: list) -> None:
        out = self.s.output
        out.write("\nStudy parameters:\n\n", "text")
        for name, val in items:
            self.line(name, val)

    def results(self, heading: str, items: list) -> None:
        self.s.output.write(f"\n{heading}\n\n", "text")
        for name, val in items:
            self.line(name, val)

    def table(self, rows: list[dict], what: str) -> None:
        out = self.s.output
        self.header(rows[0], what, self.hyp_lines(rows[0], what))
        if getattr(self, "means", None) is not None:
            out.write("   means = <matrix>\n", "text")
        cols = self.columns(rows[0], what)
        vals = [self.row(c, what) for c in rows]
        width = [max(7, len(cn)) for cn in cols]
        inner = sum(w + 1 for w in width) + 1
        out.write("\n  +" + "-" * inner + "+\n", "text")
        out.write("  | " + " ".join(f"{cn:>{w}}" for cn, w in zip(cols, width)) + " |\n", "text")
        out.write("  |" + "-" * inner + "|\n", "text")
        for d in vals:
            out.write("  | ", "text")
            out.write(" ".join(f"{P._gtab(d.get(cn)):>{w}}" for cn, w in zip(cols, width)), "result")
            out.write(" |\n", "text")
        out.write("  +" + "-" * inner + "+\n", "text")

    def single(self, c: dict, what: str) -> None:
        self.header(c, what, self.hyp_lines(c, what))
        self.show(c, what)

    def hyp_op(self, c: dict, up: bool) -> str:
        return "!=" if self.sides() == 2 else (">" if up else "<")


# ---------------------------------------------------------------------------------------------
# twoway
# ---------------------------------------------------------------------------------------------

class TwoWay(_M):
    method = "twoway"
    allowed = {"npercell", "cellweights", "nrows", "ncols", "factor", "vareffect", "varrow",
               "varcolumn", "varrowcolumn", "varerror", "showcellsizes"}
    _NAMES = {"row": "row", "column": "column", "rowcol": "row-by-column"}
    _VAR = {"row": "Var_r", "column": "Var_c", "rowcol": "Var_rc"}

    def __init__(self, s, body, o):
        super().__init__(s, body, o)
        b = body.strip()
        self.means = None
        if b:
            if re.fullmatch(r"[A-Za-z_]\w*", b):
                self.means = np.array(get_matrix(s, b).data, dtype=float)
            else:
                self.means = np.array(_grid(b), dtype=float)
            if min(self.means.shape) < 2:
                _fail("at least two rows and two columns of means are required")
            self.J, self.K = self.means.shape
        else:
            nr, nc = self.one("nrows"), self.one("ncols")
            if nr is None or nc is None:
                _fail("nrows() and ncols() are required without cell means")
            self.J, self.K = int(nr), int(nc)
            if min(self.J, self.K) < 2:
                _fail("nrows() and ncols() must be at least 2")
        fac = o.get("factor")
        which = None
        if fac not in (None, True):
            f = str(fac).strip().lower()
            which = ("rowcol" if f.startswith("rowc") else "row" if f in ("r", "ro", "row")
                     else "column" if f and "column".startswith(f) else None)
            if which is None:
                _fail("option factor() incorrectly specified")
        named = [n for n in ("varrow", "varcolumn", "varrowcolumn") if o.get(n)]
        if len(named) > 1 or (named and which is not None):
            _fail("only one of factor(), varrow(), varcolumn() or varrowcolumn() may be specified")
        if named:
            which = {"varrow": "row", "varcolumn": "column", "varrowcolumn": "rowcol"}[named[0]]
            self.veff = self.vals(named[0])
        else:
            self.veff = self.vals("vareffect")
        self.which = which or "row"
        self.cw = None
        if o.get("cellweights") not in (None, True):
            t = str(o["cellweights"]).strip()
            arr = (np.array(get_matrix(s, t).data, dtype=float) if re.fullmatch(r"[A-Za-z_]\w*", t)
                   else np.array(_grid(t), dtype=float))
            if arr.shape != (self.J, self.K) or (arr <= 0).any():
                _fail("option cellweights() incorrectly specified")
            if not self.flag("nfractional") and not np.allclose(arr, np.round(arr)):
                _fail("cellweights() must contain integers unless nfractional is specified")
            self.cw = arr
            self.we = np.array(C.standard_weights(list(arr.reshape(-1)))).reshape(arr.shape)
        if self.cw is not None and o.get("npercell"):
            _fail("npercell() and cellweights() may not be combined")

    def base(self):
        b = self.common()
        v = self.vals("npercell")
        if v is not None:
            b["npercell"] = v
        if self.veff is not None:
            b["veff"] = self.veff
        b["var_e"] = self.vals("varerror") or [1.0]
        return b

    def title(self, what):
        if what == "n":
            return "Estimated sample size for two-way ANOVA"
        if what == "power":
            return "Estimated power for two-way ANOVA"
        return f"Estimated {self._NAMES[self.which]} variance for two-way ANOVA"

    def hyp_lines(self, c, what):
        return [f"F test for {self._NAMES[self.which]} effect", "Ho: delta = 0  versus  Ha: delta != 0"]

    def compute(self, c):
        self.fix(c)
        J, K, which, a = self.J, self.K, self.which, c["alpha"]
        ve = c["var_e"]
        if self.means is not None:
            veff = C.twoway_variances(self.means)[("row", "column", "rowcol").index(which)]
        else:
            veff = c.get("veff")
        have_n = c.get("n") is not None or c.get("npercell") is not None
        if self.means is None and veff is None:
            if not have_n:
                _fail("cell means or vareffect() required")
            what = "effect"
            if c.get("power") is None:
                _fail("power() is required to compute the effect size")          # VERIFY
        else:
            # VERIFY: vareffect() given together with cell means is ignored
            if have_n:
                if c.get("power") is not None:
                    _fail("only two of n(), power() and the effect size may be specified")
                what = "power"
            else:
                c.setdefault("power", 0.8)
                what = "n"
        sw_bal = np.full((J, K), 1.0)
        w = self.we if self.cw is not None else sw_bal
        frac = w / w.sum()

        def lam(N, fr):
            if self.means is not None:
                return C.twoway_noncentrality(self.means, fr, N, which, ve)
            return N * veff / ve

        if what == "power":
            N = float(c["n"]) if c.get("n") is not None else J * K * c["npercell"]
            fr = frac if self.cw is not None else sw_bal / (J * K)
            Ne = N
            if self.cw is None and not self.flag("nfractional"):       # whole subjects per cell
                Ne = J * K * math.floor(N / (J * K) + 1e-9)
            c["power_out"] = C.twoway_power(lam(Ne, fr), J, K, Ne, which, a)
            c["N"], c["delta2"], c["Ne"] = N, lam(Ne, fr) / Ne, Ne
        elif what == "n":
            sw = w.sum()

            def pw(k):
                return C.twoway_power(lam(k * sw, frac), J, K, k * sw, which, a)
            k = self.size(pw, math.floor(J * K / sw) + 1, c["power"])
            c["N"], c["delta2"] = k * sw, lam(k * sw, frac) / (k * sw)
            c["power_out"] = pw(k)
            c["_init"] = float(J * K)
            if self.exact is not None:
                c["N_a"] = self.exact * sw
        else:
            if self.cw is not None:
                _fail("option cellweights() is not available yet in OpenDTA for the effect size")
            N = float(c["n"]) if c.get("n") is not None else J * K * c["npercell"]
            c["N"] = N
            lm = C.lambda_for_power(C.twoway_df1(J, K, which), N - J * K, a, c["power"])
            c["delta2"] = lm / N
            c["power_out"] = c["power"]
            c["_init"] = 1.0
        c["veff"] = c["delta2"] * ve if (what == "effect" or self.cw is not None) else veff
        c["balanced"] = self.cw is None
        c["cells"] = (w / w.sum() * c["N"]) if self.cw is not None else None
        c["delta"] = math.sqrt(c["delta2"])
        return what

    def sizes(self, c) -> list:
        N = c["N"]
        if c["balanced"]:
            return [("N", _I(N)), ("N per cell", _I(c.get("Ne", N) / (self.J * self.K)))]
        return [("N", _I(N)), ("Average N", P._f4(N / (self.J * self.K))), ("Nij", _MAT())]

    def show(self, c, what):
        vn = self._VAR[self.which]
        items = [("alpha", P._f4(c["alpha"]))]
        if what != "n":
            if what == "effect":
                items.append(("power", P._f4(c["power"])))
            items += self.sizes(c)
        else:
            items.append(("power", P._f4(c["power"])))
        if what != "effect":
            items.append(("delta", P._f4(c["delta"])))
        items += [("N_r", _I(self.J)), ("N_c", _I(self.K))]
        if self.means is not None:
            items.append(("means", _MAT()))
        if what != "effect":
            items.append((vn, P._f4(c["veff"])))
        items.append(("Var_e", P._f4(c["var_e"])))
        self.put(items)
        if what == "n":
            self.results("Estimated sample sizes:", self.sizes(c))
            if self.flag("showcellsizes") and not c["balanced"]:
                self.cellsizes(c)
        elif what == "power":
            self.results("Estimated power:", [("power", P._f4(c["power_out"]))])
        else:
            self.results(f"Estimated effect size and {self._NAMES[self.which]} variance:",
                         [("delta", P._f4(c["delta"])), (vn, P._f4(c["veff"]))])
        self.returns(c, what)

    def cellsizes(self, c):                       # VERIFY: layout of the cell sizes
        out = self.s.output
        out.write("\nCell sample sizes\n", "text")
        out.write("columns\n" + " ".join(f"{j + 1:>6}" for j in range(self.K)) + "\nrows\n", "text")
        for i in range(self.J):
            out.write(f"{i + 1} " + " ".join(f"{P._nfmt(v):>6}" for v in c["cells"][i]) + "\n", "result")

    def columns(self, c, what):
        cols = ["alpha", "power", "N", "N_per_cell" if c["balanced"] else "N_avg", "delta", "N_r", "N_c"]
        return cols + [self._VAR[self.which], "Var_e"]

    def row(self, c, what):
        N = c["N"]
        return {"alpha": c["alpha"], "power": c["power_out"] if what == "power" else c["power"],
                "N": N, "N_per_cell": c.get("Ne", N) / (self.J * self.K), "N_avg": N / (self.J * self.K),
                "delta": c["delta"], "N_r": self.J, "N_c": self.K, self._VAR[self.which]: c["veff"],
                "Var_e": c["var_e"]}

    _ORDER = (["Var_e", "balanced", "N_per_cell", "N_avg", "N_rc", "N_r", "N_c", "N_a", "m1_1", "m1_2",
               "delta"] + _SOLVER + ["onesided", "nfractional", "alpha", "beta", "power", "N", "m2_3",
                                     "m2_2", "m2_1", "m1_3", "type", "method", "means", "Nij"])

    def returns(self, c, what):
        pw = c["power_out"] if what == "power" else c["power"]
        N, JK = c["N"], self.J * self.K
        v = {"N": float(N), "N_avg": N / JK, "N_rc": float(JK), "N_a": c.get("N_a", float(N)),
             "power": pw, "beta": 1 - pw, "alpha": c["alpha"], "delta": c["delta"],
             "nfractional": float(self.flag("nfractional")), "onesided": 0.0,
             "balanced": float(c["balanced"]), "N_r": float(self.J), "N_c": float(self.K),
             self._VAR[self.which]: float(c["veff"]), "Var_e": float(c["var_e"])}
        if c["balanced"]:
            v["N_per_cell"] = c.get("Ne", N) / JK
        v.update(_solver_vals(c))
        v.update({"type": "test", "method": "twoway"})
        cells = np.full((self.J, self.K), c.get("Ne", N) / JK) if c["balanced"] else c["cells"]
        v["Nij"] = Matrix(cells)
        if self.means is not None:
            v["means"] = Matrix(self.means)
            for i in range(self.J):
                for j in range(self.K):
                    v[f"m{i + 1}_{j + 1}"] = float(self.means[i, j])
        _store(self.s, v, [self._VAR[self.which]] + self._ORDER)


# ---------------------------------------------------------------------------------------------
# repeated
# ---------------------------------------------------------------------------------------------

class Repeated(_M):
    method = "repeated"
    allowed = {"npergroup", "grweights", "ngroups", "nrepeated", "corr", "covmatrix", "factor",
               "vareffect", "varbetween", "varwithin", "varbwithin", "varerror"}
    _VARN = {"between": ("Var_b", "Var_be"), "within": ("Var_w", "Var_we"),
             "bwithin": ("Var_bw", "Var_bwe")}
    _TEST = {"between": "F test for between subjects",
             "within": "F test for within subject with Greenhouse-Geisser correction",
             "bwithin": "F test for between-within subjects with Greenhouse-Geisser correction"}
    _EFF = {"between": "between-subjects", "within": "within-subject",
            "bwithin": "between-within-subjects"}          # VERIFY: names in the effect-size title

    def __init__(self, s, body, o):
        super().__init__(s, body, o)
        b = body.strip()
        self.means = None
        self.cov = None
        if o.get("covmatrix") not in (None, True):
            self.cov = np.array(get_matrix(s, str(o["covmatrix"])).data, dtype=float)
        elif not o.get("corr"):
            _fail("corr() or covmatrix() is required")
        if b:
            self.means = (np.array(get_matrix(s, b).data, dtype=float)
                          if re.fullmatch(r"[A-Za-z_]\w*", b) else np.array(_grid(b), dtype=float))
            self.J, self.K = self.means.shape
            if self.K < 2:
                _fail("at least two repeated measurements are required")
        else:
            ng, nr = self.one("ngroups"), self.one("nrepeated")
            if ng is None:
                _fail("ngroups() is required without cell means")
            self.J = int(ng)
            self.K = int(nr) if nr is not None else (self.cov.shape[0] if self.cov is not None else 0)
            if self.K < 2:
                _fail("nrepeated() is required")
        if self.cov is not None and self.cov.shape != (self.K, self.K):
            _fail("covmatrix() does not conform to the number of repeated measurements")
        fac = o.get("factor")
        which = "between"
        if fac not in (None, True):
            f = str(fac).strip().lower()
            which = ("bwithin" if f.startswith("bw") else "between" if f and "between".startswith(f)
                     else "within" if f and "within".startswith(f) else None)
            if which is None:
                _fail("option factor() incorrectly specified")
        named = [n for n in ("varbetween", "varwithin", "varbwithin") if o.get(n)]
        if len(named) > 1 or (named and fac not in (None, True)):
            _fail("only one of factor(), varbetween(), varwithin() or varbwithin() may be specified")
        if named:
            which = {"varbetween": "between", "varwithin": "within", "varbwithin": "bwithin"}[named[0]]
            self.veff = self.vals(named[0])
        else:
            self.veff = self.vals("vareffect")
        self.which = which
        if self.cov is not None and o.get("varerror"):
            _fail("varerror() may not be combined with covmatrix()")
        self.gw = None
        if o.get("grweights") not in (None, True):
            w = P._values(str(o["grweights"]))
            if len(w) != self.J or min(w) <= 0:
                _fail("option grweights() incorrectly specified")
            self.gw = w
        if self.gw is not None and o.get("npergroup"):
            _fail("npergroup() and grweights() may not be combined")

    def base(self):
        b = self.common()
        v = self.vals("npergroup")
        if v is not None:
            b["npergroup"] = v
        if self.veff is not None:
            b["veff"] = self.veff
        if self.cov is None:
            b["corr"] = self.vals("corr")
            b["var_e"] = self.vals("varerror") or [1.0]
        return b

    def title(self, what):
        if what == "n":
            return "Estimated sample size for repeated-measures ANOVA"
        if what == "power":
            return "Estimated power for repeated-measures ANOVA"
        return f"Estimated {self._EFF[self.which]} variance for repeated-measures ANOVA"

    def hyp_lines(self, c, what):
        return [self._TEST[self.which], "Ho: delta = 0  versus  Ha: delta != 0"]

    def compute(self, c):
        self.fix(c)
        J, K, a = self.J, self.K, c["alpha"]
        cov = self.cov if self.cov is not None else C.cs_covariance(K, c["var_e"], c["corr"])
        ws = C.standard_weights(self.gw) if self.gw is not None else [1.0] * J
        wf = np.array(ws) / sum(ws)
        have_n = c.get("n") is not None or c.get("npergroup") is not None
        have_eff = self.means is not None or c.get("veff") is not None
        try:
            ve_arg = c.get("veff") if c.get("veff") is not None else (1.0 if self.means is None else None)
            sm = C.repeated_summary(J, K, cov, self.which, wf, self.means, ve_arg)
        except ValueError as e:
            _fail(str(e))
        if not have_eff:
            if not have_n:
                _fail("cell means or a variance option required")
            what = "effect"
            if c.get("power") is None:
                _fail("power() is required to compute the effect size")          # VERIFY
        elif have_n:
            if c.get("power") is not None:
                _fail("only two of n(), power() and the effect size may be specified")
            what = "power"
        else:
            c.setdefault("power", 0.8)
            what = "n"
        pw = lambda N, d2: C.repeated_power(d2, sm["dc"], sm["du"], J, N, sm["eps"], a)
        sw = sum(ws)
        if what == "power":
            N = float(c["n"]) if c.get("n") is not None else float(J * c["npergroup"])
            c["N"], c["power_out"], c["delta2"] = N, pw(N, sm["delta2"]), sm["delta2"]
            c["beta_out"] = C.repeated_power(sm["delta2"], sm["dc"], sm["du"], J, N, sm["eps"], a, lower=True)
        elif what == "n":
            k = self.size(lambda k: pw(k * sw, sm["delta2"]), math.floor(J / sw) + 1, c["power"])
            c["N"], c["delta2"] = k * sw, sm["delta2"]
            c["power_out"] = pw(k * sw, sm["delta2"])
            c["_init"] = float(J)
            if self.exact is not None:
                c["N_a"] = self.exact * sw
        else:
            N = float(c["n"]) if c.get("n") is not None else float(J * c["npergroup"])
            try:
                d2 = C.solve_up(lambda d: pw(N, d) - c["power"], 1e-9, 1.0)
            except ValueError:
                raise StataError(498, "failed to find a solution")
            c["N"], c["delta2"], c["power_out"] = N, d2, c["power"]
            c["_init"] = 1.0
        c["sm"], c["cov"] = sm, cov
        c["delta"] = math.sqrt(c["delta2"])
        c["var_x"] = c["delta2"] * sm["var_xe"] if what == "effect" else sm["var_x"]
        c["balanced"] = self.gw is None or len(set(ws)) == 1
        c["groups"] = [c["N"] * x for x in wf]
        return what

    def sizes(self, c):
        if c["balanced"]:
            return [("N", _I(c["N"])), ("N per group", _I(c["N"] / self.J))]
        return [("N", _I(c["N"])), ("Average N", P._f4(c["N"] / self.J))] + [
            (f"N{j + 1}", _I(g)) for j, g in enumerate(c["groups"])]

    def show(self, c, what):
        vx, vxe = self._VARN[self.which]
        items = [("alpha", P._f4(c["alpha"]))]
        if what == "effect":
            items.append(("power", P._f4(c["power"])))
        if what != "n":
            items += self.sizes(c)
        else:
            items.append(("power", P._f4(c["power"])))
        if what != "effect":
            items.append(("delta", P._f4(c["delta"])))
        items += [("N_g", _I(self.J)), ("N_rep", _I(self.K))]
        if self.means is not None:
            items.append(("means", _MAT()))
        if what != "effect":
            items.append((vx, P._f4(c["var_x"])))
        items.append((vxe, P._f4(c["sm"]["var_xe"])))
        if self.cov is not None:
            items += [("Cov", _MAT()), ("spherical", _T("true" if c["sm"]["spherical"] else "false"))]
        else:
            items += [("Var_e", P._f4(c["var_e"])), ("rho", P._f4(c["corr"]))]
        self.put(items)
        if what == "n":
            self.results("Estimated sample sizes:", self.sizes(c))
        elif what == "power":
            self.results("Estimated power:", [("power", P._f4(c["power_out"]))])
        else:
            self.results(f"Estimated effect size and {self._EFF[self.which]} variance:",
                         [("delta", P._f4(c["delta"])), (vx, P._f4(c["var_x"]))])
        self.returns(c, what)

    def columns(self, c, what):
        vx, vxe = self._VARN[self.which]
        cols = ["alpha", "power", "N", "N_per_group" if c["balanced"] else "N_avg", "delta", "N_g", "N_rep", vx]
        return cols + ([vxe] if self.cov is not None else [vxe, "rho"])

    def row(self, c, what):
        vx, vxe = self._VARN[self.which]
        return {"alpha": c["alpha"], "power": c["power_out"] if what == "power" else c["power"],
                "N": c["N"], "N_per_group": c["N"] / self.J, "N_avg": c["N"] / self.J,
                "delta": c["delta"], "N_g": self.J, "N_rep": self.K, vx: c["var_x"],
                vxe: c["sm"]["var_xe"], "rho": c.get("corr")}

    def returns(self, c, what):
        pw = c["power_out"] if what == "power" else c["power"]
        vx, vxe = self._VARN[self.which]
        N = c["N"]
        v = {"N": float(N), "N_per_group": N / self.J, "N_avg": N / self.J, "N_a": c.get("N_a", float(N)),
             "power": pw, "beta": c.get("beta_out", 1 - pw), "alpha": c["alpha"],
             "delta": c["delta"], "nfractional": float(self.flag("nfractional")), "onesided": 0.0,
             "balanced": float(c["balanced"]), "N_g": float(self.J), "N_rep": float(self.K),
             vx: float(c["var_x"]), vxe: float(c["sm"]["var_xe"]), "epsilon": float(c["sm"]["eps"]),
             "epsilon_m": float(c["sm"]["eps"]),                  # VERIFY: epsilon_m when not spherical
             "spherical": float(c["sm"]["spherical"])}
        if self.cov is None:
            v.update({"Var_e": float(c["var_e"]), "corr": float(c["corr"])})
        v.update(_solver_vals(c))
        v.update({"factor": self.which, "type": "test", "method": "repeated"})
        v["Cov"] = Matrix(c["cov"])
        if self.means is not None:
            v["means"] = Matrix(self.means)
        order = (["balanced", "N_per_group", "N_avg", "N_g", "N_rep", "N_a", "N", "power", "beta", "alpha"]
                 + _SOLVER + ["nfractional", "onesided", "Var_e", "corr", vx, "delta", "epsilon_m",
                              "epsilon", "spherical", vxe, "factor", "type", "method", "Cov", "means"])
        _store(self.s, v, order)


# ---------------------------------------------------------------------------------------------
# cmh
# ---------------------------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------------------------
# mcc
# ---------------------------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------------------------
# trend
# ---------------------------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------------------------
# cox
# ---------------------------------------------------------------------------------------------

class Cox(_M):
    method = "cox"
    iterates = False                       # closed forms: no iteration message (checked, case 1919)
    allowed = {"hratio", "sd", "r2", "eventprob", "failprob", "wdprob", "effect", "direction", "onesided"}
    _EFFN = {"coefficient": "coefficient", "hratio": "hazard ratio", "lnhratio": "log hazard-ratio"}

    def __init__(self, s, body, o):
        super().__init__(s, body, o)
        t = body.strip()
        self.b1 = P._values(t) if t else None
        if self.b1 is not None and o.get("hratio"):
            _fail("b1 and hratio() may not be combined")
        eff = str(o.get("effect") or "coefficient").lower()
        self.eff = next((e for e in self._EFFN if e.startswith(eff)), None)
        if self.eff is None:
            _fail("option effect() incorrectly specified")
        self.up_dir = str(o.get("direction") or "lower").lower().startswith("u")
        if self.vals("eventprob") and self.vals("failprob"):
            _fail("only one of eventprob() or failprob() may be specified")

    def base(self):
        b = self.common()
        if self.b1 is not None:
            b["b1"] = self.b1
        hr = self.vals("hratio")
        if hr is not None:
            b["hratio"] = hr
        for k, dflt in (("sd", [0.5]), ("r2", [0.0]), ("wdprob", [0.0])):
            b[k] = self.vals(k) or dflt
        b["pe"] = self.vals("eventprob") or self.vals("failprob") or [1.0]
        return b

    def title(self, what):
        w = {"n": "sample size", "power": "power", "effect": "target coefficient"}[what]
        return f"Estimated {w} for Cox PH regression"

    def hyp_lines(self, c, what):
        up = self.up_dir if what == "effect" else c["b1"] > 0
        extra = f"; beta1 {'>' if self.up_dir else '<'} 0" if (what == "effect" and self.sides() == 2) else ""
        return ["Wald test", f"Ho: beta1 = 0  versus  Ha: beta1 {self.hyp_op(c, up)} 0{extra}"]

    def compute(self, c):
        self.fix(c)
        a, sides = c["alpha"], self.sides()
        sd, r2, pe, wd = c["sd"], c["r2"], c["pe"], c["wdprob"]
        if sd <= 0 or not 0 <= r2 < 1 or not 0 < pe <= 1 or not 0 <= wd < 1:
            _fail("invalid argument value")                          # VERIFY
        b1 = c.get("b1")
        if b1 is None and c.get("hratio") is not None:
            b1 = math.log(c["hratio"])
        have_n = c.get("n") is not None
        if b1 is None:
            if have_n:
                what = "effect"
                c.setdefault("power", 0.8)
            else:
                b1 = math.log(0.5)                 # hratio() defaults to 0.5
                what = "n"
        if b1 is not None:
            if have_n:
                if c.get("power") is not None:
                    _fail("only two of n(), power() and the effect size may be specified")
                what = "power"
            else:
                what = "n"
        if what == "n":
            c.setdefault("power", 0.8)
            ev = C.cox_events(b1, sd, r2, a, sides, c["power"])
            n = ev / pe                          # N from the events before rounding (Stata 14: 164 with wdprob(.2))
            if wd > 0:
                n = n / (1 - wd)
            if not self.flag("nfractional"):
                n = math.ceil(n - 1e-9)
            ev = n * (1 - wd) * pe                  # events shown: from the final N (Stata 14)
            c["E"], c["N"] = float(ev if self.flag("nfractional") else math.ceil(ev - 1e-9)), float(n)
            c["power_out"] = C.cox_power(b1, sd, r2, n * (1 - wd), pe, a, sides)
        elif what == "power":
            c["N"] = float(c["n"])
            c["power_out"] = C.cox_power(b1, sd, r2, c["N"], pe, a, sides)
            c["E"] = c["N"] * pe if self.flag("nfractional") else math.ceil(c["N"] * pe - 1e-9)
        else:
            c["N"] = float(c["n"])
            mag = C.cox_effect(sd, r2, c["N"], pe, a, sides, c["power"])
            b1 = mag if self.up_dir else -mag
            c["power_out"] = c["power"]
            c["E"] = c["N"] * pe if self.flag("nfractional") else math.ceil(c["N"] * pe - 1e-9)
        c["b1"], c["hr"] = b1, math.exp(b1)
        c["delta"] = {"coefficient": b1, "hratio": math.exp(b1), "lnhratio": b1}[self.eff]
        return what

    def show(self, c, what):
        items = [("alpha", P._f4(c["alpha"]))]
        if what != "power":
            items.append(("power", P._f4(c["power"])))
        if what != "n":
            items.append(("N", _I(c["N"])))
        if what != "effect":
            items.append(("delta", P._f4(c["delta"]) + f"  ({self._EFFN[self.eff]})"))
            items.append(("hratio" if self.eff == "coefficient" else "b1",
                          P._f4(c["hr"] if self.eff == "coefficient" else c["b1"])))
        items.append(("sd", P._f4(c["sd"])))
        if self.o.get("r2"):
            items.append(("R2", P._f4(c["r2"])))
        self.put(items)
        cens = [("Pr_E", P._f4(c["pe"]))] + ([("Pr_w", P._f4(c["wdprob"]))] if self.o.get("wdprob") else [])
        if what == "n":
            self.results("Censoring and withdrawal:" if self.o.get("wdprob") else "Censoring:", cens)
            self.results("Estimated number of events and sample size:", [("E", _I(c["E"])), ("N", _I(c["N"]))])
        elif what == "power":
            self.results("Number of events and censoring:", [("E", _I(c["E"]))] + cens)
            self.results("Estimated power:", [("power", P._f4(c["power_out"]))])
        else:
            self.results("Number of events and censoring:", [("E", _I(c["E"]))] + cens)
            self.results(f"Estimated effect size and target {self._EFFN[self.eff]}:",
                         [("delta", P._f4(c["delta"]) + f"  ({self._EFFN[self.eff]})"), ("b1", P._f4(c["b1"]))])
        self.returns(c, what)

    def columns(self, c, what):
        return ["alpha", "power", "N", "E", "delta", "hratio" if self.eff == "coefficient" else "b1", "sd"] + (
            ["R2"] if self.o.get("r2") else []) + ["Pr_E"]

    def row(self, c, what):
        return {"alpha": c["alpha"], "power": c["power_out"] if what == "power" else c["power"],
                "N": c["N"], "E": c["E"], "delta": c["delta"], "b1": c["b1"], "hratio": c["hr"],
                "sd": c["sd"], "R2": c["r2"], "Pr_E": c["pe"]}

    _ORDER = ["N", "power", "beta", "alpha", "delta", "nfractional", "onesided", "E", "hratio", "b1", "sd",
              "R2", "Pr_E", "Pr_w", "type", "method", "effect"]

    def returns(self, c, what):
        pw = c["power_out"] if what == "power" else c["power"]
        v = {"N": c["N"], "power": pw, "beta": 1 - pw, "alpha": c["alpha"], "delta": c["delta"],
             "nfractional": float(self.flag("nfractional")), "onesided": float(self.sides() == 1),
             "E": c["E"], "hratio": c["hr"], "b1": c["b1"], "sd": c["sd"], "R2": c["r2"], "Pr_E": c["pe"],
             "Pr_w": c["wdprob"], "type": "test", "method": "cox", "effect": self.eff}
        _store(self.s, {k: (float(x) if not isinstance(x, str) else x) for k, x in v.items()}, self._ORDER)


# ---------------------------------------------------------------------------------------------
# exponential
# ---------------------------------------------------------------------------------------------

class Exponential(_M):
    method = "exponential"
    iterates = False
    allowed = {"time", "nratio", "hratio", "lnhratio", "hdifference", "difference", "diff",
               "loghazard", "unconditional", "onesided", "studytime", "fperiod", "aperiod",
               "aprob", "aptime", "atime", "ashape", "lossprob", "lossprob1", "lossprob2",
               "losstime", "losshaz", "losshaz1", "losshaz2", "effect"}

    def __init__(self, s, body, o):
        super().__init__(s, body, o)
        if "effect" in o:
            _fail("option effect() is not available yet in OpenDTA")
        toks = P._split_args(body)
        if len(toks) > 2:
            _fail("too many values specified")
        self.args = [P._values(t) for t in toks]
        self.surv = bool(o.get("time"))
        self.logtest = self.flag("loghazard") or not self.args
        self.unc = self.flag("unconditional")
        # accrual and follow-up (scalars)
        T, f, r = self.one("studytime"), self.one("fperiod"), self.one("aperiod")
        if r is not None and f is None and T is None:
            _fail("two of studytime(), fperiod() and aperiod() are required")    # VERIFY
        self.r = r or 0.0
        if T is not None and f is not None:
            self.r = T - f if r is None else r
            if abs(self.r + f - T) > 1e-9:
                _fail("studytime(), fperiod() and aperiod() are inconsistent")
        elif T is not None:
            f = T - self.r
        elif f is not None:
            T = self.r + f
        self.T, self.f = T, f
        if self.r < 0 or (T is not None and T <= 0):
            _fail("invalid accrual or follow-up period")
        self.gamma = self.one("ashape", 0.0)
        self.pa = self.ta = None
        if any(o.get(k) for k in ("aprob", "aptime", "atime")):
            if self.r <= 0 or o.get("ashape"):
                _fail("aprob(), aptime() and atime() require aperiod() and may not be combined with ashape()")
            self.pa = self.one("aprob", 0.5)
            self.ta = self.one("atime") if o.get("atime") else self.one("aptime", 0.5) * self.r
            if not 0 < self.ta <= self.r or not 0 < self.pa < 1:
                _fail("option atime(), aptime() or aprob() incorrectly specified")
            self.gamma = C.exp_gamma(self.r, self.ta, self.pa)
        # losses
        self.eta1 = self.eta2 = 0.0
        self.loss = None
        tl = self.one("losstime", 1.0)
        lp, lp1, lp2 = self.one("lossprob"), self.one("lossprob1"), self.one("lossprob2")
        lh, lh1, lh2 = self.one("losshaz"), self.one("losshaz1"), self.one("losshaz2")
        if (lp is not None or lp1 is not None or lp2 is not None) and (lh is not None or lh1 is not None
                                                                         or lh2 is not None):
            _fail("lossprob() and losshaz() may not be combined")
        if any(v is not None for v in (lp, lp1, lp2, lh, lh1, lh2)):
            if T is None:
                _fail("loss options require aperiod() or fperiod()")
            def eta_of(lpg, lpc, lhg, lhc):
                if lpg is not None or lpc is not None:
                    return -math.log(1 - (lpg if lpg is not None else lpc)) / tl
                return lhg if lhg is not None else (lhc or 0.0)
            self.eta1 = eta_of(lp1, lp, lh1, lh)
            self.eta2 = eta_of(lp2, lp, lh2, lh)
            if any(v is not None for v in (lp, lp1, lp2)):
                self.loss = [("lpr1", 1 - math.exp(-self.eta1 * tl)), ("lpr2", 1 - math.exp(-self.eta2 * tl)),
                             ("lt", tl)]
            else:
                self.loss = [("lh1", self.eta1), ("lh2", self.eta2)]                # VERIFY: names with losshaz()
        self.up_dir = None

    def base(self):
        b = self.common()
        names = ("s1", "s2") if self.surv else ("h1", "h2")
        for nm, v in zip(names, self.args):
            b[nm] = v
        for k in ("hratio", "lnhratio", "hdifference", "time"):
            v = self.vals(k)
            if v is not None:
                b[k] = v
        d = self.vals("diff") or self.vals("difference")
        if d is not None:
            b["hdifference"] = d
        b["nratio"] = self.vals("nratio") or [1.0]
        return b

    def title(self, what):
        return ("Estimated sample sizes for two-sample comparison of survivor functions" if what == "n"
                else "Estimated power for two-sample comparison of survivor functions")

    def hyp_lines(self, c, what):
        lt = c["logtest"]
        up = c["hr"] > 1
        if lt:
            return [c["test"], f"Ho: ln(HR) = 0  versus  Ha: ln(HR) {self.hyp_op(c, up)} 0"]
        return [c["test"], f"Ho: h2 = h1  versus  Ha: h2 {self.hyp_op(c, up)} h1"]

    def compute(self, c):
        self.fix(c)
        a, sides, R = c["alpha"], self.sides(), c["nratio"]
        if R <= 0:
            _fail("nratio() must be positive")
        l1 = l2 = None
        if self.surv:
            t = c.get("time")
            if t is None or t <= 0:
                _fail("time() is required with survival probabilities")
            if c.get("s1") is not None:
                l1 = -math.log(c["s1"]) / t
            if c.get("s2") is not None:
                l2 = -math.log(c["s2"]) / t
        else:
            l1, l2 = c.get("h1"), c.get("h2")
        eff = [k for k in ("hratio", "lnhratio", "hdifference") if c.get(k) is not None]
        if l1 is not None and l2 is not None:
            if eff:
                _fail("hratio(), lnhratio() and hdifference() may not be combined with both groups given")
        elif len(eff) > 1:
            _fail("only one of hratio(), lnhratio() and hdifference() may be specified")
        elif l1 is not None:
            if eff == ["hdifference"]:
                l2 = l1 + c["hdifference"]
            elif eff == ["lnhratio"]:
                l2 = l1 * math.exp(c["lnhratio"])
            else:
                l2 = l1 * (c["hratio"] if eff else 0.5)
        else:
            if eff == ["hdifference"]:
                _fail("hdifference() requires the control-group hazard or survival")
            hr0 = c["hratio"] if eff == ["hratio"] else (math.exp(c["lnhratio"]) if eff else 0.5)
            if self.T is not None:
                _fail("the control-group hazard or survival is required with accrual, follow-up or losses")
            l1, l2 = 1.0, hr0                   # VERIFY: without censoring only the ratio matters
            c["hz_unknown"] = True
        if min(l1, l2) <= 0:
            _fail("hazards must be positive")
        hr = l2 / l1
        logtest = self.logtest
        c.update({"l1": l1, "l2": l2, "hr": hr, "logtest": logtest})
        delta = math.log(hr) if logtest else l2 - l1
        c["delta"] = delta
        c["test"] = ("Exponential test, log hazard-ratio, " if logtest else "Exponential test, hazard difference, ") + \
            ("unconditional" if self.unc else "conditional")
        xo, xa = C.exp_xi(l1, l2, self.eta1, self.eta2, R, self.T, self.r, self.gamma, logtest, self.unc)
        have_n = c.get("n") is not None
        if have_n:
            if c.get("power") is not None:
                _fail("only two of n(), power() and the effect size may be specified")
            c["power_out"] = C.exp_power(delta, xo, xa, float(c["n"]), a, sides, self.unc)
            N = float(c["n"])
            n1 = N / (1 + R)
            n2 = N - n1
            what = "power"
        else:
            c.setdefault("power", 0.8)
            Nc = C.exp_size(delta, xo, xa, a, sides, c["power"], self.unc)
            if self.flag("nfractional"):
                n1 = Nc / (1 + R)
                n2 = R * n1
            else:
                n1 = float(math.ceil(Nc / (1 + R) - 1e-9))        # VERIFY: rounding of N1 and N2
                n2 = float(math.ceil(R * n1 - 1e-9))
            N = n1 + n2
            c["power_out"] = C.exp_power(delta, xo, xa, N, a, sides, self.unc)
            what = "n"
        c["N"], c["N1"], c["N2"] = N, n1, n2
        c["R_a"] = n2 / n1
        return what

    def show(self, c, what):
        out = self.s.output
        items = [("alpha", P._f4(c["alpha"]))]
        if what == "n":
            items.append(("power", P._f4(c["power"])))
        else:
            items += [("N", _I(c["N"])), ("N per group", _I(c["N1"])) if c["N1"] == c["N2"] else ("N1", _I(c["N1"]))]
            if c["N1"] != c["N2"]:
                items.append(("N2", _I(c["N2"])))
        lab = "log hazard-ratio" if c["logtest"] else "hazard difference"
        items.append(("delta", P._f4(c["delta"]) + f"  ({lab})"))
        if c["nratio"] != 1:
            items.append(("N2/N1", P._f4(c["nratio"])))
        self.put(items)
        if self.T is not None:
            out.write("\nAccrual and follow-up information:\n\n", "text")
            self.line("duration", P._f4(self.T))
            self.line("follow-up", P._f4(self.f))
            self.line("accrual", P._f4(self.r) + ("  (uniform)" if abs(self.gamma) < 1e-6 else "  (exponential)"))
            if self.pa is not None:
                self.line("accrual(%)", P._f4(100 * self.pa))
                self.line("by time", P._f4(self.ta))
        out.write("\nSurvival information:\n\n", "text")
        if c.get("hz_unknown"):
            self.line("hratio", P._f4(c["hr"]))
        elif self.surv:
            for h, nm in ((c["l1"], "1"), (c["l2"], "2")):
                sv = c.get("s" + nm) if c.get("s" + nm) is not None else math.exp(-h * c["time"])
                out.write(f"{'h' + nm:>{self.w}} = ", "text")
                out.write(P._f4(h), "result")
                out.write(f"{'s' + nm:>{self.w}} = ", "text")
                out.write(P._f4(sv) + "\n", "result")
            out.write(" " * (self.w + 12) + f"{'t':>{self.w}} = ", "text")
            out.write(P._f4(c["time"]) + "\n", "result")
        else:
            self.line("h1", P._f4(c["l1"]))
            self.line("h2", P._f4(c["l2"]))
        if self.loss is not None:
            out.write("\nLoss-to-follow-up information:\n\n", "text")
            for nm, v in self.loss:
                self.line(nm, P._f4(v))
        if what == "n":
            res = [("N", _I(c["N"]))]
            if c["N1"] == c["N2"]:
                res.append(("N per group", _I(c["N1"])))
            else:
                res += [("N1", _I(c["N1"])), ("N2", _I(c["N2"]))]
            self.results("Estimated sample sizes:", res)
        else:
            self.results("Estimated power:", [("power", P._f4(c["power_out"]))])
        if c.get("hz_unknown"):
            out.write("Note: Individual hazard rates are not reported because a\n"
                      "      control-group hazard rate is not specified.\n", "text")
        self.returns(c, what)

    def header(self, c, what, lines):
        if self.args:
            self.s.output.write("note: input parameters are " + ("survival probabilities" if self.surv else "hazard rates")
                                + "\n", "text")
        super().header(c, what, lines)

    def columns(self, c, what):
        return ["alpha", "power", "N", "delta", "h1", "h2", "hratio"]

    def row(self, c, what):
        return {"alpha": c["alpha"], "power": c["power_out"] if what == "power" else c["power"], "N": c["N"],
                "delta": c["delta"], "h1": c["l1"], "h2": c["l2"], "hratio": c["hr"]}

    _ORDER = ["power", "beta", "alpha", "delta", "nfractional", "onesided", "N", "N1", "N2", "nratio",
              "hratio", "lnhratio", "h1", "h2", "hdiff", "unconditional", "test", "effect", "type", "method"]

    def returns(self, c, what):
        # r(E0_vec), r(Ea_vec) and r(Pr_vec) of Stata are not stored: the manual does not define them
        pw = c["power_out"] if what == "power" else c["power"]
        v = {"power": pw, "beta": 1 - pw, "alpha": c["alpha"], "delta": c["delta"],
             "nfractional": float(self.flag("nfractional")), "onesided": float(self.sides() == 1),
             "N": c["N"], "N1": c["N1"], "N2": c["N2"], "nratio": c["nratio"], "hratio": c["hr"],
             "lnhratio": math.log(c["hr"]), "h1": c["l1"], "h2": c["l2"], "hdiff": c["l2"] - c["l1"],
             "unconditional": float(self.unc),
             "test": "log hazard-ratio" if c["logtest"] else "hazard difference",
             "effect": "lnhratio" if c["logtest"] else "difference",          # VERIFY: effect with loghazard
             "type": "test", "method": "exponential"}
        if self.surv:
            v.update({"s1": math.exp(-c["l1"] * c["time"]), "s2": math.exp(-c["l2"] * c["time"]),
                      "time": c["time"]})
        if self.T is not None:
            v.update({"aperiod": self.r, "fperiod": self.f, "studytime": self.T, "ashape": self.gamma})
        if self.eta1 or self.eta2:
            v.update({"losshaz1": self.eta1, "losshaz2": self.eta2})
        _store(self.s, {k: (float(x) if not isinstance(x, str) else x) for k, x in v.items()}, self._ORDER)


# ---------------------------------------------------------------------------------------------
# logrank
# ---------------------------------------------------------------------------------------------

class Logrank(_M):
    method = "logrank"
    iterates = False
    allowed = {"nratio", "hratio", "lnhratio", "schoenfeld", "direction", "onesided", "simpson", "wdprob"}

    def __init__(self, s, body, o):
        super().__init__(s, body, o)
        toks = P._split_args(body)
        if len(toks) > 2:
            _fail("too many values specified")
        self.args = [P._values(t) for t in toks]
        self.schoen = self.flag("schoenfeld")
        self.up_dir = str(o.get("direction") or "lower").lower().startswith("u")
        self.simp = None
        if o.get("simpson") not in (None, True):
            t = str(o["simpson"]).strip()
            v = (list(np.array(get_matrix(s, t).data, dtype=float).reshape(-1))
                 if re.fullmatch(r"[A-Za-z_]\w*", t) else P._values(t))
            if len(v) != 3 or any(not 0 < x < 1 for x in v):
                _fail("option simpson() incorrectly specified")
            self.simp = tuple(v)
            if self.args:
                _fail("simpson() may not be combined with surv1 and surv2")

    def base(self):
        b = self.common()
        for nm, v in zip(("s1", "s2"), self.args):
            b[nm] = v
        for k in ("hratio", "lnhratio", "wdprob"):
            v = self.vals(k)
            if v is not None:
                b[k] = v
        b["nratio"] = self.vals("nratio") or [1.0]
        return b

    def title(self, what):
        w = {"n": "sample sizes", "power": "power", "effect": "hazard ratio"}[what]
        return f"Estimated {w} for two-sample comparison of survivor functions"

    def hyp_lines(self, c, what):
        up = self.up_dir if what == "effect" else c["hr"] > 1
        op = self.hyp_op(c, up)
        sub = "Log-rank test, Schoenfeld method" if self.schoen else "Log-rank test, Freedman method"
        hyp = f"Ho: ln(HR) = 0  versus  Ha: ln(HR) {op} 0" if self.schoen else f"Ho: HR = 1  versus  Ha: HR {op} 1"
        if what == "effect" and self.sides() == 2:
            hyp += f"; {'ln(HR)' if self.schoen else 'HR'} {'>' if self.up_dir else '<'} {'0' if self.schoen else '1'}"
        return [sub, hyp]

    def compute(self, c):
        self.fix(c)
        a, sides, R = c["alpha"], self.sides(), c["nratio"]
        s1, s2 = c.get("s1"), c.get("s2")
        wd = c.get("wdprob", 0.0)
        eff = [k for k in ("hratio", "lnhratio") if c.get(k) is not None]
        if len(eff) > 1 or (eff and s1 is not None and s2 is not None):
            _fail("hratio() and lnhratio() may not be combined with each other or with surv2")
        have_n = c.get("n") is not None
        if have_n and c.get("power") is not None:
            what = "effect"
            if eff or s2 is not None:
                _fail("the hazard ratio may not be specified when computing the effect size")
        elif have_n:
            what = "power"
        else:
            what = "n"
            c.setdefault("power", 0.8)
        sch = self.schoen
        simp = self.simp

        def pe_of(hr):
            if simp is not None:
                return C.logrank_pe(hr, R, None, None, simp)
            return C.logrank_pe(hr, R, s1, None)

        if what == "effect":
            N = float(c["n"])
            sgn = 1.0 if self.up_dir else -1.0
            try:
                t = C.solve_up(lambda t: C.logrank_power(math.exp(sgn * t), R, N, pe_of(math.exp(sgn * t)),
                                                         a, sides, sch) - c["power"], 1e-9, 1.0)
            except ValueError:
                raise StataError(498, "failed to find a solution")
            hr = math.exp(sgn * t)
            c["_init"] = 1.0
        else:
            if s1 is not None and s2 is not None:
                hr = math.log(s2) / math.log(s1)
            elif eff:
                hr = c["hratio"] if eff == ["hratio"] else math.exp(c["lnhratio"])
            else:
                hr = 0.5
        if hr <= 0 or hr == 1:
            _fail("the hazard ratio must be positive and different from 1")
        pe = pe_of(hr)
        c.update({"hr": hr, "pe": pe, "wd": wd})
        c["delta"] = math.log(hr) if sch else hr
        if what == "n":
            ev = C.logrank_events(hr, R, a, sides, c["power"], sch)
            n_c = ev / pe / (1 - wd)                # N from the events before rounding (Stata 14, wdprob(.1))
            if self.flag("nfractional"):
                n1 = n_c / (1 + R)
                n2 = R * n1
            else:
                n1 = float(math.ceil(n_c / (1 + R) - 1e-9))
                n2 = float(math.ceil(R * n1 - 1e-9))
            N = n1 + n2
            ev = N * (1 - wd) * pe                  # events shown: from the final N (Stata 14)
            c["E"] = float(ev if self.flag("nfractional") else math.ceil(ev - 1e-9))
            c["power_out"] = C.logrank_power(hr, R, N * (1 - wd), pe, a, sides, sch)
        else:
            N = float(c["n"])
            n1 = N / (1 + R)
            n2 = N - n1
            ev = N * (1 - wd) * pe
            c["E"] = ev if self.flag("nfractional") else math.ceil(ev - 1e-9)             # Stata 14 shows 81 for 80.46
            c["power_out"] = c["power"] if what == "effect" else C.logrank_power(hr, R, N, pe, a, sides, sch)
        c["N"], c["N1"], c["N2"] = N, n1, n2
        return what

    def header(self, c, what, lines):
        out = self.s.output
        if self.simp is not None:
            s1 = self.simp
            s2 = [x ** c["hr"] for x in s1]
            out.write("note: probability of an event is computed using Simpson's rule with\n", "text")
            for g, v in (("1", s1), ("2", s2)):
                out.write(f"      S{g}(f) = {v[0]:.2f}, S{g}(f+r/2) = {v[1]:.2f}, S{g}(T) = {v[2]:.2f}\n", "text")
        self.iterates = what == "effect"              # the other cases are closed forms
        super().header(c, what, lines)

    def show(self, c, what):
        items = [("alpha", P._f4(c["alpha"]))]
        if what != "power":
            items.append(("power", P._f4(c["power"])))
        if what != "n":
            items += [("N", _I(c["N"])), ("N per group", _I(c["N1"]))] if c["N1"] == c["N2"] else [
                ("N", _I(c["N"])), ("N1", _I(c["N1"])), ("N2", _I(c["N2"]))]
        lab = "log hazard-ratio" if self.schoen else "hazard ratio"
        if what != "effect":
            items.append(("delta", P._f4(c["delta"]) + f"  ({lab})"))
            items.append(("hratio", P._f4(c["hr"])))
        if c["nratio"] != 1:
            items.append(("N2/N1", P._f4(c["nratio"])))
        self.put(items)
        cens = []
        if c.get("s1") is not None:
            cens += [("s1", P._f4(c["s1"])), ("s2", P._f4(c["s1"] ** c["hr"]))]
        cens.append(("Pr_E", P._f4(c["pe"])))
        if c["wd"]:
            cens.append(("Pr_w", P._f4(c["wd"])))
        if what == "n":
            self.results("Censoring and withdrawal:" if c["wd"] else "Censoring:", cens)
            res = [("E", _I(c["E"])), ("N", _I(c["N"]))]
            res += ([("N per group", _I(c["N1"]))] if c["N1"] == c["N2"] else
                    [("N1", _I(c["N1"])), ("N2", _I(c["N2"]))])
            self.results("Estimated number of events and sample sizes:", res)
        elif what == "power":
            self.results("Number of events and censoring:", [("E", _I(c["E"]))] + cens)
            self.results("Estimated power:", [("power", P._f4(c["power_out"]))])
        else:
            self.results("Number of events and censoring:", [("E", _I(c["E"]))] + cens)
            self.results("Estimated effect size and hazard ratio:",
                         [("delta", P._f4(c["delta"]) + f"  ({lab})"), ("hratio", P._f4(c["hr"]))])
        self.returns(c, what)

    def columns(self, c, what):
        return ["alpha", "power", "N", "N1", "N2", "E", "delta", "hratio", "Pr_E"]

    def row(self, c, what):
        return {"alpha": c["alpha"], "power": c["power_out"] if what == "power" else c["power"],
                "N": c["N"], "N1": c["N1"], "N2": c["N2"], "E": c["E"], "delta": c["delta"],
                "hratio": c["hr"], "Pr_E": c["pe"]}

    _ORDER = (["power", "beta", "alpha"] + _SOLVER + ["delta", "nfractional", "onesided", "N", "N1", "N2",
              "nratio", "hratio", "E", "s1", "s2", "Pr_E", "Pr_w", "type", "method", "test", "effect"])

    def returns(self, c, what):
        pw = c["power_out"] if what == "power" else c["power"]
        v = {"power": pw, "beta": 1 - pw, "alpha": c["alpha"], "delta": c["delta"],
             "nfractional": float(self.flag("nfractional")), "onesided": float(self.sides() == 1),
             "N": c["N"], "N1": c["N1"], "N2": c["N2"], "nratio": c["nratio"], "hratio": c["hr"],
             "E": c["E"], "Pr_E": c["pe"], "Pr_w": c["wd"]}
        if c.get("s1") is not None:
            v.update({"s1": c["s1"], "s2": c["s1"] ** c["hr"]})
        v.update(_solver_vals(c))
        v.update({"type": "test", "method": "logrank", "test": "Schoenfeld" if self.schoen else "Freedman",
                  "effect": "hratio"})                                  # VERIFY: effect with schoenfeld
        _store(self.s, {k: (float(x) if not isinstance(x, str) else x) for k, x in v.items()}, self._ORDER)


_CLASSES = {"twoway": TwoWay, "repeated": Repeated, "cmh": CMH, "mcc": MCC, "trend": Trend,
            "cox": Cox, "exponential": Exponential, "logrank": Logrank}


def cmd(s: "Session", method: str, body: str, o: dict) -> None:
    M = _CLASSES[method](s, body, o)
    configs, multi = _combos(M.base(), bool(o.get("parallel")))
    rows, what = [], None
    for c in configs:
        c = {k: v for k, v in c.items()}
        w = M.compute(c)
        if what is None:
            what = w
        elif w != what:
            _fail("the values given ask for different computations")        # VERIFY
        rows.append(c)
    if (multi or o.get("table")) and not o.get("notable"):
        M.table(rows, what)
        M.returns(rows[-1], what)
    else:
        M.single(rows[-1], what)


def _install() -> None:
    P._METHODS.update(METHODS)
    P._OPTS.update({k: v for k, v in _NEW_OPTS.items() if k not in P._OPTS})
    P._RAW.update(METHODS)


_install()

"""Shared parts of the power methods of part 3 (commands/power3.py): the column formats of the
tables, the common options, the solver's r() items, the order of r() and the base class _M, which
reads the options and writes the output of every method."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

from ..core import power_suites as C
from ..core.errors import StataError
from . import power as P

if TYPE_CHECKING:
    from ..session import Session


_COMMON = {"alpha", "power", "beta", "n", "nfractional", "parallel", "table", "notable", "nolog",
           "log", "dots", "nodots", "init", "iterate", "tolerance", "ftolerance", "notitle"}


_NOT_YET = {"graph", "saving", "show", "st1", "showmatrices", "showmeans", "showasmatrix",
            "n1", "n2", "effect_exp", "gridlines"}


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

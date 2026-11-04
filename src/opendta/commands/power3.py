"""power ([PSS] power), part 3: twoway, repeated, cmh, mcc, trend, cox, exponential and logrank.

The formulas are in core/power_suites.py (the Methods and formulas of each [PSS] page); this
module dispatches each method to its class (power3_anova, power3_tables, power3_surv, on the base
class of power3_base), which reads the arguments and options, chooses what to compute (sample
size, power or effect size), writes the output and stores r(). The commands are added to the method tables of
commands/power.py when this module is imported.

The layout of the output, the order of r() and the rounding of the sizes are only partly given
by the manual and are marked VERIFY (case compat/do/1919_other_suites2.do). Options that are
not implemented raise "not available yet": graph(), saving(), show, st1(), n#()."""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

from . import power as P
from .power3_anova import Repeated, TwoWay
from .power3_base import _fail
from .power3_surv import Cox, Exponential, Logrank
from .power3_tables import CMH, MCC, Trend

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


# ---------------------------------------------------------------------------------------------
# logrank
# ---------------------------------------------------------------------------------------------


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

"""signestimationsample, checkestimationsample, makecns and matcproc ([P] signestimationsample,
[P] makecns): tools for programs that write estimation commands."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core.constraints import constraint_transform
from ..core.dataset import NAME_RE, Dataset, Variable
from ..core.errors import StataError
from ..lang.syntax import match_options
from .datacmds import _signature
from .matrix import Matrix, _store, get_matrix
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def _sample_signature(s: "Session", names: list[str]) -> str:
    """The signature of `names` over the observations in e(sample)."""
    ds = s.data
    mask = ds.esample
    if mask is None:
        raise StataError(301, "last estimates not found")
    sub = Dataset()
    sub.nobs = int(mask.sum())
    for n in names:
        v = ds.get(n)
        nv = Variable(v.name, v.vtype, np.empty(0))
        nv.raw = v.raw[mask]
        sub.vars.append(nv)
    # e(datasignature) has no variable-type part: 30:3:chk1:chk2 (compat 1801)
    return re.sub(r"\(\d+\)", "", _signature(sub), count=1)


@command("signestimationsample")
def cmd_signestimationsample(s: "Session", args: str) -> None:
    if not args.strip():
        raise StataError(100, "varlist required")
    names = s.expand_varlist(args.strip())
    wexp = s.e.get("wexp")
    if isinstance(wexp, str):
        # the variables of the weight expression enter the signature as well
        for w in re.findall(r"[A-Za-z_]\w*", wexp):
            if s.data.has(w) and w not in names:
                names.append(w)
    sig = _sample_signature(s, names)
    s.e["datasignaturevars"] = " ".join(names)
    s.e["datasignature"] = sig


@command("checkestimationsample")
def cmd_checkestimationsample(s: "Session", args: str) -> None:
    if args.strip():
        raise StataError(198, "invalid syntax")
    names = str(s.e.get("datasignaturevars", "")).split()
    if not names or "datasignature" not in s.e:
        raise StataError(301, "last estimates not found")                # VERIFY
    if any(not s.data.has(n) for n in names) or _sample_signature(s, names) != s.e["datasignature"]:
        raise StataError(459, "data have changed since estimation")


def _constraint_text(row: np.ndarray, names: list[str]) -> str:
    """'x1 - x2 = 0': the coefficients of the row (the last element is the right-hand side)."""
    from ..core.formats import format_value
    text = ""
    for coef, name in zip(row[:-1], names):
        if coef == 0:
            continue
        mag = abs(float(coef))
        term = name if mag == 1 else f"{format_value(mag, '%9.0g', pad=False).strip()}*{name}"
        if not text:
            text = term if coef > 0 else f"-{term}"
        else:
            text += f" + {term}" if coef > 0 else f" - {term}"
    return f"{text or '0'} = {format_value(float(row[-1]), '%9.0g', pad=False).strip()}"   # VERIFY


@command("makecns")
def cmd_makecns(s: "Session", args: str) -> None:
    head, comma, opts = args.partition(",")
    o = match_options(opts, {"nocnsnotes": 5, "displaycns": 4, "r": 1}) if comma else {}
    name = head.strip()
    if not name:
        return                                    # VERIFY: Stata does nothing without arguments
    if not NAME_RE.match(name):
        raise StataError(198, "option makecns with a list of constraint numbers is not available yet "
                              "in OpenDTA")
    if o.get("displaycns") and not o.get("r"):
        raise StataError(198, "option displaycns is not available yet in OpenDTA")
    b = s.e.get("b")
    if not isinstance(b, Matrix):
        raise StataError(301, "last estimates not found")
    try:
        C = get_matrix(s, name)
        if C.cols != b.cols + 1:
            raise StataError(503, "conformability error")
    except StataError as e:
        # a bad constraint matrix is only reported (VERIFY: the return code afterwards)
        s.output.write(f"(note: constraint matrix '{name}' caused error r({e.rc}))\n", "text")
        return
    s.e["Cns"] = Matrix(C.data.copy(), list(C.rownames), list(b.colnames) + ["_r"])
    if o.get("r"):
        # return list shows the results in the reverse order of this dictionary
        s.r = {"k_autoCns": 0.0, "k": float(C.rows)}
        for k in range(C.rows):
            s.r[f"cns{k + 1}"] = _constraint_text(C.data[k], list(b.colnames))
        s.r["clist"] = name


@command("matcproc")
def cmd_matcproc(s: "Session", args: str) -> None:
    words = args.split()
    if len(words) < 3:
        raise StataError(100, "namelist:  too few specified")     # VERIFY: return code
    if len(words) != 3 or not all(NAME_RE.match(w) for w in words):
        raise StataError(198, "invalid syntax")
    cns = s.e.get("Cns")
    if not isinstance(cns, Matrix):
        raise StataError(301, "last estimates not found")                # VERIFY
    b = s.e["b"]
    T, a = constraint_transform(cns.data)
    store = _store(s)
    store[words[0]] = Matrix(T, list(b.colnames), [f"c{j}" for j in range(1, T.shape[1] + 1)])
    store[words[1]] = Matrix(a, ["r1"], list(b.colnames))
    store[words[2]] = cns.copy()

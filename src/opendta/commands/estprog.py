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


def _b_est(b: Matrix):
    """The coefficient names of e(b), in the form the constraint parser of postest reads."""
    from types import SimpleNamespace
    eqs = ["" if q == "_" else q for q in (b.coleq or ["_"] * b.cols)]
    return SimpleNamespace(names=list(b.colnames), eqnames=eqs if any(eqs) else [])


def constraint_lines(b: Matrix, C: np.ndarray) -> list[str]:
    """'x1 - x2 = 0' for each row of C = (R, r), with [eq] prefixes in multiple-equation models."""
    from ..estimation.postest import _ctext
    est = _b_est(b)
    return [_ctext(est, {j: float(v) for j, v in enumerate(row[:-1]) if v != 0}, float(row[-1]))
            for row in C]


def _from_clist(s: "Session", b: Matrix, clist: str, notes: bool) -> tuple[np.ndarray, list[int]]:
    """The rows of the defined constraints in clist that can be applied, in order. A constraint that
    cannot be read, or that is redundant with or contradicts the earlier ones, is left out with a
    note. VERIFY: the wording of the notes and r(412) for a redundant constraint."""
    from ..estimation.postest import _constraint
    from ..estimation.sysmodels import _constraints, _numlist
    try:
        numbers = _numlist(clist)
    except ValueError:
        raise StataError(198, "invalid numlist")
    defined = _constraints(s)
    est = _b_est(b)
    rows: list[np.ndarray] = []
    used: list[int] = []
    for n in numbers:
        try:
            if n not in defined:
                raise StataError(111, f"constraint {n} not found")
            coef, rhs = _constraint(s, est, defined[n])
            row = np.zeros(b.cols + 1)
            for j, v in coef.items():
                row[j] += v
            row[-1] = rhs
            if not row[:-1].any() or (rows and np.linalg.matrix_rank(np.array(rows + [row])[:, :-1])
                                      == np.linalg.matrix_rank(np.array(rows)[:, :-1])):
                raise StataError(412, "redundant or inconsistent constraints")
        except StataError as e:
            if notes:
                s.output.write(f"(note: constraint number {n} caused error r({e.rc}))\n", "text")
            continue
        rows.append(row)
        used.append(n)
    return np.array(rows).reshape(-1, b.cols + 1), used


@command("makecns")
def cmd_makecns(s: "Session", args: str) -> None:
    head, comma, opts = args.partition(",")
    o = match_options(opts, {"nocnsnotes": 5, "displaycns": 4, "r": 1, "suminfo": 3}) if comma else {}
    name = head.strip()
    if not name:
        return                                    # nothing to do (compat 1801)
    b = s.e.get("b")
    if not isinstance(b, Matrix):
        raise StataError(301, "last estimates not found")
    if NAME_RE.match(name):
        try:
            C = get_matrix(s, name)
            if C.cols != b.cols + 1:
                raise StataError(503, "conformability error")
        except StataError as e:
            # a bad constraint matrix is only reported (VERIFY: the return code afterwards)
            s.output.write(f"(note: constraint matrix '{name}' caused error r({e.rc}))\n", "text")
            return
        rows, clist, rownames = C.data.copy(), name, list(C.rownames)
    else:
        rows, used = _from_clist(s, b, name, not o.get("nocnsnotes"))
        clist, rownames = " ".join(map(str, used)), [f"r{k}" for k in range(1, len(used) + 1)]
    if len(rows):
        s.e["Cns"] = Matrix(rows, rownames, list(b.colnames) + ["_r"], ["_"] * len(rows),
                            list(b.coleq or ["_"] * b.cols) + ["_"])
    else:
        s.e.pop("Cns", None)
    texts = constraint_lines(b, rows)
    if o.get("displaycns"):
        _display_cns(s, texts)
    if o.get("r"):
        # return list shows the results in the reverse order of this dictionary
        s.r = {"k_autoCns": 0.0, "k": float(len(rows))}
        for k, t in enumerate(texts, 1):
            s.r[f"cns{k}"] = t
        s.r["clist"] = clist


def _display_cns(s: "Session", texts: list[str]) -> None:
    """' ( 1)  x1 - x2 = 0', as in the header of estimation commands with constraints."""
    for i, t in enumerate(texts, 1):
        s.output.write(f" ({i:>2})  {t}\n", "text")


def dispcns(s: "Session", rest: str) -> None:
    """matrix dispCns [, r]: the constraints of e(Cns)."""
    head, comma, opts = rest.partition(",")
    if head.strip():
        raise StataError(198, "invalid syntax")
    o = match_options(opts, {"r": 1}) if comma else {}
    cns, b = s.e.get("Cns"), s.e.get("b")
    if not isinstance(b, Matrix):
        raise StataError(301, "last estimates not found")
    texts = constraint_lines(b, cns.data) if isinstance(cns, Matrix) else []
    if o.get("r"):
        s.r = {"k": float(len(texts))}
        for k, t in enumerate(texts, 1):
            s.r[f"cns{k}"] = t
        return
    if not texts:
        s.output.write("(no constraints)\n", "text")              # VERIFY
        return
    _display_cns(s, texts)


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

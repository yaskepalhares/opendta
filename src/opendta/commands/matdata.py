"""mkmat and svmat ([P] matrix mkmat): variables to matrices and back."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.dataset import Variable
from ..core.errors import StataError
from ..core.varlist import expand
from ..lang.syntax import match_options, parse_standard
from ..lang.words import strip_outer_quotes
from ._util import touse
from .matrix import Matrix, _store, get_matrix
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def _name_values(s: "Session", name: str, rows: np.ndarray, nchar: int) -> list[str]:
    var = s.data.get(expand(s.data, name)[0])
    vals = var.data[rows]
    if var.is_string:
        out = [str(v) for v in vals]
    else:
        out = [str(int(v)) if v == int(v) else repr(float(v)) if v < M.SYSMISS else "." for v in vals]
    return [x.replace(" ", "_")[:nchar] or "_" for x in out]


@command("mkmat")
def cmd_mkmat(s: "Session", args: str) -> None:
    """mkmat varlist [if] [in] [, matrix() nomissing rownames() roweq() rowprefix() obs nchar()].
    VERIFY: the row names of string rownames() with spaces and of obs."""
    p = parse_standard(args)
    if not p.varlist.strip():
        raise StataError(100, "varlist required")   # compat 2005
    o = match_options(p.options, {"matrix": 3, "nomissing": 5, "rownames": 5, "roweq": 5,
                                  "rowprefix": 4, "obs": 3, "nchar": 2})
    ds = s.data
    names = expand(ds, p.varlist)
    for n in names:
        if ds.get(n).is_string:
            raise StataError(109, "type mismatch")
    rows = np.flatnonzero(touse(s, p))
    X = np.column_stack([ds.get(n).data[rows].astype(float) for n in names]) if len(rows) \
        else np.zeros((0, len(names)))
    if o.get("nomissing") and len(rows):
        keep = ~(X >= M.SYSMISS).any(axis=1)
        X, rows = X[keep], rows[keep]
    if len(rows) == 0:
        raise StataError(2000, "no observations")   # VERIFY
    nchar = int(o["nchar"]) if o.get("nchar") else 32
    if o.get("rownames"):
        rnames = _name_values(s, str(o["rownames"]), rows, nchar)
    elif o.get("obs"):
        rnames = [str(i + 1) for i in rows]
    else:
        # rowprefix() replaces the r of r1, r2... (compat 2008: p1 p2)
        pre = strip_outer_quotes(str(o["rowprefix"])) if o.get("rowprefix") else "r"
        rnames = [f"{pre}{i}"[:nchar] for i in range(1, len(rows) + 1)]
    if o.get("rowprefix") and (o.get("rownames") or o.get("obs")):
        pre = strip_outer_quotes(str(o["rowprefix"]))
        rnames = [(pre + r)[:nchar] for r in rnames]   # VERIFY
    req = _name_values(s, str(o["roweq"]), rows, nchar) if o.get("roweq") else []
    store = _store(s)
    if o.get("matrix"):
        store[str(o["matrix"]).strip()] = Matrix(X, rnames, list(names), req)
    else:
        for j, n in enumerate(names):
            store[n] = Matrix(X[:, [j]], rnames, [n], req)


@command("svmat")
def cmd_svmat(s: "Session", args: str) -> None:
    """svmat [type] A [, names(col|eqcol|matcol|string)]: one variable per column of A."""
    p = parse_standard(args)
    words = p.varlist.split()
    vtype = "float"
    if words and words[0] in ("byte", "int", "long", "float", "double"):
        vtype = words.pop(0)
    if len(words) > 1:
        raise StataError(198, "invalid syntax")
    name = words[0] if words else ""
    if not name:
        raise StataError(111, "matrix  not found")   # compat 2005
    o = match_options(p.options, {"names": 1})
    A = get_matrix(s, name)
    how = strip_outer_quotes(str(o["names"])) if o.get("names") else ""
    if how == "col":
        vnames = list(A.colnames)
    elif how == "eqcol":
        # compat 2008: e1u, Du (the parts joined with nothing between)
        vnames = [(eq if eq != "_" else "") + c for eq, c in zip(A.coleq, A.colnames)]
    elif how == "matcol":
        vnames = [name + c for c in A.colnames]
    else:
        stub = how or name
        vnames = [f"{stub}{j}" for j in range(1, A.cols + 1)]
    ds = s.data
    for v in vnames:
        if ds.has(v):
            raise StataError(110, f"{v} already defined")
    if A.rows > ds.nobs:
        ds.set_obs(A.rows)
    for j, v in enumerate(vnames):
        col = np.full(ds.nobs, M.SYSMISS)
        col[:A.rows] = A.data[:, j]
        ds.add(Variable(v, vtype, col))
    ds.changed = True
    s.notify_state()

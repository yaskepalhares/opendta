"""tabdisp ([R] tabdisp): displays the data in memory as a table.

    tabdisp rowvar [colvar [supercolvar]] [if] [in], cellvar(varnames) [by(superrowvars) format(%fmt)
            center left concise missing totals dotz cellwidth(#) csepwidth(#) scsepwidth(#)
            stubwidth(#)]

It calculates nothing: each cell is the value of the cell variables in the first observation of
that combination of row, column, supercolumn and superrows. The layout follows `table`.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..core.layout import render_tabdisp
from ..lang.syntax import match_options, parse_standard
from ._util import touse
from .registry import command
from .summarize import by_groups
from .tabulate import _categories, _codes

if TYPE_CHECKING:
    from ..session import Session

_OPTIONS = {"cellvar": 1, "by": 2, "format": 3, "center": 3, "left": 4, "concise": 3, "missing": 4,
            "totals": 3, "dotz": 4, "cellwidth": 5, "csepwidth": 5, "scsepwidth": 6,
            "stubwidth": 5}


def _number(s: "Session", o: dict, name: str, default: int | None = None) -> int | None:
    if o.get(name) in (None, True):
        return default
    x = float(s.eval(str(o[name])))
    if x != int(x) or x < 0:
        raise StataError(198, f"option {name}() incorrectly specified")
    return int(x)


def _classes(s: "Session", v, mask: np.ndarray, o: dict):
    """Categories of a classifying variable: labels, codes per observation and the position of the
    category that stands for the total."""
    cats, labels, _ = _categories(s, v, mask, missing=True, nolabel=False)
    dotz = bool(o.get("dotz"))
    total = M.EXTENDED["z"] if dotz else M.SYSMISS
    if not v.is_string:
        if dotz:
            # . and .z swap roles: .z sorts last and . is an ordinary missing value
            key = lambda c: (c == M.EXTENDED["z"], c if c != M.EXTENDED["z"] else 0.0)  # noqa: E731
            order = sorted(range(len(cats)), key=lambda k: key(cats[k]))
            cats, labels = [cats[k] for k in order], [labels[k] for k in order]
        if o.get("totals"):
            labels = ["Total" if c == total else lab for c, lab in zip(cats, labels)]
    return cats, labels, _codes(v, cats)


@command("tabdisp", byable=True)
def cmd_tabdisp(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, _OPTIONS) if p.options.strip() else {}
    names = s.expand_varlist(p.varlist)
    if not 1 <= len(names) <= 3:
        raise StataError(198, "invalid syntax")
    if o.get("cellvar") in (None, True):
        raise StataError(198, "cellvar() required")
    cellvars = s.expand_varlist(str(o["cellvar"]))
    if not 1 <= len(cellvars) <= 5:
        typed = re.search(r"(?<![\w])(c(?:e(?:l(?:l(?:v(?:a(?:r)?)?)?)?)?)?)\s*\(", p.options)
        raise StataError(198, "too many variables specified\noption "
                              f"{typed.group(1) if typed else 'cellvar'}() incorrectly specified")
    byvars = s.expand_varlist(str(o["by"])) if o.get("by") not in (None, True) else []
    if len(byvars) > 4:
        raise StataError(198, "option by() incorrectly specified")
    fmt = str(o.get("format", "%10.0g")).strip()   # VERIFY: the default (1.235e+08, compat 1801)
    if not re.fullmatch(r"%[-0,~+]*\d*(\.\d+)?[a-zA-Z]+", fmt):
        raise StataError(198, "option format() incorrectly specified")
    ds = s.data
    base = touse(s, p)
    if not base.any():
        raise StataError(2000, "no observations")
    for sel, _ in by_groups(s):
        _table(s, o, names, cellvars, byvars, fmt, base & sel)


def _cell(s: "Session", name: str, i: int, fmt: str, show_missing: bool) -> str:
    v = s.data.get(name)
    x = v.raw[i] if v.is_string else float(v.data[i])
    if v.is_string:
        return str(x)
    if x >= M.SYSMISS:
        return M.missing_name(x) if show_missing else ""
    return format_value(x, fmt, pad=False).strip()


def _table(s: "Session", o: dict, names: list[str], cellvars: list[str], byvars: list[str],
           fmt: str, mask: np.ndarray) -> None:
    ds, out = s.data, s.output
    n = ds.nobs
    rv = ds.get(names[0])
    cv = ds.get(names[1]) if len(names) > 1 else None
    sv = ds.get(names[2]) if len(names) > 2 else None
    bvs = [ds.get(b) for b in byvars]
    rcats, rlabels, rc = _classes(s, rv, mask, o)
    ccats, clabels, cc = _classes(s, cv, mask, o) if cv is not None else ([None], [""], np.zeros(n, dtype=np.int64))
    scats, slabels, sc = _classes(s, sv, mask, o) if sv is not None else ([None], [""], np.zeros(n, dtype=np.int64))
    bclasses = [_classes(s, b, mask, o) for b in bvs]
    bc = [np.asarray(c[2]) for c in bclasses]
    # first observation of each combination
    first: dict[tuple, int] = {}
    for i in np.flatnonzero(mask):
        key = (tuple(int(c[i]) for c in bc), int(sc[i]), int(rc[i]), int(cc[i]))
        first.setdefault(key, int(i))
    show_missing = bool(o.get("missing"))
    ncell = len(cellvars)
    # the columns: with no column variable, the cell variables are the columns
    if cv is None:
        columns = [(0, k) for k in range(ncell)]
        heads = [ds.get(c).label or c for c in cellvars]
    else:
        columns = [(c, None) for c in range(len(ccats))]
        heads = list(clabels)
    groups = list(range(len(scats)))
    # superrow combinations present, in order
    bkeys = sorted({k[0] for k in first})

    def text_of(key: tuple, line: int, col: int, cellk: int | None) -> str:
        i = first.get(key)
        if i is None:
            return ""
        return _cell(s, cellvars[line if cellk is None else cellk], i, fmt, show_missing)

    # the body: a heading for each superrow combination and the rows with their lines of cells
    body: list[tuple] = []
    multi_line = cv is not None and ncell > 1
    for bk in bkeys:
        if bvs:
            title = " ".join(bclasses[j][1][bk[j]] for j in range(len(bvs)))
            body.append(("group", title))
        for r in range(len(rcats)):
            grid: list[list[str]] = []   # line -> cells across supercolumns and columns
            for line in range(ncell if multi_line else 1):
                row = []
                for g in groups:
                    for col, cellk in columns:
                        key = (bk, g, r, col)
                        row.append(text_of(key, line, col, cellk))
                grid.append(row)
            if o.get("concise") and not any(any(x for x in row) for row in grid):
                continue
            if multi_line and body and body[-1][0] == "row":
                body.append(("blank",))
            body.append(("row", rlabels[r], grid))
    cw = _number(s, o, "cellwidth")
    # the widest value of the cell variables in the whole sample, not only in the cells shown
    # (VERIFY: compat 1801)
    contents = [_cell(s, cvn, int(i), fmt, show_missing)
                for cvn in cellvars for i in np.flatnonzero(mask)]
    longest = max([1] + [len(x) for x in contents + heads])
    if cw is not None:
        width = cw
    else:
        # VERIFY: the minimum width of a cell
        width = max(longest, 4 if cv is not None else 10)
    csep = _number(s, o, "csepwidth", 2)
    scsep = _number(s, o, "scsepwidth", 4)
    # titles of the stub: the labels of the by variables and of the row variable
    rtitle = " and ".join([b.label or b.name for b in bvs] + [rv.label or rv.name])
    stub = max([9] + [len(item[1]) for item in body if item[0] != "blank"])
    stub = _number(s, o, "stubwidth", stub)
    group_title = (cv.label or cv.name) if cv is not None else ""
    super_title = " and ".join([sv.label or sv.name, group_title]) if sv is not None else ""
    out.write("\n", "text")
    for line in render_tabdisp(stub_title=rtitle, body=body, heads=heads,
                               groups=slabels if sv is not None else None,
                               group_title=group_title, super_title=super_title, width=width,
                               csep=csep, scsep=scsep, stub=stub, center=bool(o.get("center")),
                               left=bool(o.get("left")),
                               block=0 if any(ds.get(c).is_string for c in cellvars) else max([0] + [len(x) for x in contents])):
        out.write(line + "\n", "text")

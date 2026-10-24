"""matlist ([P] matlist): displays a matrix as a table.

    matlist matrix_exp [, style_options general_options]
    matlist matrix_exp, cspec(cspec) rspec(rspec) [general_options]

The layout is in core/layout.py. Equation names: the column equations go in a line above the
column names (showcoleq()), and each row equation in a line of its own above its rows. VERIFY the
whole layout with equations against Stata (compat 1910).
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value, parse_format
from ..core.layout import Field, Sep, align, render_table
from ..lang.syntax import find_top, match_options
from ..lang.words import strip_outer_quotes
from .matrix import Matrix, evaluate_matrix
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

_OPTIONS = {"lines": 3, "border": 3, "format": 3, "twidth": 2, "left": 3, "right": 3,
            "title": 3, "tindent": 3, "rowtitle": 4, "names": 5, "nonames": 3, "noblank": 3,
            "nohalf": 3, "nodotz": 3, "underscore": 5, "aligncolnames": 5, "linesize": 5,
            "cspec": 5, "rspec": 5, "showcoleq": 5, "roweqonly": 5, "coleqonly": 5,
            "colorcoleq": 6, "keepcoleq": 5}
_SHOWCOLEQ = {"first": "R", "each": "R", "combined": "C", "lcombined": "L", "rcombined": "R"}
_LINES = ("oneline", "eq", "rowtotal", "coltotal", "rctotal", "rows", "columns", "cells", "none")


def _int(s: "Session", text, name: str, low: int = 0) -> int:
    try:
        x = float(s.eval(str(text)))
    except StataError:
        raise StataError(198, f"option {name}() incorrectly specified")
    if x != int(x) or x < low:
        raise StataError(198, f"option {name}() incorrectly specified")
    return int(x)


def _border(text: str) -> set[str]:
    words = set(text.split()) if text else {"all"}
    allowed = {"none", "all", "rows", "columns", "left", "right", "top", "bottom"}
    if not words <= allowed:
        raise StataError(198, "option border() incorrectly specified")
    out: set[str] = set()
    for w in words:
        out |= {"all": {"left", "right", "top", "bottom"}, "rows": {"top", "bottom"},
                "columns": {"left", "right"}, "none": set()}.get(w, {w})
    return out


def _parse_cspec(text: str, nfields: int) -> tuple[list[Field], list[Sep], list[str]]:
    """cspec(): separators (& or |, with o# spaces before and after) and the specs of the fields.
    Returns the fields, the separators and the formats of the fields."""
    tokens = re.findall(r"o\d+|&|\||%[-~]?\d*(?:\.\d+)?[a-zA-Z]+|w\d+|[sbitecLRC]", text)
    seps: list[dict] = [{"line": None, "before": None, "after": None}]
    fields: list[dict] = []
    cur_field: dict | None = None
    for tok in tokens:
        sep = seps[-1]
        if tok[0] == "o" and tok[1:].isdigit():
            if cur_field is not None and cur_field.get("fmt"):
                cur_field = None
                seps.append({"line": None, "before": None, "after": None})
                sep = seps[-1]
            n = int(tok[1:])
            if sep["line"] is None:
                sep["before"] = n
            else:
                sep["after"] = n
        elif tok in "&|":
            if cur_field is not None and cur_field.get("fmt"):
                cur_field = None
                seps.append({"line": None, "before": None, "after": None})
                sep = seps[-1]
            sep["line"] = tok == "|"
        else:
            if cur_field is None or cur_field.get("fmt"):
                cur_field = {"fmt": "", "align": "", "width": None}
                fields.append(cur_field)
            if tok.startswith("%"):
                cur_field["fmt"] = tok
            elif tok in ("L", "R", "C"):
                cur_field["align"] = tok
            elif tok.startswith("w"):
                cur_field["width"] = int(tok[1:])
    if len(fields) != nfields or len(seps) != nfields + 1:
        raise StataError(198, "option cspec() incorrectly specified")      # VERIFY
    out_fields, out_seps, fmts = [], [], []
    for k, f in enumerate(fields):
        fmt = f["fmt"] or "%9.0g"
        w = f["width"] or parse_format(fmt).width
        out_fields.append(Field(w, f["align"] or ("L" if "-" in fmt else "R")))
        fmts.append(fmt)
    last = len(seps) - 1
    for k, sp in enumerate(seps):
        line = bool(sp["line"])
        # one space on each side of a line, none before the first separator or after the last
        before = sp["before"] if sp["before"] is not None else (0 if k == 0 or (k == last and not line) else 1)
        after = sp["after"] if sp["after"] is not None else (0 if k == last or (k == 0 and not line) else 1)
        out_seps.append(Sep(line, before, after))
    return out_fields, out_seps, fmts


@command("matlist")
def cmd_matlist(s: "Session", args: str) -> None:
    i = find_top(args, ",")
    head, opts = (args, "") if i < 0 else (args[:i], args[i + 1:])
    if not head.strip():
        raise StataError(198, "invalid syntax")
    mat: Matrix = evaluate_matrix(s, head.strip())
    o = match_options(opts, _OPTIONS) if opts.strip() else {}
    showcoleq = str(o.get("showcoleq", "combined")).strip() if o.get("showcoleq") not in (None, True) \
        else "combined"
    if showcoleq not in _SHOWCOLEQ:
        raise StataError(198, "option showcoleq() incorrectly specified")         # VERIFY
    if o.get("colorcoleq") not in (None, True) and str(o["colorcoleq"]).strip() not in ("txt", "res"):
        raise StataError(198, "option colorcoleq() incorrectly specified")        # VERIFY
    fmt = str(o.get("format", "%9.0g")).strip()
    twidth = _int(s, o["twidth"], "twidth", 1) if o.get("twidth") not in (None, True) else 12
    left = _int(s, o["left"], "left") if o.get("left") not in (None, True) else 0
    right = _int(s, o["right"], "right") if o.get("right") not in (None, True) else 0
    tindent = _int(s, o["tindent"], "tindent") if o.get("tindent") not in (None, True) else 0
    title = strip_outer_quotes(str(o["title"]).strip()) if o.get("title") not in (None, True) else ""
    rowtitle = strip_outer_quotes(str(o["rowtitle"]).strip()) if o.get("rowtitle") not in (None, True) else ""
    names = str(o.get("names", "all")).strip()
    if o.get("nonames"):
        names = "none"
    if names not in ("rows", "columns", "all", "none"):
        raise StataError(198, "option names() incorrectly specified")
    show_rows, show_cols = names in ("rows", "all"), names in ("columns", "all")
    if rowtitle and not (show_rows and show_cols):
        raise StataError(198, "option rowtitle() requires both row and column names")
    aname = str(o.get("aligncolnames", "ralign")).strip()
    how = {"ralign": "R", "lalign": "L", "center": "C"}.get(aname)
    if how is None:
        raise StataError(198, f"option {aname} not allowed")
    lines = str(o.get("lines", "oneline")).strip()
    if lines not in _LINES:
        raise StataError(198, "option lines() incorrectly specified")
    border = _border(str(o["border"]).strip()) if o.get("border") else set()
    linesize = _int(s, o["linesize"], "linesize", 1) if o.get("linesize") not in (None, True) \
        else int(float(s.settings.get("linesize", 80)))

    rown = [r.replace("_", " ") for r in mat.rownames] if o.get("underscore") else list(mat.rownames)
    coln = [c.replace("_", " ") for c in mat.colnames] if o.get("underscore") else list(mat.colnames)
    req = [q if q != "_" else "" for q in (mat.roweq or ["_"] * mat.rows)]
    ceq = [q if q != "_" else "" for q in (mat.coleq or ["_"] * mat.cols)]
    has_req, has_ceq = any(req), any(ceq)
    sym = mat.rows == mat.cols and mat.rows > 1 and np.array_equal(mat.data, mat.data.T)
    half = sym and not o.get("nohalf")

    def text(x: float, f: str) -> str:
        if x >= M.SYSMISS:
            return "" if o.get("nodotz") and x == M.EXTENDED["z"] else M.missing_name(float(x))
        return format_value(float(x), f, pad=False).strip()

    if o.get("cspec") or o.get("rspec"):
        if not (o.get("cspec") and o.get("rspec")):
            raise StataError(198, "options cspec() and rspec() must be specified together")
        fields, seps, fmts = _parse_cspec(str(o["cspec"]), mat.cols + (1 if show_rows else 0))
        numfmts = fmts[1:] if show_rows else fmts
        numfields = fields[1:] if show_rows else fields
    else:
        numfmts = [fmt] * mat.cols
        wfmt = parse_format(fmt).width
        numfields = [Field(max(wfmt, len(c)) if show_cols else wfmt, "R") for c in coln]
        fields = ([Field(twidth, "R")] if show_rows else []) + numfields
        n = len(fields)
        seps = [Sep() for _ in range(n + 1)]
        seps[n] = Sep(False, 0, 1)        # a blank after the last column (compat 1801)
        for k in range(1, n):
            seps[k] = Sep(False, 1, 1)
        vline = {"oneline", "eq", "rowtotal", "coltotal", "rctotal", "rows", "columns", "cells"}
        if show_rows and lines in vline and mat.cols:
            seps[1].line = True
        if lines in ("columns", "cells"):
            for k in range(1 if show_rows else 0, n):
                if k > 0:
                    seps[k].line = True
        if lines in ("coltotal", "rctotal") and mat.cols > 1:
            seps[n - 1].line = True
        if lines == "eq" and has_ceq:
            off = 1 if show_rows else 0
            for c in range(1, mat.cols):
                if ceq[c] != ceq[c - 1]:
                    seps[off + c].line = True
        if "left" in border:
            seps[0] = Sep(True, 0, 1)
        if "right" in border:
            seps[n] = Sep(True, 1, 0)

    rows: list[list[str]] = []
    spans: dict[int, list[tuple[int, int, str, str]]] = {}
    off = 1 if show_rows else 0
    eq_header = show_cols and has_ceq
    if eq_header:
        # the column equations: over each group of columns, or over its first column or each one
        cells = [""] * (off + mat.cols)
        start = 0
        for c in range(1, mat.cols + 1):
            if c == mat.cols or ceq[c] != ceq[start]:
                if showcoleq == "each":
                    for k in range(start, c):
                        cells[off + k] = ceq[k]
                elif showcoleq == "first":
                    cells[off + start] = ceq[start]
                else:
                    spans.setdefault(0, []).append((off + start, off + c - 1, ceq[start],
                                                    _SHOWCOLEQ[showcoleq]))
                start = c
        rows.append(cells)
    if show_cols and not (o.get("coleqonly") and has_ceq):
        header = [rowtitle] if show_rows else []
        rows.append(header + [c for c in coln])
    nh = len(rows)
    eq_rows: list[int] = []
    for r in range(mat.rows):
        if show_rows and has_req and not o.get("roweqonly") and (r == 0 or req[r] != req[r - 1]):
            eq_rows.append(len(rows))
            rows.append([req[r].ljust(fields[0].width)[:fields[0].width]] + [""] * mat.cols)
        stub = req[r] if o.get("roweqonly") and has_req else rown[r]
        cells = [stub[:fields[0].width]] if show_rows else []
        for c in range(mat.cols):
            cells.append("" if half and c > r else text(mat.data[r, c], numfmts[c]))
        rows.append(cells)
    # widths of the columns from the cspec/format are kept; names longer than the field are cut
    rules = [False] * (len(rows) + 1)
    if o.get("rspec"):
        spec = re.sub(r"\s", "", str(o["rspec"]))
        if len(spec) != len(rows) + 1 or set(spec) - set("-|&"):
            raise StataError(198, "option rspec() incorrectly specified")     # VERIFY
        rules = [ch in "-|" for ch in spec]
    else:
        if "top" in border:
            rules[0] = True
        if "bottom" in border:
            rules[-1] = True
        if lines in ("oneline", "eq", "rowtotal", "coltotal", "rctotal", "columns", "cells",
                     "rows") and nh and len(rows) > nh:
            rules[nh] = True
        if lines == "eq":
            for k in eq_rows[1:]:
                rules[k] = True
        if lines in ("rows", "cells"):
            for k in range(1 + nh, len(rows)):
                rules[k] = True
        if lines in ("rowtotal", "rctotal") and len(rows) - nh > 1:
            rules[len(rows) - 1] = True
    # column names are aligned as asked (default: right)
    if show_cols:
        for k in range(mat.cols):
            rows[nh - 1][off + k] = align(rows[nh - 1][off + k], fields[off + k].width, how)
    out = s.output
    if not o.get("noblank"):
        out.write("\n", "text")
    if title:
        out.write(" " * tindent + title + "\n", "text")
    for line in render_table(rows, fields, seps, rules, left=left,
                             limit=linesize - right, repeat=1 if show_rows else 0, spans=spans,
                             groups=[None] * off + ceq if o.get("keepcoleq") else None):
        out.write(line + "\n", "text")

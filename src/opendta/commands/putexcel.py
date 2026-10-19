"""putexcel ([P] putexcel): writes text, expressions and matrices to
cells of an .xlsx spreadsheet (openpyxl) or of an .xls one (xlrd and xlwt, io/excel.py).

Accepts both Stata 14 syntaxes:
    putexcel set file [, sheet(name) replace modify]
    putexcel A1 = "text" B2 = (expression) C3 = matrix(name), names
and the 14.0 one, with the file in the command itself:
    putexcel A1=("text") B1=(expr) using file [, sheet() replace modify]
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..lang.syntax import find_top, match_options, parse_options
from ..lang.words import strip_outer_quotes
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

_CELL = re.compile(r"^([A-Za-z]{1,3})(\d+)$")


def _state(s: "Session") -> dict:
    st = getattr(s, "_putexcel", None)
    if st is None:
        st = {}
        s._putexcel = st
    return st


def _path(text: str) -> Path:
    t = strip_outer_quotes(text.strip())
    if not t:
        raise StataError(198, "invalid file specification")
    p = Path(t).expanduser()
    if p.suffix == "":
        p = p.with_suffix(".xlsx")
    return p


def _split_cells(text: str) -> list[tuple[str, str]]:
    """'A1 = "x" B2 = (1+1) C3=matrix(m), names' -> [(cell, value)]."""
    items: list[tuple[str, str]] = []
    i, n = 0, len(text)
    while i < n:
        while i < n and text[i].isspace():
            i += 1
        if i >= n:
            break
        m = re.match(r"([A-Za-z]{1,3}\d+)\s*=\s*", text[i:])
        if not m:
            raise StataError(198, "invalid syntax")
        cell = m.group(1).upper()
        i += m.end()
        j = i
        depth = 0
        while j < n:
            c = text[j]
            if c == '"':
                k = text.find('"', j + 1)
                j = n if k == -1 else k + 1
                continue
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
            elif c.isspace() and depth == 0:
                rest = text[j:].lstrip()
                if re.match(r"[A-Za-z]{1,3}\d+\s*=", rest):
                    break
            j += 1
        items.append((cell, text[i:j].strip()))
        i = j
    return items


def _col_index(letters: str) -> int:
    k = 0
    for ch in letters.upper():
        k = k * 26 + (ord(ch) - 64)
    return k


@command("putexcel", "putexcel")
def cmd_putexcel(s: "Session", args: str) -> None:
    t = args.strip()
    st = _state(s)
    head = t.split(None, 1)[0] if t else ""
    if head == "set":
        rest = t[3:].strip()
        c = find_top(rest, ",")
        body, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
        o = match_options(opt, {"sheet": 2, "replace": 3, "modify": 3}) if opt.strip() else {}
        st.clear()
        st["path"] = _path(body)
        st["sheet"] = strip_outer_quotes(str(o["sheet"]).split(",")[0]) if isinstance(o.get("sheet"), str) else ""
        st["sheet_replace"] = isinstance(o.get("sheet"), str) and "replace" in str(o["sheet"])
        st["mode"] = "replace" if o.get("replace") else "modify"
        st["fresh"] = bool(o.get("replace"))
        return
    if head == "clear":
        st.clear()
        return
    if head == "describe":
        out = s.output
        if not st:
            out.write("  (no file set)\n", "text")   # VERIFY
            return
        out.write(f"  File: {st['path']}\n  Sheet: {st.get('sheet') or 'Sheet1'}\n", "text")
        return
    # 14.0 syntax: ... using file
    m = re.search(r"\busing\b", t)
    c = find_top(t, ",")
    body, opt = (t, "") if c == -1 else (t[:c], t[c + 1:])
    o = {}
    for nm, a in parse_options(opt) if opt.strip() else []:
        o[nm] = a if a is not None else True
    if m and (c == -1 or m.start() < c):
        cells_txt = body[:m.start()]
        path = _path(body[m.end():])
        sheet = strip_outer_quotes(str(o.get("sheet", "")).split(",")[0]) if isinstance(o.get("sheet"), str) else ""
        fresh = "replace" in o
    else:
        if opt.strip():
            # with the putexcel set syntax, Stata 14.0 does not accept options at the end of the
            # command (matrix options go inside matrix(name, names)); compat 0804: "option names
            # not allowed"
            raise StataError(198, f"option {next(iter(o))} not allowed")
        if not st:
            raise StataError(198, "file not set; use putexcel set")   # VERIFY
        cells_txt = body
        path, sheet = st["path"], st.get("sheet", "")
        fresh = st.get("fresh", False)
    items = _split_cells(cells_txt)
    if not items:
        raise StataError(198, "invalid syntax")
    if path.suffix.lower() == ".xls":
        from ..io.excel import XlsBook
        wb = XlsBook(path if path.exists() and not fresh else None)
    else:
        import openpyxl
        if path.exists() and not fresh:
            wb = openpyxl.load_workbook(path)
        else:
            wb = openpyxl.Workbook()
            if sheet:
                wb.active.title = sheet
    ws = wb[sheet] if sheet and sheet in wb.sheetnames else (wb.create_sheet(sheet) if sheet else wb.active)
    for cell, val in items:
        mm = _CELL.match(cell)
        col, row = _col_index(mm.group(1)), int(mm.group(2))
        _write_value(s, ws, row, col, val, o)
    try:
        wb.save(path)
    except OSError:
        raise StataError(603, f"file {path} could not be opened")
    if st:
        st["fresh"] = False    # after the first write, the following ones modify
    s.output.write(f"file {path} saved\n", "text")


def _write_value(s: "Session", ws, row: int, col: int, val: str, o: dict) -> None:
    v = val.strip()
    mm = re.match(r"^matrix\s*\((.*)\)\s*(,.*)?$", v, re.S)
    if mm:
        inner = mm.group(1)
        c = find_top(inner, ",")
        name, sub = (inner, "") if c == -1 else (inner[:c], inner[c + 1:])
        flags = {n for n, _ in parse_options(sub)} if sub.strip() else set()
        if mm.group(2):
            flags |= {n for n, _ in parse_options(mm.group(2)[1:])}
        flags |= {k for k in o if k in ("names", "rownames", "colnames")}
        _write_matrix(s, ws, row, col, name.strip(), flags)
        return
    if v.startswith('"') or v.startswith('`"'):
        ws.cell(row=row, column=col, value=strip_outer_quotes(v))
        return
    if v.startswith("(") and v.endswith(")"):
        v = v[1:-1]
    if v.startswith('"'):
        ws.cell(row=row, column=col, value=strip_outer_quotes(v))
        return
    res = s.eval(v)
    if isinstance(res, str):
        ws.cell(row=row, column=col, value=res)
    else:
        from ..core.missing import SYSMISS
        ws.cell(row=row, column=col, value=None if res >= SYSMISS else float(res))


def _write_matrix(s: "Session", ws, row: int, col: int, name: str, flags: set) -> None:
    from .matrix import get_matrix
    mat = get_matrix(s, name)
    data = mat.data
    rn = "names" in flags or "rownames" in flags
    cn = "names" in flags or "colnames" in flags
    r0 = row + (1 if cn else 0)
    c0 = col + (1 if rn else 0)
    if cn:
        for j, nm in enumerate(mat.colnames):
            ws.cell(row=row, column=c0 + j, value=nm)
    if rn:
        for i, nm in enumerate(mat.rownames):
            ws.cell(row=r0 + i, column=col, value=nm)
    from ..core.missing import SYSMISS
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            x = float(data[i, j])
            ws.cell(row=r0 + i, column=c0 + j, value=None if x >= SYSMISS else x)

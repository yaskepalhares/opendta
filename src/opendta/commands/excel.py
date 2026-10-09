"""import excel e export excel ([D] import excel).

Mensagens e regras não documentadas estão marcadas VERIFICAR, com casos em
compat/do/0106_excel.do.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from ..core.errors import StataError
from ..core.varlist import expand, unique
from ..io.excel import ExcelOptions, ExcelWriteOptions, describe_excel, read_excel, write_excel
from ..lang.syntax import Parsed, match_options, parse_standard
from ..lang.words import strip_outer_quotes
from ._util import touse

if TYPE_CHECKING:
    from ..session import Session


def excel_path(text: str, *, must_exist: bool) -> Path:
    t = strip_outer_quotes(text.strip())
    if not t:
        raise StataError(198, "invalid file specification")
    p = Path(t).expanduser()
    if p.suffix == "":
        # VERIFICAR: extensão padrão (o Stata 14 usa .xls; o OpenDTA grava .xlsx)
        candidates = [p.with_suffix(".xlsx"), p.with_suffix(".xls")]
        if must_exist:
            for c in candidates:
                if c.exists():
                    return c
        return candidates[0]
    return p


def _sheet_arg(arg: str) -> tuple[str, str]:
    """'"Plan 1", modify' -> ('Plan 1', 'modify')."""
    text = arg.strip()
    mode = ""
    if text.startswith('"'):
        end = text.find('"', 1)
        name, rest = text[1:end], text[end + 1:]
    else:
        name, _, rest = text.partition(",")
        name = name.strip()
        rest = "," + rest if _ else ""
    rest = rest.strip().lstrip(",").strip()
    if rest:
        if rest not in ("modify", "replace"):
            raise StataError(198, f"option sheet() invalid: {rest}")
        mode = rest
    return name, mode


def import_excel(s: "Session", args: str) -> None:
    p = parse_standard(args)
    if p.using is not None:
        file = p.using       # varlist antes de using (colunas): VERIFICAR, ainda ignorado
    else:
        file = p.varlist
    o = match_options(p.options, {"sheet": 5, "cellrange": 9, "firstrow": 8, "case": 4,
                                  "allstring": 9, "clear": 5, "describe": 8, "locale": 6})
    path = excel_path(file, must_exist=True)
    if not path.exists():
        raise StataError(601, f"file {path} not found")
    if o.get("describe"):
        rows = describe_excel(path)
        out = s.output
        width = max([len(n) + 4 for n, _ in rows] + [10])
        out.write(f"\n{'Sheet':>{width}} | Range\n", "text")
        out.write(f"  {'-' * (width - 1)}+{'-' * 9}\n", "text")
        for name, rng in rows:
            out.write(f"{name:>{width}} | {rng}\n", "text")
        # gravados do último para o primeiro: return list mostra worksheet_1 antes
        s.r = {}
        for k, (name, rng) in reversed(list(enumerate(rows, start=1))):
            s.r[f"range_{k}"] = rng
            s.r[f"worksheet_{k}"] = name
        s.r["N_worksheet"] = float(len(rows))
        return
    ds = s.data
    if ds.changed and ds.nvars and not o.get("clear"):
        raise StataError(4, "no; data in memory would be lost")
    opts = ExcelOptions()
    if o.get("sheet"):
        opts.sheet = strip_outer_quotes(str(o["sheet"]).strip())
    if o.get("cellrange"):
        opts.cellrange = str(o["cellrange"]).strip()
    opts.firstrow = bool(o.get("firstrow"))
    if o.get("case"):
        c = str(o["case"]).strip()
        if c not in ("preserve", "lower", "upper"):
            raise StataError(198, "option case() invalid")
        opts.case = c
    opts.allstring = bool(o.get("allstring"))
    new = read_excel(path, opts)
    s.data = new
    s.notify_state()


def export_excel(s: "Session", args: str) -> None:
    p = parse_standard(args)
    if p.using is not None:
        names_text, file = p.varlist, p.using
    else:
        names_text, file = "", p.varlist
    o = match_options(p.options, {"sheet": 5, "cell": 4, "firstrow": 8, "nolabel": 5,
                                  "replace": 7, "sheetmodify": 11, "sheetreplace": 12,
                                  "datestring": 10, "missing": 7, "locale": 6})
    ds = s.data
    if ds.nvars == 0:
        raise StataError(102, "no variables defined")   # VERIFICAR
    names = unique(expand(ds, names_text)) if names_text.strip() else list(ds.names)
    path = excel_path(file, must_exist=False)
    opts = ExcelWriteOptions()
    sheet_mode = ""
    if o.get("sheet"):
        opts.sheet, sheet_mode = _sheet_arg(str(o["sheet"]))
    if o.get("cell"):
        opts.cell = str(o["cell"]).strip()
    if o.get("firstrow"):
        fr = str(o["firstrow"]).strip()
        if fr not in ("variables", "varlabels"):
            raise StataError(198, "option firstrow() must be variables or varlabels")
        opts.firstrow = fr
    opts.nolabel = bool(o.get("nolabel"))
    if o.get("sheetmodify") or sheet_mode == "modify":
        opts.mode = "modify"
    elif o.get("sheetreplace") or sheet_mode == "replace":
        opts.mode = "sheetreplace"
    elif o.get("replace"):
        opts.mode = "replace"
    if path.exists() and opts.mode == "new":
        raise StataError(602, f"file {path} already exists")
    rows = np.flatnonzero(touse(s, Parsed(if_=p.if_, in_=p.in_)))
    try:
        write_excel(ds, path, names, rows, opts)
    except OSError as e:
        raise StataError(603, f"file {path} could not be opened ({e.strerror})")
    s.output.write(f"file {path} saved\n", "text")

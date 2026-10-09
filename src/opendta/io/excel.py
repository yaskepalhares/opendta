"""Planilhas do Excel: leitura (.xlsx, e .xls com xlrd) e gravação (.xlsx).

Usado por import excel / export excel. As regras seguem o manual
[D] import excel do Stata 14; os pontos não documentados estão marcados com
VERIFICAR e têm casos em compat/do/0106_excel.do.
"""

from __future__ import annotations

import datetime as _dt
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..core import missing as M
from ..core import storage
from ..core.dataset import (NAME_RE, RESERVED, Dataset, Variable, smallest_type_for,
                            str_len, str_type_for)
from ..core.errors import StataError

EPOCH = _dt.datetime(1960, 1, 1)
_CELL = re.compile(r"^([A-Za-z]{1,3})(\d+)$")


def col_letters(n: int) -> str:
    """1 -> A, 27 -> AA."""
    s = ""
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def col_number(letters: str) -> int:
    n = 0
    for ch in letters.upper():
        n = n * 26 + ord(ch) - 64
    return n


def parse_cell(text: str) -> tuple[int, int]:
    """'B3' -> (linha 3, coluna 2)."""
    m = _CELL.match(text.strip())
    if not m:
        raise StataError(198, f"{text} invalid cell name")
    return int(m.group(2)), col_number(m.group(1))


def parse_range(text: str) -> tuple[tuple[int, int] | None, tuple[int, int] | None]:
    """'A2:C10', 'A2', ':C10', 'A2:' -> ((linha, coluna) | None, (linha, coluna) | None)."""
    a, colon, b = text.strip().partition(":")
    start = parse_cell(a) if a.strip() else None
    end = parse_cell(b) if colon and b.strip() else None
    return start, end


# ---------------------------------------------------------------------------
# leitura
# ---------------------------------------------------------------------------

@dataclass
class ExcelOptions:
    sheet: str | None = None
    cellrange: str | None = None
    firstrow: bool = False
    case: str = "preserve"          # VERIFICAR: padrão do import excel
    allstring: bool = False


def _open(path: Path):
    """Devolve (nomes das planilhas, função que lê uma planilha em linhas)."""
    suffix = path.suffix.lower()
    if suffix == ".xls":
        try:
            import xlrd
        except ImportError:
            raise StataError(603, "reading .xls files requires the xlrd package (pip install xlrd)")
        book = xlrd.open_workbook(str(path))

        def rows_xls(name: str):
            sh = book.sheet_by_name(name)
            for r in range(sh.nrows):
                out = []
                for c in range(sh.ncols):
                    cell = sh.cell(r, c)
                    if cell.ctype == xlrd.XL_CELL_DATE:
                        out.append(xlrd.xldate.xldate_as_datetime(cell.value, book.datemode))
                    elif cell.ctype in (xlrd.XL_CELL_EMPTY, xlrd.XL_CELL_BLANK):
                        out.append(None)
                    else:
                        out.append(cell.value)
                yield out
        return book.sheet_names(), rows_xls, None
    import openpyxl
    try:
        book = openpyxl.load_workbook(str(path), read_only=True, data_only=True)
    except Exception as e:  # noqa: BLE001
        raise StataError(603, f"file {path} could not be read as an Excel file ({e})")

    def rows_xlsx(name: str):
        for row in book[name].iter_rows(values_only=True):
            yield list(row)
    return book.sheetnames, rows_xlsx, book


def _trim(rows: list[list]) -> list[list]:
    """Remove linhas e colunas vazias no fim."""
    while rows and all(v is None or v == "" for v in rows[-1]):
        rows.pop()
    width = 0
    for r in rows:
        for j in range(len(r) - 1, -1, -1):
            if r[j] is not None and r[j] != "":
                width = max(width, j + 1)
                break
    return [list(r[:width]) + [None] * (width - len(r[:width])) for r in rows]


def describe_excel(path: Path) -> list[tuple[str, str]]:
    names, reader, book = _open(path)
    out = []
    for name in names:
        rows = _trim(list(reader(name)))
        if rows and rows[0]:
            out.append((name, f"A1:{col_letters(len(rows[0]))}{len(rows)}"))
        else:
            out.append((name, ""))
    if book is not None:
        book.close()
    return out


def _cell_text(v) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "TRUE" if v else "FALSE"
    if isinstance(v, float) and v.is_integer() and abs(v) < 1e15:
        return str(int(v))
    if isinstance(v, _dt.datetime):
        # allstring: a data como o Excel a mostra (m/d/aaaa; observado no Stata 14)
        if not (v.hour or v.minute or v.second or v.microsecond):
            return f"{v.month}/{v.day}/{v.year}"
        return f"{v.month}/{v.day}/{v.year} {v.hour}:{v.minute:02d}:{v.second:02d}"   # VERIFICAR
    if isinstance(v, _dt.date):
        return f"{v.month}/{v.day}/{v.year}"
    return str(v)


def _make_name(raw: str, case: str) -> str | None:
    t = raw.strip()
    if case == "lower":
        t = t.lower()
    elif case == "upper":
        t = t.upper()
    # caracteres inválidos são removidos ("Renda mensal" → Rendamensal,
    # observado no Stata 14)
    t = re.sub(r"[^A-Za-z0-9_]", "", t)[:32]
    if not t or not NAME_RE.match(t) or t in RESERVED:
        return None
    return t


def read_excel(path: str | Path, opts: ExcelOptions | None = None) -> Dataset:
    opts = opts or ExcelOptions()
    path = Path(path)
    names, reader, book = _open(path)
    sheet = opts.sheet if opts.sheet is not None else names[0]
    if sheet not in names:
        raise StataError(601, f"worksheet {sheet} not found")
    rows = list(reader(sheet))
    if book is not None:
        book.close()

    if opts.cellrange:
        start, end = parse_range(opts.cellrange)
        r0, c0 = start or (1, 1)
        r1, c1 = end or (len(rows), max((len(r) for r in rows), default=0))
        rows = [r[c0 - 1:c1] + [None] * max(0, c1 - max(len(r), c0 - 1)) for r in rows[r0 - 1:r1]]
        first_col = c0
    else:
        first_col = 1
    rows = _trim(rows)
    if not opts.cellrange:
        # VERIFICAR: linhas vazias no início também são descartadas
        while rows and all(v is None or v == "" for v in rows[0]):
            rows.pop(0)

    header = None
    if opts.firstrow and rows:
        header = rows[0]
        rows = rows[1:]
    ncols = len(rows[0]) if rows else (len(header) if header else 0)

    ds = Dataset()
    ds.nobs = len(rows)
    used: set[str] = set()
    for j in range(ncols):
        letters = col_letters(first_col + j)
        raw_name = _cell_text(header[j]) if header is not None and j < len(header) else ""
        name = _make_name(raw_name, opts.case) if raw_name.strip() else None
        if name is None or name in used:
            name = letters
        used.add(name)
        # VERIFICAR: o cabeçalho vira rótulo da variável
        label = raw_name.strip() if header is not None else ""
        col = [r[j] if j < len(r) else None for r in rows]
        filled = [v for v in col if v is not None and v != ""]
        is_num = (not opts.allstring and bool(filled)
                  and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in filled))
        is_date = (not opts.allstring and bool(filled)
                   and all(isinstance(v, (_dt.datetime, _dt.date)) for v in filled))
        if is_date:
            has_time = any(isinstance(v, _dt.datetime) and (v.hour or v.minute or v.second
                                                             or v.microsecond) for v in filled)
            vals = np.full(len(col), M.SYSMISS)
            for i, v in enumerate(col):
                if v is None or v == "":
                    continue
                d = v if isinstance(v, _dt.datetime) else _dt.datetime(v.year, v.month, v.day)
                delta = d - EPOCH
                vals[i] = (delta.days * 86_400_000 + delta.seconds * 1000
                           + delta.microseconds // 1000) if has_time else delta.days
            vtype = "double" if has_time else smallest_type_for(vals)
            # VERIFICAR: formato de data-hora (o de data foi observado no Stata 14)
            var = Variable(name, vtype, vals,
                           fmt="%tcnn/dd/CCYY_hh:MM:SS" if has_time else "%tdnn/dd/CCYY")
        elif is_num:
            vals = np.array([M.SYSMISS if v is None or v == "" else float(v) for v in col])
            vtype = smallest_type_for(vals)
            if vtype == "float":
                vtype = "double"        # VERIFICAR: decimais do Excel ficam double
            var = Variable(name, vtype, vals)
        else:
            # datas como texto ficam alinhadas à direita em 10 posições: o
            # Stata cria str10 e lista "1/1/2020" numa coluna de 10 (VERIFICAR
            # se o espaço fica mesmo no texto)
            texts = [_cell_text(v).rjust(10) if isinstance(v, (_dt.datetime, _dt.date)) else _cell_text(v)
                     for v in col]
            width = max((str_len(s) for s in texts), default=1) or 1
            var = Variable(name, str_type_for(width), texts)
        var.label = label
        ds.vars.append(var)
    ds.changed = True
    return ds


# ---------------------------------------------------------------------------
# gravação
# ---------------------------------------------------------------------------

@dataclass
class ExcelWriteOptions:
    sheet: str = "Sheet1"
    mode: str = "new"            # new | replace (arquivo) | modify | sheetreplace
    cell: str = "A1"
    firstrow: str | None = None  # None | variables | varlabels
    nolabel: bool = False


def _date_kind(v: Variable) -> str:
    fmt = v.fmt.lstrip("%-")
    if fmt.startswith(("td", "d")):
        return "date"
    if fmt.startswith(("tc", "tC")):
        return "datetime"
    return ""


def _excel_values(ds: Dataset, v: Variable, idx: np.ndarray, nolabel: bool) -> list:
    """Valores de uma coluna para a planilha (None = célula vazia)."""
    if v.is_string:
        return [s if s != "" else None for s in v.raw[idx].tolist()]
    dec = storage.decode(v.raw[idx], v.vtype)
    miss = (dec >= M.SYSMISS).tolist()
    vals = dec.tolist()
    kind = _date_kind(v)
    if kind == "date":
        return [None if m else (EPOCH + _dt.timedelta(days=x)).date() for x, m in zip(vals, miss)]
    if kind == "datetime":
        return [None if m else EPOCH + _dt.timedelta(milliseconds=x) for x, m in zip(vals, miss)]
    lab = ds.value_labels.get(v.value_label, {}) if v.value_label and not nolabel else {}
    out = []
    for x, m in zip(vals, miss):
        if m:
            out.append(None)
        elif x.is_integer() and abs(x) < 1e15:
            k = int(x)
            out.append(lab.get(k, k))
        else:
            out.append(x)
    return out


_NUMBER_FORMAT = {"date": "m/d/yyyy", "datetime": "m/d/yyyy h:mm:ss"}  # VERIFICAR


def write_excel(ds: Dataset, path: str | Path, names: list[str], rows: np.ndarray,
                opts: ExcelWriteOptions | None = None) -> None:
    """Arquivo novo: modo write_only do openpyxl (rápido, linha a linha).
    sheetmodify/sheetreplace: abre a pasta de trabalho existente."""
    import openpyxl
    from openpyxl.cell import WriteOnlyCell

    opts = opts or ExcelWriteOptions()
    path = Path(path)
    if path.suffix.lower() == ".xls":
        raise StataError(198, "export to .xls is not supported by OpenDTA; use .xlsx")
    vars_ = [ds.get(n) for n in names]
    kinds = [_date_kind(v) for v in vars_]
    r0, c0 = parse_cell(opts.cell)
    header = None
    if opts.firstrow:
        header = [(v.label or v.name) if opts.firstrow == "varlabels" else v.name for v in vars_]
    rows = np.asarray(rows)

    def data_rows():
        step = 50_000
        for start in range(0, len(rows), step):
            idx = rows[start:start + step]
            cols = [_excel_values(ds, v, idx, opts.nolabel) for v in vars_]
            yield from zip(*cols)

    if path.exists() and opts.mode in ("modify", "sheetreplace"):
        book = openpyxl.load_workbook(str(path))
        if opts.sheet in book.sheetnames and opts.mode == "sheetreplace":
            pos = book.sheetnames.index(opts.sheet)
            del book[opts.sheet]
            ws = book.create_sheet(opts.sheet, pos)
        elif opts.sheet in book.sheetnames:
            ws = book[opts.sheet]
        else:
            ws = book.create_sheet(opts.sheet)
        r = r0
        if header:
            for j, text in enumerate(header):
                ws.cell(row=r, column=c0 + j, value=text)
            r += 1
        for values in data_rows():
            for j, x in enumerate(values):
                if x is not None:
                    c = ws.cell(row=r, column=c0 + j, value=x)
                    if kinds[j]:
                        c.number_format = _NUMBER_FORMAT[kinds[j]]
            r += 1
    else:
        book = openpyxl.Workbook(write_only=True)
        ws = book.create_sheet(opts.sheet)
        pad = [None] * (c0 - 1)
        for _ in range(r0 - 1):
            ws.append([])
        if header:
            ws.append(pad + header)
        for values in data_rows():
            row = list(pad)
            for j, x in enumerate(values):
                if x is not None and kinds[j]:
                    cell = WriteOnlyCell(ws, value=x)
                    cell.number_format = _NUMBER_FORMAT[kinds[j]]
                    row.append(cell)
                else:
                    row.append(x)
            ws.append(row)
    tmp = path.with_name(path.name + ".opendta-tmp")
    try:
        book.save(str(tmp))
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    tmp.replace(path)

"""Arquivos de texto delimitado (CSV, TSV): leitura e gravação.

Usado por import delimited / export delimited e pelos antigos insheet /
outsheet. As regras de inferência seguem o manual [D] import delimited do
Stata 14; os pontos não documentados estão marcados com VERIFICAR e têm
casos em compat/do/0105_delimitado.do.
"""

from __future__ import annotations

import csv
import io
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..core import missing as M
from ..core.dataset import (NAME_RE, RESERVED, Dataset, Variable, smallest_type_for,
                            str_len, str_type_for)

_NUM_RE = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?$")


@dataclass
class ReadOptions:
    delimiter: str | None = None        # None: detecta vírgula ou tabulação
    varnames: int | None = None         # linha dos nomes (1-based); 0 = nonames; None = detectar
    case: str = "lower"                 # lower | upper | preserve
    asdouble: bool = False
    stringcols: set[int] = field(default_factory=set)    # 1-based
    numericcols: set[int] = field(default_factory=set)
    rowrange: tuple[int | None, int | None] = (None, None)
    colrange: tuple[int | None, int | None] = (None, None)
    encoding: str | None = None
    bindquote: str = "loose"            # loose | strict | nobind
    stripquotes: str = "default"        # default | yes | no


def _decode(raw: bytes, encoding: str | None) -> str:
    if encoding:
        return raw.decode(encoding)
    if raw.startswith(b"\xef\xbb\xbf"):
        return raw[3:].decode("utf-8")
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


def _sniff(text: str) -> str:
    # VERIFICAR: o Stata detecta vírgula ou tabulação pela primeira linha
    first = text.split("\n", 1)[0]
    return "\t" if first.count("\t") > first.count(",") else ","


def split_rows(text: str, delimiter: str, bindquote: str = "loose") -> list[list[str]]:
    if bindquote == "nobind":
        return [line.split(delimiter) for line in text.splitlines()]
    reader = csv.reader(io.StringIO(text, newline=""), delimiter=delimiter,
                        quotechar='"', doublequote=True, strict=False)
    return [row for row in reader]


def _is_number(s: str) -> bool:
    return bool(_NUM_RE.match(s))


def _to_number(s: str) -> float:
    s = s.strip()
    if s in ("", "."):
        return M.SYSMISS
    if len(s) == 2 and s[0] == "." and s[1].isalpha() and s[1].islower():
        return M.missing_code(s)
    return float(s.replace("d", "e").replace("D", "e"))


def _numeric_cell(s: str) -> bool:
    s = s.strip()
    return s in ("", ".") or _is_number(s) or (len(s) == 2 and s[0] == "." and "a" <= s[1] <= "z")


def make_name(raw: str, case: str) -> str | None:
    """Nome de variável a partir do cabeçalho; None se não der um nome válido."""
    t = raw.strip()
    if case == "lower":
        t = t.lower()
    elif case == "upper":
        t = t.upper()
    # VERIFICAR: troca de caracteres inválidos por _ e corte em 32
    t = re.sub(r"[^A-Za-z0-9_]", "_", t)[:32]
    if not t or not NAME_RE.match(t) or t in RESERVED:
        return None
    return t


def read_delimited(path: str | Path, opts: ReadOptions | None = None) -> Dataset:
    opts = opts or ReadOptions()
    text = _decode(Path(path).read_bytes(), opts.encoding)
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    delim = opts.delimiter or _sniff(text)
    rows = split_rows(text, delim, opts.bindquote)
    while rows and all(c.strip() == "" for c in rows[-1]):
        rows.pop()
    ncols = max((len(r) for r in rows), default=0)
    rows = [r + [""] * (ncols - len(r)) for r in rows]

    c0, c1 = opts.colrange
    c0 = (c0 or 1) - 1
    c1 = c1 or ncols
    rows = [r[c0:c1] for r in rows]
    ncols = max(0, min(c1, ncols) - c0)

    # linha de nomes
    header: list[str] | None = None
    if opts.varnames is None:
        # VERIFICAR: heurística de detecção dos nomes na 1ª linha (aqui: todas
        # as células preenchidas da 1ª linha são texto não numérico)
        if rows:
            filled = [c for c in rows[0] if c.strip()]
            if filled and not any(_numeric_cell(c) for c in filled):
                header = rows[0]
                rows = rows[1:]
    elif opts.varnames > 0:
        k = opts.varnames - 1
        header = rows[k] if k < len(rows) else None
        rows = rows[k + 1:]

    r0, r1 = opts.rowrange
    if r0 or r1:
        # VERIFICAR: aqui rowrange conta as linhas de dados (depois do cabeçalho)
        rows = rows[(r0 or 1) - 1:(r1 or len(rows))]

    ds = Dataset()
    ds.nobs = len(rows)
    used: set[str] = set()
    for j in range(ncols):
        col = [r[j] for r in rows]
        raw_name = header[j] if header is not None and j < len(header) else ""
        name = make_name(raw_name, opts.case) if raw_name.strip() else None
        if name is None or name in used:
            name = f"v{j + 1}"
        used.add(name)
        # VERIFICAR: o cabeçalho original vira rótulo quando o nome teve de mudar
        # (além da caixa)
        label = raw_name.strip() if raw_name.strip() and name.lower() != raw_name.strip().lower() else ""
        force_str = (j + 1) in opts.stringcols
        force_num = (j + 1) in opts.numericcols
        numeric = not force_str and (force_num or all(_numeric_cell(c) for c in col))
        if numeric:
            vals = np.array([_to_number(c) if _numeric_cell(c) else M.SYSMISS for c in col],
                            dtype=np.float64)
            vtype = smallest_type_for(vals)
            if vtype == "float" and opts.asdouble:
                vtype = "double"
            elif vtype == "double" and not opts.asdouble:
                nm = vals[vals < M.SYSMISS]
                # inteiros grandes ficam double; decimais viram float (VERIFICAR)
                if not (nm.size and np.all(nm == np.trunc(nm))):
                    vtype = "float"
            if vtype == "float":
                ok = vals < M.SYSMISS
                vals[ok] = vals[ok].astype(np.float32).astype(np.float64)
            var = Variable(name, vtype, vals)
        else:
            vals = np.array(col, dtype=object)
            width = max((str_len(s) for s in col), default=1)
            var = Variable(name, str_type_for(width), vals)
        var.label = label
        ds.vars.append(var)
    ds.changed = True
    return ds


# ---------------------------------------------------------------------------
# gravação
# ---------------------------------------------------------------------------

@dataclass
class WriteOptions:
    delimiter: str = ","
    varnames: bool = True
    nolabel: bool = False
    datafmt: bool = False
    quote: bool = False          # aspas em todas as strings
    leading_zero: bool = False   # VERIFICAR: o Stata grava .5 em vez de 0.5


def number_text(x: float, vtype: str, fmt: str | None = None, *, leading_zero: bool = False) -> str:
    if x >= M.SYSMISS:
        return "" if M.missing_name(x) == "." else M.missing_name(x)
    if fmt:
        from ..core.formats import format_value
        return format_value(x, fmt, pad=False).strip()
    if vtype in ("byte", "int", "long"):
        return str(int(x))
    v = np.float32(x) if vtype == "float" else np.float64(x)
    if v == 0 or 1e-5 <= abs(v) < 1e16:
        t = np.format_float_positional(v, trim="-")
    else:
        t = np.format_float_scientific(v, trim="-")
    if not leading_zero:
        if t.startswith("0."):
            t = t[1:]
        elif t.startswith("-0."):
            t = "-" + t[2:]
    return t


def write_delimited(ds: Dataset, path: str | Path, names: list[str], rows: np.ndarray,
                    opts: WriteOptions | None = None) -> None:
    opts = opts or WriteOptions()
    d = opts.delimiter
    vars_ = [ds.get(n) for n in names]

    def q(s: str) -> str:
        if opts.quote or d in s or '"' in s or "\n" in s:
            return '"' + s.replace('"', '""') + '"'
        return s

    lines: list[str] = []
    if opts.varnames:
        lines.append(d.join(v.name for v in vars_))
    for i in rows:
        cells = []
        for v in vars_:
            x = v.data[i]
            if v.is_string:
                cells.append(q(str(x)))
                continue
            x = float(x)
            if not opts.nolabel and v.value_label and x < M.SYSMISS and x == int(x):
                lab = ds.value_labels.get(v.value_label, {})
                if int(x) in lab:
                    cells.append(q(lab[int(x)]))
                    continue
            cells.append(number_text(x, v.vtype, v.fmt if opts.datafmt else None,
                                     leading_zero=opts.leading_zero))
        lines.append(d.join(cells))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")

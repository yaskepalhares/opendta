"""Arquivos de texto delimitado (CSV, TSV): leitura e gravação.

Usado por import delimited / export delimited e pelos antigos insheet /
outsheet. As regras de inferência seguem o manual [D] import delimited do
Stata 14; os pontos não documentados estão marcados com VERIFICAR e têm
casos em compat/do/0105_delimitado.do.
"""

from __future__ import annotations

import csv
import io
import itertools
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..core import missing as M
from ..core import storage
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
    # insheet: se a 1ª linha tem menos campos que os dados, ela não é cabeçalho
    header_must_cover: bool = False
    trim_cells: bool = False            # insheet: tira espaços das pontas dos campos


def _sniff(first_line: str) -> str:
    # VERIFICAR: o Stata detecta vírgula ou tabulação pela primeira linha
    return "\t" if first_line.count("\t") > first_line.count(",") else ","


def split_rows(text: str, delimiter: str, bindquote: str = "loose") -> list[list[str]]:
    """Divide um texto inteiro em linhas e células (usado em testes e trechos curtos)."""
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


CHUNK_ROWS = 50_000           # linhas convertidas por vez
_STR_CACHE_MAX = 50_000       # valores distintos guardados por coluna de texto


class _Column:
    """Coluna em construção. Começa numérica; vira texto na primeira célula
    que não é número (a não ser com numericcols)."""

    def __init__(self, index: int, offset: int, force_str: bool, force_num: bool):
        self.index = index
        self.offset = offset             # linhas anteriores ao aparecimento da coluna
        self.numeric = not force_str
        self.force_num = force_num
        self.num_chunks: list[np.ndarray] = []
        self.text: list[str] = []
        self.maxlen = 1
        self.cache: dict[str, str] | None = {}
        self.flipped_late = False        # virou texto depois do 1º bloco: precisa reler

    def add_text(self, cells) -> None:
        cache = self.cache
        if cache is None:
            self.text.extend(cells)
            longest = max((str_len(s) for s in cells), default=0)
        else:
            longest = 0
            for s in cells:
                hit = cache.get(s)
                if hit is None:
                    hit = cache[s] = s
                    n = str_len(s)
                    if n > longest:
                        longest = n
                self.text.append(hit)
            if len(cache) > _STR_CACHE_MAX:
                self.cache = None
        self.maxlen = max(self.maxlen, longest)

    def add(self, cells: list[str]) -> None:
        if not self.numeric:
            self.add_text(cells)
            return
        vals = _parse_numbers(cells, strict=not self.force_num)
        if vals is None:                 # apareceu texto
            self.numeric = False
            if self.num_chunks or self.offset:
                self.flipped_late = True
                self.num_chunks = []
            else:
                self.add_text(cells)
            return
        self.num_chunks.append(vals)


def _parse_numbers(cells: list[str], strict: bool) -> np.ndarray | None:
    """Converte um bloco de células em float64. None se houver texto (strict)."""
    if not cells:
        return np.empty(0)
    s = np.array(cells, dtype=str)
    s = np.char.strip(s)
    miss = (s == "") | (s == ".")
    # missing estendidos .a–.z
    ext = (np.char.str_len(s) == 2) & np.char.startswith(s, ".") & ~miss
    if ext.any():
        ext &= np.char.isalpha(np.char.lstrip(s, "."))
    try:
        if (np.char.find(s, "_") >= 0).any():
            raise ValueError
        vals = np.where(miss | ext, "0", s).astype(np.float64)
        if not np.isfinite(vals).all():
            raise ValueError
        vals[miss] = M.SYSMISS
        if ext.any():
            vals[ext] = [M.missing_code(x) for x in s[ext].tolist()]
        return vals
    except (ValueError, KeyError):
        pass
    # caminho lento: missing estendidos (.a), expoente com d, texto
    out = np.empty(len(cells))
    for i, c in enumerate(cells):
        if _numeric_cell(c):
            out[i] = _to_number(c)
        elif strict:
            return None
        else:
            out[i] = M.SYSMISS
    return out


def _open_rows(path: Path, opts: ReadOptions, encoding: str):
    """Arquivo aberto e iterador de linhas (listas de células)."""
    f = open(path, encoding=encoding, newline="")
    first = f.readline()
    f.seek(0)
    delim = opts.delimiter or _sniff(first)
    if opts.bindquote == "nobind":
        rows = (line.rstrip("\r\n").split(delim) for line in f)
    else:
        rows = csv.reader(f, delimiter=delim, quotechar='"', doublequote=True, strict=False)
    if opts.trim_cells:
        # insheet tira os espaços das pontas de cada campo (observado no Stata 14)
        rows = ([c.strip() for c in r] for r in rows)
    return f, rows


def _data_rows(rows, opts: ReadOptions, header_out: list):
    """Aplica cabeçalho, colrange e rowrange; descarta linhas em branco no fim.
    O cabeçalho encontrado é posto em header_out[0]."""
    c0, c1 = opts.colrange
    c0 = (c0 or 1) - 1
    r0, r1 = opts.rowrange
    it = iter(rows)

    def cut(r: list[str]) -> list[str]:
        return r[c0:c1] if (c0 or c1) else r

    header_out.append(None)
    lead: list[list[str]] = []
    if opts.varnames is None:
        # VERIFICAR: heurística de detecção dos nomes na 1ª linha (aqui: todas
        # as células preenchidas da 1ª linha são texto não numérico)
        first = next(it, None)
        if first is None:
            return
        first = cut(first)
        filled = [c for c in first if c.strip()]
        if filled and not any(_numeric_cell(c) for c in filled):
            header_out[0] = first
        else:
            lead = [first]
    elif opts.varnames > 0:
        for _ in range(opts.varnames - 1):
            next(it, None)
        h = next(it, None)
        header_out[0] = cut(h) if h is not None else None

    k = 0          # linha de dados (1-based), para o rowrange
    pending: list[list[str]] = []
    for r in itertools.chain(lead, (cut(x) for x in it)):
        if not ((r and r[0].strip()) or any(c.strip() for c in r)):
            pending.append(r)
            continue
        for b in pending + [r]:
            k += 1
            # VERIFICAR: aqui rowrange conta as linhas de dados (depois do cabeçalho)
            if r0 and k < r0:
                continue
            if r1 and k > r1:
                return
            yield b
        pending = []


def read_delimited(path: str | Path, opts: ReadOptions | None = None) -> Dataset:
    opts = opts or ReadOptions()
    encodings = [opts.encoding] if opts.encoding else ["utf-8-sig", "latin-1"]
    last: Exception | None = None
    for enc in encodings:
        try:
            return _read(Path(path), opts, enc)
        except UnicodeDecodeError as e:
            last = e
    raise last  # type: ignore[misc]


def _read(path: Path, opts: ReadOptions, encoding: str) -> Dataset:
    header_box: list = []
    cols: list[_Column] = []
    nrows = 0
    f, rows = _open_rows(path, opts, encoding)
    with f:
        buf: list[list[str]] = []

        def flush() -> None:
            nonlocal nrows
            if not buf:
                return
            width = max(len(r) for r in buf)
            while len(cols) < width:
                j = len(cols)
                cols.append(_Column(j, nrows, (j + 1) in opts.stringcols, (j + 1) in opts.numericcols))
            padded = [r if len(r) == width else r + [""] * (width - len(r)) for r in buf]
            for c, cells in zip(cols, zip(*padded)):
                c.add(list(cells))
            nrows += len(buf)
            buf.clear()

        for r in _data_rows(rows, opts, header_box):
            buf.append(r)
            if len(buf) >= CHUNK_ROWS:
                flush()
        flush()

    header = header_box[0] if header_box else None
    if header is not None and opts.header_must_cover and len(header) < len(cols):
        from dataclasses import replace
        return _read(path, replace(opts, varnames=0, header_must_cover=False), encoding)
    if header is not None:
        while len(cols) < len(header):
            j = len(cols)
            cols.append(_Column(j, nrows, (j + 1) in opts.stringcols, (j + 1) in opts.numericcols))

    # colunas que viraram texto depois do 1º bloco: relê só elas
    late = [c for c in cols if c.flipped_late]
    if late:
        for c in late:
            c.text, c.cache, c.maxlen, c.offset = [], {}, 1, 0
        f, rows = _open_rows(path, opts, encoding)
        with f:
            chunk: list[list[str]] = []
            for r in _data_rows(rows, opts, []):
                chunk.append(r)
                if len(chunk) >= CHUNK_ROWS:
                    for c in late:
                        c.add_text([x[c.index] if c.index < len(x) else "" for x in chunk])
                    chunk = []
            for c in late:
                c.add_text([x[c.index] if c.index < len(x) else "" for x in chunk])

    ds = Dataset()
    ds.nobs = nrows
    used: set[str] = set()
    for c in cols:
        j = c.index
        raw_name = header[j] if header is not None and j < len(header) else ""
        name = make_name(raw_name, opts.case) if raw_name.strip() else None
        if name is None or name in used:
            name = f"v{j + 1}"
        used.add(name)
        # VERIFICAR: o cabeçalho original vira rótulo quando o nome teve de mudar
        # (além da caixa)
        label = raw_name.strip() if raw_name.strip() and name.lower() != raw_name.strip().lower() else ""
        if c.numeric:
            parts = ([np.full(c.offset, M.SYSMISS)] if c.offset else []) + c.num_chunks
            vals = np.concatenate(parts) if parts else np.full(nrows, M.SYSMISS)
            c.num_chunks = []
            vtype = smallest_type_for(vals)
            if vtype == "float" and opts.asdouble:
                vtype = "double"
            elif vtype == "double" and not opts.asdouble:
                nm = vals[vals < M.SYSMISS]
                # inteiros grandes ficam double; decimais viram float (VERIFICAR)
                if not (nm.size and np.all(nm == np.trunc(nm))):
                    vtype = "float"
            var = Variable(name, vtype, vals)       # float: arredonda para precisão simples
            del vals
        else:
            text = ([""] * c.offset if c.offset else []) + c.text
            c.text = []
            var = Variable(name, str_type_for(c.maxlen), text)
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
    never_quote: bool = False    # outsheet, noquote: texto cru, mesmo com o delimitador


def number_text(x: float, vtype: str, fmt: str | None = None, *, leading_zero: bool = False) -> str:
    """Texto de um valor, igual ao que export delimited grava."""
    if fmt and x < M.SYSMISS:
        from ..core.formats import format_value, plain_numbers
        with plain_numbers():
            return format_value(x, fmt, pad=False).strip()
    v = Variable("_", vtype, np.array([x], dtype=np.float64))
    return str(_number_cells(v, v.raw, leading_zero)[0])


def _number_cells(v: Variable, raw: np.ndarray, leading_zero: bool) -> np.ndarray:
    """Texto de cada valor numérico (vetor de str), sem laço em Python."""
    dec = storage.decode(raw, v.vtype)
    miss = dec >= M.SYSMISS
    out = np.empty(len(dec), dtype=object)
    if v.vtype in ("byte", "int", "long"):
        txt = raw.astype(np.int64).astype(str)
    else:
        vals = raw if v.vtype == "float" else dec
        whole = ~miss & (dec == np.trunc(dec)) & (np.abs(dec) < 1e15)
        txt = np.where(miss, 0, vals).astype(str)
        if whole.any():
            txt = txt.astype(object)
            txt[whole] = dec[whole].astype(np.int64).astype(str)
            txt = txt.astype(str)
        if not leading_zero:
            m = np.char.startswith(txt, "0.")
            if m.any():
                txt[m] = np.char.replace(txt[m], "0.", ".", count=1)
            m = np.char.startswith(txt, "-0.")
            if m.any():
                txt[m] = np.char.replace(txt[m], "-0.", "-.", count=1)
    out[:] = txt
    if miss.any():
        out[miss] = [("" if M.missing_name(x) == "." else M.missing_name(x)) for x in dec[miss]]
    return out


def write_delimited(ds: Dataset, path: str | Path, names: list[str], rows: np.ndarray,
                    opts: WriteOptions | None = None) -> None:
    """Grava em blocos de linhas, montando cada coluna de uma vez."""
    opts = opts or WriteOptions()
    d = opts.delimiter
    vars_ = [ds.get(n) for n in names]
    special = (d, '"', "\n", "\r")

    def q(s: str) -> str:
        if opts.never_quote:
            return s
        if opts.quote or any(ch in s for ch in special):
            return '"' + s.replace('"', '""') + '"'
        return s

    rows = np.asarray(rows)
    with open(path, "w", encoding="utf-8", newline="") as f:
        if opts.varnames:
            f.write(d.join(v.name for v in vars_) + "\n")
        for start in range(0, len(rows), CHUNK_ROWS):
            idx = rows[start:start + CHUNK_ROWS]
            columns = []
            for v in vars_:
                part = v.raw[idx]
                if v.is_string:
                    columns.append([q(s) for s in part.tolist()])
                    continue
                if opts.datafmt:
                    from ..core.formats import format_value, plain_numbers
                    dec = storage.decode(part, v.vtype)
                    with plain_numbers():
                        cells = np.array([("" if x >= M.SYSMISS and M.missing_name(x) == "."
                                           else M.missing_name(x) if x >= M.SYSMISS
                                           else format_value(x, v.fmt, pad=False).strip()) for x in dec],
                                         dtype=object)
                else:
                    cells = _number_cells(v, part, opts.leading_zero)
                lab = ds.value_labels.get(v.value_label, {}) if v.value_label and not opts.nolabel else {}
                if lab:
                    dec = storage.decode(part, v.vtype)
                    for value, text in lab.items():
                        hit = dec == value
                        if hit.any():
                            cells[hit] = q(text)
                columns.append(cells.tolist())
            f.write("".join(d.join(r) + "\n" for r in zip(*columns)) if columns else "")

"""Texto em formato livre (infile) e em colunas fixas (infix e dicionários).

Formato livre: valores separados por espaços ou vírgulas, lidos em
sequência sem olhar as quebras de linha ([D] infile (free format)).

Colunas fixas: cada variável ocupa uma faixa de colunas numa das linhas da
observação ([D] infix, [D] infile (fixed format)). Quando todas as posições
são conhecidas de antemão, a leitura é vetorizada: o arquivo vira uma
matriz de bytes e cada variável é um recorte de colunas.

Pontos não documentados estão marcados com VERIFICAR e têm casos em
compat/do/0107_infile_infix.do.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..core import missing as M
from ..core.dataset import Dataset, Variable, str_type_for
from ..core.errors import StataError

_NUMERIC = ("byte", "int", "long", "float", "double")


def read_text(path: Path) -> str:
    raw = path.read_bytes()
    if raw.startswith(b"\xef\xbb\xbf"):
        raw = raw[3:]
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


def to_number(text: str) -> float | None:
    """Número do Stata a partir do texto; None se não for número."""
    t = text.strip()
    if t in ("", "."):
        return M.SYSMISS
    if len(t) == 2 and t[0] == "." and "a" <= t[1] <= "z":
        return M.EXTENDED[t[1]]
    try:
        x = float(t.replace("d", "e").replace("D", "e"))
    except ValueError:
        return None
    if x != x or x in (float("inf"), float("-inf")) or "_" in t:
        return None
    return x


# ---------------------------------------------------------------------------
# formato livre
# ---------------------------------------------------------------------------

def tokenize_free(text: str) -> list[str]:
    """Separa valores por espaços e vírgulas; aspas (simples ou duplas)
    delimitam strings; duas vírgulas seguidas indicam valor ausente."""
    out: list[str] = []
    i, n = 0, len(text)
    last_comma = False
    while i < n:
        ch = text[i]
        if ch in " \t\r\n":
            i += 1
            continue
        if ch == ",":
            if last_comma:
                out.append("")
            last_comma = True
            i += 1
            continue
        last_comma = False
        if ch in "\"'":
            j = text.find(ch, i + 1)
            if j == -1:
                j = n
            out.append(text[i + 1:j])
            i = j + 1
            continue
        j = i
        while j < n and text[j] not in " \t\r\n,":
            j += 1
        out.append(text[i:j])
        i = j
    return out


@dataclass
class FreeVar:
    name: str            # "" para _skip
    vtype: str
    skip: int = 0
    value_label: str = ""    # nome:rótulo; com automatic, textos viram códigos


def read_free(text: str, spec: list[FreeVar], *, automatic: bool = False
              ) -> tuple[Dataset, list[tuple[str, str, int]], bool]:
    """Devolve (dados, avisos (texto, variável, registro 0-based), terminou no
    meio de uma observação). `automatic` só vale para variáveis declaradas
    com :rótulo (observado no Stata 14)."""
    tokens = tokenize_free(text)
    width = sum(v.skip if not v.name else 1 for v in spec)
    if width == 0:
        raise StataError(198, "invalid syntax")
    nobs, rest = divmod(len(tokens), width)
    partial = rest != 0
    if partial:
        nobs += 1
        tokens += [""] * (width - rest)
    warnings: list[tuple[str, str, int]] = []
    ds = Dataset()
    ds.nobs = nobs
    pos = 0
    for v in spec:
        if not v.name:
            pos += v.skip
            continue
        col = tokens[pos::width][:nobs]
        pos += 1
        if v.vtype.startswith("str"):
            if v.vtype in ("str", "strL"):
                longest = max((len(s.encode("utf-8")) for s in col), default=1)
                vtype = "strL" if v.vtype == "strL" else str_type_for(longest)
            else:
                w = int(v.vtype[3:])
                col = [s.encode("utf-8")[:w].decode("utf-8", "ignore") for s in col]
                vtype = v.vtype
            ds.vars.append(Variable(v.name, vtype, col))
            continue
        vals = np.empty(nobs)
        labels: dict[str, int] = {}
        for i, s in enumerate(col):
            x = to_number(s)
            if x is None:
                if automatic and v.value_label:
                    x = labels.setdefault(s, len(labels) + 1)
                else:
                    warnings.append((s, v.name, i))
                    x = M.SYSMISS
            vals[i] = x
        var = Variable(v.name, v.vtype, vals)
        if v.value_label:
            if labels:
                ds.value_labels[v.value_label] = {code: text for text, code in labels.items()}
            var.value_label = v.value_label
        ds.vars.append(var)
    ds.changed = True
    return ds, warnings, partial


# ---------------------------------------------------------------------------
# colunas fixas
# ---------------------------------------------------------------------------

@dataclass
class FixedField:
    name: str
    vtype: str
    line: int                    # 1-based dentro da observação
    start: int                   # coluna inicial, 1-based
    width: int
    decimals: int = 0            # casas decimais implícitas (%#.#f)
    label: str = ""
    value_label: str = ""


@dataclass
class FixedLayout:
    fields: list[FixedField] = field(default_factory=list)
    lines: int = 1               # linhas por observação
    first: int = 1               # primeira linha de dados no arquivo


def read_fixed(text: str, layout: FixedLayout) -> Dataset:
    """Lê colunas fixas de forma vetorizada."""
    raw_lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    if raw_lines and raw_lines[-1] == "":
        raw_lines.pop()
    raw_lines = raw_lines[layout.first - 1:]
    L = max(1, layout.lines)
    nobs = len(raw_lines) // L
    # VERIFICAR: observação incompleta no fim é descartada
    raw_lines = raw_lines[:nobs * L]
    ds = Dataset()
    ds.nobs = nobs
    if nobs == 0:
        for f in layout.fields:
            ds.vars.append(Variable(f.name, f.vtype if f.vtype != "str" else "str1", []))
        return ds
    by_line: list[np.ndarray] = []
    for k in range(L):
        sel = raw_lines[k::L]
        enc = [s.encode("utf-8") for s in sel]
        need = max([f.start + f.width - 1 for f in layout.fields if f.line == k + 1] + [1])
        width = max(need, max((len(b) for b in enc), default=0))
        arr = np.array([b.ljust(width) for b in enc], dtype=f"S{width}")
        by_line.append(arr.view(np.uint8).reshape(nobs, width))
    for f in layout.fields:
        mat = by_line[f.line - 1]
        cut = np.ascontiguousarray(mat[:, f.start - 1:f.start - 1 + f.width])
        cells = cut.view(f"S{f.width}").reshape(nobs)
        if f.vtype.startswith("str"):
            texts = [b.decode("utf-8", "replace").rstrip() for b in cells.tolist()]
            texts = [t.lstrip() for t in texts]        # VERIFICAR: espaços à esquerda
            vtype = f.vtype
            if vtype == "str":
                vtype = str_type_for(max(1, f.width))
            var = Variable(f.name, vtype, texts)
        else:
            var = Variable(f.name, f.vtype, _numbers(cells, f.decimals))
        var.label = f.label
        var.value_label = f.value_label
        ds.vars.append(var)
    ds.changed = True
    return ds


def _numbers(cells: np.ndarray, decimals: int) -> np.ndarray:
    s = np.char.strip(cells)
    miss = (s == b"") | (s == b".")
    try:
        vals = np.where(miss, b"0", s).astype(np.float64)
        bad = ~np.isfinite(vals)
        if bad.any():
            raise ValueError
    except ValueError:
        out = np.empty(len(s))
        has_point = np.zeros(len(s), dtype=bool)
        for i, b in enumerate(s.tolist()):
            t = b.decode("latin-1")
            x = to_number(t)
            out[i] = M.SYSMISS if x is None else x   # VERIFICAR: texto vira missing
            has_point[i] = "." in t
        if decimals:
            scale = ~has_point & (out < M.SYSMISS)
            out[scale] = out[scale] / 10 ** decimals
        return out
    vals[miss] = M.SYSMISS
    if decimals:
        scale = ~miss & (np.char.find(s, b".") < 0)
        vals[scale] = vals[scale] / 10 ** decimals
    return vals


# ---------------------------------------------------------------------------
# especificação do infix: [#lines] [#:] [type] var start[-end] ...
# ---------------------------------------------------------------------------

_INFIX_TOKEN = re.compile(r"\d+:|\d+-\d+|\d+|/|\S+")


def parse_infix_spec(text: str) -> FixedLayout:
    layout = FixedLayout()
    toks = _INFIX_TOKEN.findall(text)
    i = 0
    line = 1
    pending_type = "float"
    while i < len(toks):
        t = toks[i]
        if t.isdigit() and i + 1 < len(toks) and toks[i + 1] in ("lines", "line"):
            layout.lines = int(t)
            i += 2
            continue
        if t.isdigit() and i + 1 < len(toks) and toks[i + 1] == "firstlineoffile":
            layout.first = int(t)
            i += 2
            continue
        if t.endswith(":") and t[:-1].isdigit():
            line = int(t[:-1])
            i += 1
            continue
        if t == "/":
            line += 1
            i += 1
            continue
        if t in _NUMERIC or t == "str" or re.fullmatch(r"str\d+", t):
            pending_type = t
            i += 1
            continue
        if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", t):
            names = [t]
            i += 1
            # varlist com uma única posição não é permitido; cada nome tem a sua
            if i >= len(toks):
                raise StataError(198, "invalid syntax")
            pos = toks[i]
            if pos.endswith(":") and pos[:-1].isdigit():
                line = int(pos[:-1])
                i += 1
                pos = toks[i]
            m = re.fullmatch(r"(\d+)(?:-(\d+))?", pos)
            if not m:
                raise StataError(198, f"invalid column specification for {names[0]}")
            a = int(m.group(1))
            b = int(m.group(2)) if m.group(2) else a
            if b < a:
                raise StataError(198, f"invalid column range {pos}")
            vtype = pending_type
            if vtype.startswith("str") and vtype != "str":
                vtype = vtype
            layout.fields.append(FixedField(names[0], vtype, line, a, b - a + 1))
            pending_type = "float"
            i += 1
            continue
        raise StataError(198, f"invalid syntax: {t}")
    layout.lines = max([layout.lines] + [f.line for f in layout.fields])
    return layout


# ---------------------------------------------------------------------------
# dicionários do infile e do infix
# ---------------------------------------------------------------------------

_FMT = re.compile(r"^%(\d+)(?:\.(\d+))?([fgesS])$")
_DIRECTIVE = re.compile(r"^_(column|skip|line|lines|newline|firstlineoffile|first)(?:\((\d+)\))?$")


@dataclass
class Dictionary:
    layout: FixedLayout
    datafile: str | None
    inline: str | None           # dados depois da chave de fechamento
    free: bool = False           # algum campo sem largura: leitura livre


def _strip_comments(line: str) -> str:
    s = line.strip()
    if s.startswith("*"):
        return ""
    if "//" in s:
        s = s.split("//", 1)[0]
    return s.strip()


def parse_dictionary(text: str, kind: str = "infile") -> Dictionary:
    """Lê um arquivo de dicionário. kind = infile | infix."""
    m = re.search(r"\b(?:infile\s+|infix\s+)?dictionary\b(.*?)\{", text, re.S)
    if not m:
        raise StataError(198, "dictionary invalid")
    head = m.group(1)
    um = re.search(r"\busing\s+(\"[^\"]+\"|\S+)", head)
    datafile = um.group(1).strip('"') if um else None
    close = text.find("}", m.end())
    if close == -1:
        raise StataError(198, "dictionary invalid: '}' not found")
    body = text[m.end():close]
    after = text[close + 1:]
    inline = after[1:] if after.startswith("\n") else after.lstrip("\r\n")
    lines = [_strip_comments(ln) for ln in body.splitlines()]
    if kind == "infix":
        layout = parse_infix_spec(" ".join(ln for ln in lines if ln))
        return Dictionary(layout, datafile, inline or None)

    layout = FixedLayout()
    col, line = 1, 1
    free = False
    for raw in lines:
        if not raw:
            continue
        toks = re.findall(r'"[^"]*"|\S+', raw)
        # diretivas no começo da linha (_column(5) long id %6f ...)
        while toks and _DIRECTIVE.match(toks[0]):
            d = _DIRECTIVE.match(toks.pop(0))
            name, arg = d.group(1), int(d.group(2)) if d.group(2) else None
            if name == "column":
                col = arg or 1
            elif name == "skip":
                col += arg if arg is not None else 1
            elif name == "line":
                line, col = arg or 1, 1
            elif name == "newline":
                line, col = line + (arg or 1), 1
            elif name == "lines":
                layout.lines = arg or 1
            elif name in ("firstlineoffile", "first"):
                layout.first = arg or 1
        if not toks:
            continue
        i = 0
        vtype = "float"
        if toks[0] in _NUMERIC or re.fullmatch(r"str\d+|strL?", toks[0]):
            vtype = toks[0]
            i = 1
        if i >= len(toks):
            raise StataError(198, f"dictionary invalid: {raw}")
        name = toks[i]
        i += 1
        value_label = ""
        if i < len(toks) and toks[i].startswith(":"):
            value_label = toks[i][1:]
            i += 1
        width, decimals = None, 0
        if i < len(toks) and toks[i].startswith("%"):
            fm = _FMT.match(toks[i])
            if not fm:
                raise StataError(198, f"invalid format {toks[i]}")
            width = int(fm.group(1))
            decimals = int(fm.group(2) or 0) if fm.group(3) == "f" else 0
            i += 1
        label = toks[i].strip('"') if i < len(toks) and toks[i].startswith('"') else ""
        if width is None:          # sem formato: leitura livre (str# só limita o tamanho)
            free = True
            width = 0
        layout.fields.append(FixedField(name, vtype, line, col, width, decimals, label, value_label))
        col += width
    layout.lines = max([layout.lines] + [f.line for f in layout.fields])
    return Dictionary(layout, datafile, inline or None, free)


def read_dictionary_free(text: str, layout: FixedLayout) -> Dataset:
    """Dicionário com campos sem largura: lê campo a campo (mais lento).
    Campos sem largura pegam o próximo valor separado por espaço a partir
    da coluna atual."""
    raw_lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    if raw_lines and raw_lines[-1] == "":
        raw_lines.pop()
    raw_lines = raw_lines[layout.first - 1:]
    L = max(1, layout.lines)
    nobs = len(raw_lines) // L
    cols: list[list] = [[] for _ in layout.fields]
    for o in range(nobs):
        obs_lines = raw_lines[o * L:(o + 1) * L]
        cursor = {}
        for j, f in enumerate(layout.fields):
            ln = obs_lines[f.line - 1] if f.line - 1 < len(obs_lines) else ""
            pos = max(f.start - 1, cursor.get(f.line, 0)) if f.width == 0 else f.start - 1
            if f.width:
                cell = ln[pos:pos + f.width]
                cursor[f.line] = pos + f.width
            else:
                k = pos
                while k < len(ln) and ln[k] in " \t,":
                    k += 1
                e = k
                if k < len(ln) and ln[k] == '"':
                    e = ln.find('"', k + 1)
                    e = len(ln) if e == -1 else e
                    cell = ln[k + 1:e]
                    e += 1
                else:
                    while e < len(ln) and ln[e] not in " \t,":
                        e += 1
                    cell = ln[k:e]
                cursor[f.line] = e
            cols[j].append(cell)
    ds = Dataset()
    ds.nobs = nobs
    for f, col in zip(layout.fields, cols):
        if f.vtype.startswith("str"):
            texts = [c.strip() for c in col]
            if f.vtype in ("str", "strL"):
                vtype = "strL" if f.vtype == "strL" else str_type_for(
                    max((len(t.encode()) for t in texts), default=1))
            else:
                vtype = f.vtype
                w = int(vtype[3:])
                texts = [x.encode("utf-8")[:w].decode("utf-8", "ignore") for x in texts]
            var = Variable(f.name, vtype, texts)
        else:
            vals = np.array([(lambda x: M.SYSMISS if x is None else x)(to_number(c)) for c in col])
            if f.decimals:
                scale = np.array(["." not in c for c in col]) & (vals < M.SYSMISS)
                vals[scale] = vals[scale] / 10 ** f.decimals
            var = Variable(f.name, f.vtype, vals)
        var.label = f.label
        var.value_label = f.value_label
        ds.vars.append(var)
    ds.changed = True
    return ds

"""Sintaxe padrão dos comandos (manual [U] 11):

    [by varlist:] command [varlist] [=exp] [if exp] [in range] [weight] [using filename] [, options]

`parse_standard` separa as partes respeitando aspas, parênteses e colchetes.
`parse_options` e `match_options` tratam as opções com abreviações.
`parse_in` converte um intervalo `in` em índices.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from ..core.errors import StataError

_WEIGHT = re.compile(r"^\[\s*(fweight|fw|frequency|aweight|aw|cellsize|pweight|pw|iweight|iw|weight)\s*=\s*(.+?)\s*\]$")
_WEIGHT_NAMES = {"fw": "fweight", "fweight": "fweight", "frequency": "fweight",
                 "aw": "aweight", "aweight": "aweight", "cellsize": "aweight",
                 "pw": "pweight", "pweight": "pweight", "iw": "iweight", "iweight": "iweight",
                 "weight": "weight"}


@dataclass
class Parsed:
    varlist: str = ""
    exp: str | None = None
    if_: str | None = None
    in_: str | None = None
    weight: tuple[str, str] | None = None
    using: str | None = None
    options: str = ""
    has_comma: bool = False


def _scan(text: str):
    """Gera (posição, caractere, profundidade) fora de aspas."""
    depth = 0
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if text.startswith('`"', i):
            d, j = 1, i + 2
            while j < n and d:
                if text.startswith('`"', j):
                    d += 1
                    j += 2
                elif text.startswith("\"'", j):
                    d -= 1
                    j += 2
                else:
                    j += 1
            i = j
            continue
        if c == '"':
            j = text.find('"', i + 1)
            i = n if j == -1 else j + 1
            continue
        if c in "([":
            depth += 1
        elif c in ")]":
            depth -= 1
        yield i, c, depth
        i += 1


def find_top(text: str, ch: str) -> int:
    for i, c, d in _scan(text):
        if c == ch and d == 0:
            return i
    return -1


def _find_keyword(text: str, word: str) -> int:
    """Posição de `word` como palavra solta no nível superior (ou -1)."""
    for i, c, d in _scan(text):
        if d != 0 or c != word[0]:
            continue
        if not text.startswith(word, i):
            continue
        before = text[i - 1] if i > 0 else " "
        after = text[i + len(word)] if i + len(word) < len(text) else " "
        if before in " \t" and after in " \t(":
            return i
        if i == 0 and after in " \t":
            return i
    return -1


def _find_weight(text: str) -> tuple[int, int]:
    for i, c, d in _scan(text):
        if c == "[" and d == 1:
            j = i
            depth = 0
            for k in range(i, len(text)):
                if text[k] == "[":
                    depth += 1
                elif text[k] == "]":
                    depth -= 1
                    if depth == 0:
                        j = k
                        break
            if _WEIGHT.match(text[i:j + 1]):
                before = text[:i].rstrip()
                # x[_n-1] é subscrito, não peso
                if before and (before[-1].isalnum() or before[-1] == "_") and not text[i - 1].isspace():
                    continue
                return i, j + 1
    return -1, -1


def _find_assign(text: str) -> int:
    for i, c, d in _scan(text):
        if c != "=" or d != 0:
            continue
        prev = text[i - 1] if i > 0 else ""
        nxt = text[i + 1] if i + 1 < len(text) else ""
        if prev in "=!<>~" or nxt == "=":
            continue
        return i
    return -1


def parse_standard(text: str) -> Parsed:
    p = Parsed()
    t = text.strip()
    comma = find_top(t, ",")
    if comma != -1:
        p.options = t[comma + 1:].strip()
        p.has_comma = True
        t = t[:comma].rstrip()

    # peso
    a, b = _find_weight(t)
    if a != -1:
        m = _WEIGHT.match(t[a:b])
        assert m
        p.weight = (_WEIGHT_NAMES[m.group(1)], m.group(2))
        t = (t[:a] + " " + t[b:]).strip()

    # using, if e in: localiza e recorta do fim para o começo
    marks = []
    for kw in ("using", "if", "in"):
        pos = _find_keyword(" " + t, kw)
        if pos != -1:
            marks.append((pos - 1, kw))
    marks.sort()
    for idx, (pos, kw) in enumerate(marks):
        end = marks[idx + 1][0] if idx + 1 < len(marks) else len(t)
        value = t[pos + len(kw):end].strip()
        if not value:
            raise StataError(198, "invalid syntax")
        if kw == "if":
            p.if_ = value
        elif kw == "in":
            p.in_ = value
        else:
            p.using = value
    head = t[:marks[0][0]].strip() if marks else t

    eq = _find_assign(head)
    if eq != -1:
        p.exp = head[eq + 1:].strip()
        if not p.exp:
            raise StataError(198, "invalid syntax")
        head = head[:eq].strip()
    p.varlist = head
    return p


# ---------------------------------------------------------------------------
# Opções
# ---------------------------------------------------------------------------

def parse_options(text: str) -> list[tuple[str, str | None]]:
    out: list[tuple[str, str | None]] = []
    i, n = 0, len(text)
    while i < n:
        while i < n and text[i] in " \t":
            i += 1
        if i >= n:
            break
        m = re.match(r"[A-Za-z_][A-Za-z0-9_]*", text[i:])
        if not m:
            raise StataError(198, "invalid syntax")
        name = m.group(0)
        i += m.end()
        arg = None
        if i < n and text[i] == "(":
            depth, j = 0, i
            while j < n:
                if text[j] == "(":
                    depth += 1
                elif text[j] == ")":
                    depth -= 1
                    if depth == 0:
                        break
                elif text[j] == '"':
                    k = text.find('"', j + 1)
                    j = n - 1 if k == -1 else k
                j += 1
            if depth != 0:
                raise StataError(198, "invalid syntax")
            arg = text[i + 1:j]
            i = j + 1
        out.append((name, arg))
    return out


def match_options(text: str, spec: dict[str, int]) -> dict[str, str | bool]:
    """`spec`: nome completo -> tamanho mínimo da abreviação.
    Devolve {nome completo: argumento ou True}."""
    result: dict[str, str | bool] = {}
    for name, arg in parse_options(text):
        for full, minimum in spec.items():
            if len(name) >= minimum and full.startswith(name):
                result[full] = arg if arg is not None else True
                break
        else:
            raise StataError(198, f"option {name} not allowed")
    return result


# ---------------------------------------------------------------------------
# in range
# ---------------------------------------------------------------------------

def parse_in(text: str, nobs: int) -> tuple[int, int]:
    """Converte '1/10', 'f/l', '-5/l', '5' em (início, fim) 0-based, fim exclusivo."""
    t = text.replace(" ", "")

    def one(tok: str) -> int:
        if tok == "f":
            return 1
        if tok == "l":
            return nobs
        try:
            k = int(float(tok))
        except ValueError:
            raise StataError(198, "invalid syntax")
        if k < 0:
            k = nobs + k + 1
        return k

    if "/" in t:
        a, b = t.split("/", 1)
        lo, hi = one(a), one(b)
    else:
        lo = hi = one(t)
    if lo < 1 or hi > nobs or lo > hi:
        if nobs == 0 and lo == 1:
            return 0, 0
        raise StataError(198, "Obs. nos. out of range")
    return lo - 1, hi


@dataclass
class ByPrefix:
    keys: list[str] = field(default_factory=list)
    sort_extra: list[str] = field(default_factory=list)
    sort: bool = False
    rc0: bool = False


def parse_by(text: str) -> ByPrefix:
    """Analisa o que vem entre 'by' e ':' — 'g1 g2 (t)', ', sort'."""
    t = text.strip()
    bp = ByPrefix()
    comma = find_top(t, ",")
    if comma != -1:
        opts = match_options(t[comma + 1:], {"sort": 1, "rc0": 3})
        bp.sort = bool(opts.get("sort"))
        bp.rc0 = bool(opts.get("rc0"))
        t = t[:comma]
    m = re.search(r"\(([^)]*)\)\s*$", t)
    if m:
        bp.sort_extra = m.group(1).split()
        t = t[:m.start()]
    bp.keys = t.split()
    if not bp.keys:
        raise StataError(100, "varlist required")
    return bp

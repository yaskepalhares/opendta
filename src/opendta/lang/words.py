"""Divisão de texto em palavras (com aspas) e listas numéricas (numlist)."""

from __future__ import annotations

import math
import re

from ..core.errors import StataError


def split_words(text: str, *, keep_quotes: bool = False) -> list[str]:
    """Divide por espaços respeitando "..." e `"..."'.

    Com keep_quotes=False as aspas externas de cada palavra são removidas
    (comportamento de foreach ... in e tokenize).
    """
    words: list[str] = []
    i, n = 0, len(text)
    while i < n:
        while i < n and text[i] in " \t":
            i += 1
        if i >= n:
            break
        if text.startswith('`"', i):
            depth, j = 1, i + 2
            while j < n and depth:
                if text.startswith('`"', j):
                    depth += 1
                    j += 2
                elif text.startswith("\"'", j):
                    depth -= 1
                    j += 2
                else:
                    j += 1
            word = text[i:j]
            words.append(word if keep_quotes else word[2:-2])
            i = j
            continue
        if text[i] == '"':
            j = text.find('"', i + 1)
            j = n if j == -1 else j + 1
            word = text[i:j]
            words.append(word if keep_quotes else word[1:-1] if word.endswith('"') and len(word) > 1 else word[1:])
            i = j
            continue
        j = i
        while j < n and text[j] not in " \t":
            if text[j] == '"':
                k = text.find('"', j + 1)
                j = n if k == -1 else k + 1
                continue
            j += 1
        words.append(text[i:j])
        i = j
    return words


def strip_outer_quotes(text: str) -> str:
    """Remove uma camada de aspas externas, como `local x "abc"` faz."""
    t = text.strip()
    if t.startswith('`"') and t.endswith("\"'") and len(t) >= 4:
        return t[2:-2]
    if len(t) >= 2 and t[0] == '"' and t[-1] == '"':
        return t[1:-1]
    return text.strip()


_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
_RANGE_SLASH = re.compile(rf"^({_NUM})/({_NUM})$")
_RANGE_STEP = re.compile(rf"^({_NUM})[\(\[]({_NUM})[\)\]]({_NUM})$")


def _frange(a: float, step: float, b: float) -> list[float]:
    if step == 0:
        raise StataError(121, "invalid numlist")
    out: list[float] = []
    k = 0
    eps = abs(step) * 1e-9
    while True:
        v = a + k * step
        if (step > 0 and v > b + eps) or (step < 0 and v < b - eps):
            break
        out.append(round(v, 12))
        k += 1
        if k > 10_000_000:
            raise StataError(121, "invalid numlist has too many elements")
    return out


def parse_numlist(text: str) -> list[float]:
    """Expande uma numlist: '1/5', '1(2)9', '10(-2)0', '1 2 3', '1[2]9'."""
    # junta "1 (2) 9" e "1 / 5" em tokens únicos
    t = re.sub(r"\s*/\s*", "/", text.strip())
    t = re.sub(r"\s*([\(\[])\s*", r"\1", t)
    t = re.sub(r"\s*([\)\]])\s*", r"\1", t)
    out: list[float] = []
    for tok in t.split():
        m = _RANGE_SLASH.match(tok)
        if m:
            a, b = float(m.group(1)), float(m.group(2))
            step = 1.0 if b >= a else -1.0
            out.extend(_frange(a, step, b))
            continue
        m = _RANGE_STEP.match(tok)
        if m:
            out.extend(_frange(float(m.group(1)), float(m.group(2)), float(m.group(3))))
            continue
        if tok == ".":
            from ..core.missing import SYSMISS
            out.append(SYSMISS)
            continue
        try:
            out.append(float(tok))
        except ValueError:
            raise StataError(121, "invalid numlist")
    return out


def num_str(x: float) -> str:
    """Número como aparece numa macro gerada por forvalues/numlist."""
    if math.isfinite(x) and x == int(x):
        return str(int(x))
    s = repr(x)
    return s[1:] if s.startswith("0.") else ("-" + s[2:] if s.startswith("-0.") else s)

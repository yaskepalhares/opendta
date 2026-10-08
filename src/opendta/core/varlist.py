"""Expansão de varlists (manual [U] 11.4).

    mpg            nome exato
    mp             abreviação (se `set varabbrev on`, padrão)
    m*  *g  m?g    curingas: * (zero ou mais), ? (exatamente um)
    m~             como *, mas deve corresponder a uma única variável
    make-weight    intervalo, na ordem do conjunto de dados
    _all  ou  *    todas as variáveis
"""

from __future__ import annotations

import re

from .dataset import Dataset
from .errors import StataError

_WILD = re.compile(r"[*?~]")


def _pattern(token: str) -> re.Pattern[str]:
    out = []
    for ch in token:
        if ch in "*~":
            out.append(".*")
        elif ch == "?":
            out.append(".")
        else:
            out.append(re.escape(ch))
    return re.compile("^" + "".join(out) + "$")


def resolve_name(ds: Dataset, token: str, *, abbrev: bool = True) -> str:
    names = ds.names
    if token in names:
        return token
    if abbrev:
        matches = [n for n in names if n.startswith(token)]
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise StataError(111, f"{token} ambiguous abbreviation")
    raise StataError(111, f"variable {token} not found")


def expand(ds: Dataset, text: str, *, abbrev: bool = True, allow_empty: bool = True) -> list[str]:
    t = re.sub(r"\s*-\s*", "-", text.strip())
    if not t:
        if allow_empty:
            return []
        raise StataError(100, "varlist required")
    out: list[str] = []
    for token in t.split():
        if token in ("_all", "*"):
            out.extend(ds.names)
            continue
        if "-" in token and not token.startswith("-"):
            a, b = token.split("-", 1)
            ia = ds.index(resolve_name(ds, a, abbrev=abbrev))
            ib = ds.index(resolve_name(ds, b, abbrev=abbrev))
            if ib < ia:
                raise StataError(198, f"{token}: variables out of order")
            out.extend(ds.names[ia:ib + 1])
            continue
        if _WILD.search(token):
            pat = _pattern(token)
            matches = [n for n in ds.names if pat.match(n)]
            if not matches:
                raise StataError(111, f"variable {token} not found")
            if "~" in token and len(matches) > 1:
                raise StataError(111, f"{token} ambiguous abbreviation")
            out.extend(matches)
            continue
        out.append(resolve_name(ds, token, abbrev=abbrev))
    return out


def unique(names: list[str]) -> list[str]:
    seen: list[str] = []
    for n in names:
        if n not in seen:
            seen.append(n)
    return seen

"""What every explanation rule shares: the explanation itself, the context of the failed
command and the table of rules."""

from __future__ import annotations

import difflib
import re
from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable

import numpy as np

from ...i18n import L

if TYPE_CHECKING:
    from ...session import Session
    from ..errors import StataError


@dataclass
class Explanation:
    """what: what went wrong in this command, with the names and numbers of the data at hand;
    why: the idea behind the rule, for someone new to Stata; fix: what to type or do instead;
    help: a help topic for the "learn more" link (the command's own page when empty)."""
    what: str
    why: str = ""
    fix: str = ""
    help: str = ""


@dataclass
class Ctx:
    """The failed command: its error and its text, with helpers that look at the data."""
    err: "StataError"
    s: "Session"
    command: str          # command name ("" when unknown)
    text: str             # the command as typed

    @property
    def msg(self) -> str:
        return self.err.message or ""

    @property
    def rc(self) -> int:
        return self.err.rc

    @property
    def ds(self):
        return self.s.data

    @property
    def names(self) -> list[str]:
        return self.s.data.names

    @property
    def args(self) -> str:
        """The text after the command word (and after any prefix such as by ...:)."""
        body = self.text.split(":", 1)[1] if re.match(r"\s*(by|bys|bysort)\b[^:]*:", self.text) else self.text
        parts = body.strip().split(None, 1)
        return parts[1] if len(parts) > 1 else ""

    def parsed(self):
        from ...lang.syntax import parse_standard
        try:
            return parse_standard(self.args)
        except Exception:  # noqa: BLE001 - a command the standard syntax cannot read
            return None

    def variables(self, text: str | None = None) -> list[str]:
        """Variables named in the text (in order, without repetition; abbreviations resolved)."""
        from ...core.varlist import resolve_name
        out: list[str] = []
        src = self.args if text is None else text
        for tok in re.findall(r'"[^"]*"|\b[A-Za-z_][A-Za-z0-9_]*\b', src):
            if tok.startswith('"'):
                continue
            try:
                name = resolve_name(self.ds, tok)
            except Exception:  # noqa: BLE001 - not a variable
                continue
            if name not in out:
                out.append(name)
        return out

    def sample(self) -> np.ndarray:
        """Observations selected by the command's if and in (all of them if neither)."""
        from ...commands._util import touse
        p = self.parsed()
        try:
            return touse(self.s, p) if p is not None else np.ones(self.ds.nobs, bool)
        except Exception:  # noqa: BLE001 - an if condition that cannot be evaluated
            return np.ones(self.ds.nobs, bool)

    def is_string(self, name: str) -> bool:
        return self.ds.has(name) and self.ds.get(name).is_string


Rule = Callable[[Ctx, "re.Match | None"], "Explanation | None"]
RULES: list[tuple[re.Pattern | None, int | None, Rule]] = []


def rule(pattern: str | None = None, *, rc: int | None = None):
    """Registers an explanation for the messages that fully match `pattern` (and, when given,
    the return code rc). The function returns None to let the next rule try."""
    rx = re.compile(pattern, re.S) if pattern is not None else None

    def deco(fn: Rule) -> Rule:
        RULES.append((rx, rc, fn))
        return fn
    return deco


# ---------------------------------------------------------------------------
# small text helpers
# ---------------------------------------------------------------------------

def close(word: str, candidates, n: int = 3) -> list[str]:
    return difflib.get_close_matches(word, list(candidates), n=n, cutoff=0.6)


def did_you_mean(word: str, candidates) -> str:
    near = close(word, candidates)
    if not near:
        return ""
    return L("did you mean ", "você quis dizer ") + ", ".join(f'"{x}"' for x in near) + "?"


def join(items: list[str], limit: int = 8) -> str:
    shown = items[:limit]
    more = (L(f" (and {len(items) - limit} more)", f" (e mais {len(items) - limit})")
            if len(items) > limit else "")
    return ", ".join(shown) + more


def and_join(items: list[str]) -> str:
    """'a, b and c' ('a, b e c')."""
    if len(items) < 2:
        return "".join(items)
    return ", ".join(items[:-1]) + L(" and ", " e ") + items[-1]


def count(n: int, en: str, pt: str, en_pl: str | None = None, pt_pl: str | None = None) -> str:
    """'1 observation', '3 observations' (with the Portuguese forms)."""
    one = n == 1
    word = L(en if one else (en_pl or en + "s"), pt if one else (pt_pl or pt + "s"))
    return f"{n:,} {word}"


def obs(n: int) -> str:
    return count(n, "observation", "observação", "observations", "observações")


def values_summary(x: np.ndarray, limit: int = 6) -> str:
    """The distinct values of a variable with their counts: '0 (12), 1 (30)'."""
    vals, cnt = np.unique(x, return_counts=True)
    parts = [f"{fmt_value(v)} ({c:,})" for v, c in zip(vals[:limit], cnt[:limit])]
    if len(vals) > limit:
        parts.append(L(f"and {len(vals) - limit} more values", f"e mais {len(vals) - limit} valores"))
    return ", ".join(parts)


def fmt_value(v) -> str:
    if isinstance(v, str):
        return f'"{v}"'
    v = float(v)
    from ..missing import SYSMISS
    if v >= SYSMISS:
        return "."
    return f"{int(v)}" if v == int(v) and abs(v) < 1e15 else f"{v:g}"

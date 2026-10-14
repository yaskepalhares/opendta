"""Command registry and abbreviation resolution.

Each command declares its full name and minimal abbreviation (the
underlined part in the manual's syntax diagram). E.g.: display -> "di".
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable

if TYPE_CHECKING:
    from ..session import Session

CommandFn = Callable[["Session", str], None]


@dataclass(frozen=True)
class CommandSpec:
    name: str
    min_abbrev: int
    fn: CommandFn
    prefix: bool = False  # quietly, noisily, capture, version: take a command
    byable: bool = False  # aceita o prefixo by


REGISTRY: dict[str, CommandSpec] = {}


def command(name: str, abbrev: str | None = None, *, prefix: bool = False, byable: bool = False):
    """Registers a command. `abbrev` is the minimal accepted form (default: full name)."""
    minimal = abbrev or name
    if not name.startswith(minimal):
        raise ValueError(f"abbreviation {minimal!r} is not a prefix of {name!r}")

    def deco(fn: CommandFn) -> CommandFn:
        REGISTRY[name] = CommandSpec(name, len(minimal), fn, prefix, byable)
        return fn
    return deco


def lookup(word: str) -> CommandSpec | None:
    spec = REGISTRY.get(word)
    if spec is not None:
        return spec
    # with more than one candidate (sca: scalar or scatter), the one with the longest minimal
    # abbreviation wins: it is the one the manual reserves for that form
    best = None
    for spec in REGISTRY.values():
        if len(word) >= spec.min_abbrev and spec.name.startswith(word):
            if best is None or spec.min_abbrev > best.min_abbrev:
                best = spec
    return best


def all_commands() -> list[str]:
    return sorted(REGISTRY)

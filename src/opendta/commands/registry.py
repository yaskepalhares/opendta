"""Registro de comandos e resolução de abreviações.

Cada comando declara o nome completo e a abreviação mínima (a parte
sublinhada no diagrama de sintaxe do manual). Ex.: display -> "di".
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
    prefix: bool = False  # quietly, noisily, capture, version: recebem um comando
    byable: bool = False  # aceita o prefixo by


REGISTRY: dict[str, CommandSpec] = {}


def command(name: str, abbrev: str | None = None, *, prefix: bool = False, byable: bool = False):
    """Registra um comando. `abbrev` é a forma mínima aceita (padrão: nome completo)."""
    minimal = abbrev or name
    if not name.startswith(minimal):
        raise ValueError(f"abreviação {minimal!r} não é prefixo de {name!r}")

    def deco(fn: CommandFn) -> CommandFn:
        REGISTRY[name] = CommandSpec(name, len(minimal), fn, prefix, byable)
        return fn
    return deco


def lookup(word: str) -> CommandSpec | None:
    spec = REGISTRY.get(word)
    if spec is not None:
        return spec
    for spec in REGISTRY.values():
        if len(word) >= spec.min_abbrev and spec.name.startswith(word):
            return spec
    return None


def all_commands() -> list[str]:
    return sorted(REGISTRY)

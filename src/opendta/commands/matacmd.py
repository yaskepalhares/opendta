"""Comando mata numa linha só: `mata: instrução`, `mata instrução` e os
subcomandos `mata clear`, `mata describe`, `mata drop`. O bloco
`mata ... end` é tratado pelo interpretador (lang/interpreter.py)."""

from __future__ import annotations

from typing import TYPE_CHECKING

from .registry import command

if TYPE_CHECKING:
    from ..session import Session

_SUBCOMMANDS = {"clear", "describe", "d", "drop", "set", "mlib", "mosave", "which", "memory",
                "query", "rename"}


@command("mata")
def cmd_mata(s: "Session", args: str) -> None:
    text = args.strip()
    if text.startswith(":"):
        text = text[1:].strip()
    if not text:
        return
    first = text.split()[0]
    from ..mata.interp import MataLeave
    try:
        if first in _SUBCOMMANDS:
            s.mata.subcommand(text)
        else:
            s.mata.run_statement(text)
    except MataLeave:
        pass

"""Comandos que abrem janelas da interface: doedit."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ..lang.words import strip_outer_quotes
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


@command("doedit", "doed")
def cmd_doedit(s: "Session", args: str) -> None:
    """Abre o Do-file Editor (sem interface, não faz nada, como o Stata em
    modo console)."""
    hook = s.ui_hooks.get("doedit")
    if hook is not None:
        hook(strip_outer_quotes(args.strip()))

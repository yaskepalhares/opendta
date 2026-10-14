"""Commands that open interface windows: doedit."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ..lang.words import strip_outer_quotes
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


@command("doedit", "doed")
def cmd_doedit(s: "Session", args: str) -> None:
    """Opens the Do-file Editor (without the interface it does nothing, as Stata in
    console mode)."""
    hook = s.ui_hooks.get("doedit")
    if hook is not None:
        hook(strip_outer_quotes(args.strip()))

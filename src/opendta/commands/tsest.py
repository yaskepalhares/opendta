"""Time-series estimation and test commands of phase 18, loaded on first use (start-up time)."""

from __future__ import annotations

from typing import TYPE_CHECKING

from .registry import command

if TYPE_CHECKING:
    from ..session import Session


@command("dfgls", "dfgls")
def cmd_dfgls(s: "Session", args: str) -> None:
    from ..estimation.tsunit import cmd_dfgls as run
    run(s, args)


@command("pperron", "pperron")
def cmd_pperron(s: "Session", args: str) -> None:
    from ..estimation.tsunit import cmd_pperron as run
    run(s, args)


@command("prais", "prais")
def cmd_prais(s: "Session", args: str) -> None:
    from ..estimation.prais import cmd_prais as run
    run(s, args)

"""Survival-analysis commands of the survival block (phase 15): the logic lives in
opendta/commands/stdata.py, strates.py, ctdata.py and opendta/estimation/stcrreg.py; here only the
names are registered."""

from __future__ import annotations

from typing import TYPE_CHECKING

from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def _lazy(module: str, name: str):
    def run(s: "Session", args: str) -> None:
        import importlib
        pkg = "opendta.estimation" if module in ("stcrreg",) else "opendta.commands"
        getattr(importlib.import_module(f"{pkg}.{module}"), name)(s, args)
    return run


for _cmd, _mod, _fn in (
    ("stsplit", "stdata", "cmd_stsplit"), ("stjoin", "stdata", "cmd_stjoin"),
    ("stfill", "stdata", "cmd_stfill"), ("stgen", "stdata", "cmd_stgen"),
    ("stvary", "stdata", "cmd_stvary"), ("stbase", "stdata", "cmd_stbase"),
    ("snapspan", "stdata", "cmd_snapspan"),
):
    command(_cmd, _cmd)(_lazy(_mod, _fn))

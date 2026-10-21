"""Registration of the model-fitting commands of the Models block (scobit, hetprobit, fracreg,
zip, zinb, tpoisson, tnbreg, gnbreg, binreg...). The logic lives in opendta/estimation; here only
the names and abbreviations are registered."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING

from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def _register(module: str, names: list[str]) -> None:
    for name in names:
        def run(s: "Session", args: str, _name=name, _module=module) -> None:
            getattr(import_module(f"opendta.estimation.{_module}"), f"cmd_{_name}")(s, args)
        run.__name__ = f"cmd_{name}"
        command(name, name)(run)


_register("binmodels", ["scobit", "hetprobit", "fracreg", "zip", "zinb", "tpoisson", "tnbreg", "gnbreg", "binreg", "glogit"])

_register("censmodels", ["intreg", "truncreg", "heckman", "heckprobit", "heckoprobit", "biprobit"])

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


def _varpost(name):
    def run(s: "Session", args: str) -> None:
        from ..estimation import varpost
        getattr(varpost, f"cmd_{name}")(s, args)
    run.__name__ = f"cmd_{name}"
    return run


for _n in ("varlmar", "varnorm", "varstable", "varwle"):
    command(_n, _n)(_varpost(_n))


command("fcast", "fcast")(_varpost("fcast"))

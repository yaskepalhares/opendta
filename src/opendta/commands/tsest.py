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


for _n in ("varlmar", "varnorm", "varstable", "varwle", "varbasic"):
    command(_n, _n)(_varpost(_n))


command("fcast", "fcast")(_varpost("fcast"))


def _vec(name):
    def run(s: "Session", args: str) -> None:
        from ..estimation import vec
        getattr(vec, f"cmd_{name}")(s, args)
    run.__name__ = f"cmd_{name}"
    return run


for _n in ("vec", "vecrank"):
    command(_n, _n)(_vec(_n))


@command("svar", "svar")
def cmd_svar(s: "Session", args: str) -> None:
    from ..estimation.svar import cmd_svar as run
    run(s, args)


@command("arch", "arch")
def cmd_arch(s: "Session", args: str) -> None:
    from ..estimation.archmodel import cmd_arch as run
    run(s, args)


@command("arfima", "arfima")
def cmd_arfima(s: "Session", args: str) -> None:
    from ..estimation.arfima import cmd_arfima as run
    run(s, args)


@command("ucm", "ucm")
def cmd_ucm(s: "Session", args: str) -> None:
    from ..estimation.ucm import cmd_ucm as run
    run(s, args)


@command("sspace", "sspace")
def cmd_sspace(s: "Session", args: str) -> None:
    from ..estimation.sspace import cmd_sspace as run
    run(s, args)


@command("dfactor", "dfactor")
def cmd_dfactor(s: "Session", args: str) -> None:
    from ..estimation.sspace import cmd_dfactor as run
    run(s, args)


@command("mgarch", "mgarch")
def cmd_mgarch(s: "Session", args: str) -> None:
    from ..estimation.mgarch import cmd_mgarch as run
    run(s, args)

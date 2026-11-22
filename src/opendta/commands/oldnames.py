"""Names from earlier Stata versions that Stata 14 still runs (compat 2005): hettest, ovtest,
lstat, vif, treatreg, svydes and svylc."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ..core.errors import StataError
from .registry import REGISTRY, command

if TYPE_CHECKING:
    from ..session import Session


def _run(s: "Session", name: str, args: str) -> None:
    REGISTRY[name].fn(s, args)


@command("hettest")
def cmd_hettest(s: "Session", args: str) -> None:
    _run(s, "estat", "hettest " + args)


@command("ovtest")
def cmd_ovtest(s: "Session", args: str) -> None:
    _run(s, "estat", "ovtest " + args)


@command("lstat")
def cmd_lstat(s: "Session", args: str) -> None:
    _run(s, "estat", "classification " + args)


@command("vif")
def cmd_vif(s: "Session", args: str) -> None:
    # without estimates Stata's vif reports the missing constant (compat 2005)
    if not s.e:
        raise StataError(301, "not appropriate after regress, nocons;\n"
                              "use option uncentered to get uncentered VIFs")
    _run(s, "estat", "vif " + args)


@command("treatreg")
def cmd_treatreg(s: "Session", args: str) -> None:
    _run(s, "etregress", args)


@command("svydes")
def cmd_svydes(s: "Session", args: str) -> None:
    _run(s, "svydescribe", args)


@command("svylc")
def cmd_svylc(s: "Session", args: str) -> None:
    """svylc exp: lincom after svy. VERIFY: Stata 14 runs it as lincom."""
    if not args.strip() or args.strip().startswith(","):
        raise StataError(198, "invalid syntax")   # compat 2005
    _run(s, "lincom", args)

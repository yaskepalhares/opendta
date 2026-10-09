"""Comandos de estimação e pós-estimação (fase 5). A lógica fica em
opendta/estimation; aqui só o registro dos nomes e abreviações."""

from __future__ import annotations

from typing import TYPE_CHECKING

from .registry import command

if TYPE_CHECKING:
    from ..session import Session


@command("tsset", "tsset")
def cmd_tsset(s: "Session", args: str) -> None:
    from ..estimation.tsops import cmd_tsset as run
    run(s, args)


@command("xtset", "xtset")
def cmd_xtset(s: "Session", args: str) -> None:
    from ..estimation.tsops import cmd_xtset as run
    run(s, args)


@command("regress", "reg")
def cmd_regress(s: "Session", args: str) -> None:
    from ..estimation.regress import cmd_regress as run
    run(s, args)


@command("predict", "predict")
def cmd_predict(s: "Session", args: str) -> None:
    from ..estimation.postest import cmd_predict as run
    run(s, args)


@command("test", "test")
def cmd_test(s: "Session", args: str) -> None:
    from ..estimation.postest import cmd_test as run
    run(s, args)


@command("testparm", "testparm")
def cmd_testparm(s: "Session", args: str) -> None:
    from ..estimation.postest import cmd_testparm as run
    run(s, args)


@command("lincom", "lincom")
def cmd_lincom(s: "Session", args: str) -> None:
    from ..estimation.postest import cmd_lincom as run
    run(s, args)


@command("estimates", "est")
def cmd_estimates(s: "Session", args: str) -> None:
    from ..estimation.postest import cmd_estimates as run
    run(s, args)


@command("logit", "logit")
def cmd_logit(s: "Session", args: str) -> None:
    from ..estimation.models import cmd_logit as run
    run(s, args)


@command("logistic", "logistic")
def cmd_logistic(s: "Session", args: str) -> None:
    from ..estimation.models import cmd_logistic as run
    run(s, args)


@command("probit", "probit")
def cmd_probit(s: "Session", args: str) -> None:
    from ..estimation.models import cmd_probit as run
    run(s, args)


@command("poisson", "poisson")
def cmd_poisson(s: "Session", args: str) -> None:
    from ..estimation.models import cmd_poisson as run
    run(s, args)

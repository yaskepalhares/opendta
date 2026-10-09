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


@command("ologit", "ologit")
def cmd_ologit(s: "Session", args: str) -> None:
    from ..estimation.models2 import cmd_ologit as run
    run(s, args)


@command("oprobit", "oprobit")
def cmd_oprobit(s: "Session", args: str) -> None:
    from ..estimation.models2 import cmd_oprobit as run
    run(s, args)


@command("mlogit", "mlogit")
def cmd_mlogit(s: "Session", args: str) -> None:
    from ..estimation.models2 import cmd_mlogit as run
    run(s, args)


@command("nbreg", "nbreg")
def cmd_nbreg(s: "Session", args: str) -> None:
    from ..estimation.models2 import cmd_nbreg as run
    run(s, args)


@command("tobit", "tobit")
def cmd_tobit(s: "Session", args: str) -> None:
    from ..estimation.models2 import cmd_tobit as run
    run(s, args)


@command("glm", "glm")
def cmd_glm(s: "Session", args: str) -> None:
    from ..estimation.glm import cmd_glm as run
    run(s, args)


@command("xtreg", "xtreg")
def cmd_xtreg(s: "Session", args: str) -> None:
    from ..estimation.panel import cmd_xtreg as run
    run(s, args)


@command("areg", "areg")
def cmd_areg(s: "Session", args: str) -> None:
    from ..estimation.panel import cmd_areg as run
    run(s, args)


@command("ivregress", "ivregress")
def cmd_ivregress(s: "Session", args: str) -> None:
    from ..estimation.panel import cmd_ivregress as run
    run(s, args)


@command("margins", "margins")
def cmd_margins(s: "Session", args: str) -> None:
    from ..estimation.margins import cmd_margins as run
    run(s, args)


@command("nlcom", "nlcom")
def cmd_nlcom(s: "Session", args: str) -> None:
    from ..estimation.margins import cmd_nlcom as run
    run(s, args)


@command("svyset", "svyset")
def cmd_svyset(s: "Session", args: str) -> None:
    from ..estimation.svy import cmd_svyset as run
    run(s, args)


@command("svydescribe", "svydescribe")
def cmd_svydescribe(s: "Session", args: str) -> None:
    from ..estimation.svy import cmd_svydescribe as run
    run(s, args)


@command("mean", "mean")
def cmd_mean(s: "Session", args: str) -> None:
    from ..estimation.svy import cmd_mean as run
    run(s, args)


@command("proportion", "proportion")
def cmd_proportion(s: "Session", args: str) -> None:
    from ..estimation.svy import cmd_proportion as run
    run(s, args)


@command("total", "total")
def cmd_total(s: "Session", args: str) -> None:
    from ..estimation.svy import cmd_total as run
    run(s, args)


@command("ratio", "ratio")
def cmd_ratio(s: "Session", args: str) -> None:
    from ..estimation.svy import cmd_ratio as run
    run(s, args)


@command("svy", "svy", prefix=True)
def cmd_svy(s: "Session", args: str) -> None:
    from ..estimation.svy import cmd_svy as run
    run(s, args)

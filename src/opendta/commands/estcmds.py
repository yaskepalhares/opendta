"""Estimation and postestimation commands (phase 5). The logic lives in
opendta/estimation; here only the names and abbreviations are registered."""

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


@command("estat", "estat")
def cmd_estat(s: "Session", args: str) -> None:
    from ..estimation.estat import cmd_estat as run
    run(s, args)


@command("oneway", "onew")
def cmd_oneway(s: "Session", args: str) -> None:
    from ..estimation.anova import cmd_oneway as run
    run(s, args)


@command("anova", "anova")
def cmd_anova(s: "Session", args: str) -> None:
    from ..estimation.anova import cmd_anova as run
    run(s, args)


@command("bootstrap", "bootstrap", prefix=True)
def cmd_bootstrap(s: "Session", args: str) -> None:
    from ..estimation.bootstrap import cmd_bootstrap as run
    run(s, args)


@command("lrtest", "lrtest")
def cmd_lrtest(s: "Session", args: str) -> None:
    from ..estimation.extra import cmd_lrtest as run
    run(s, args)


@command("hausman", "hausman")
def cmd_hausman(s: "Session", args: str) -> None:
    from ..estimation.extra import cmd_hausman as run
    run(s, args)


@command("xtsum", "xtsum")
def cmd_xtsum(s: "Session", args: str) -> None:
    from ..estimation.extra import cmd_xtsum as run
    run(s, args)


@command("xtdescribe", "xtdes")
def cmd_xtdescribe(s: "Session", args: str) -> None:
    from ..estimation.extra import cmd_xtdescribe as run
    run(s, args)


@command("qreg", "qreg")
def cmd_qreg(s: "Session", args: str) -> None:
    from ..estimation.extra import cmd_qreg as run
    run(s, args)


@command("xtpoisson", "xtpoisson")
def cmd_xtpoisson(s: "Session", args: str) -> None:
    from ..estimation.panel2 import cmd_xtpoisson as run
    run(s, args)


@command("xtlogit", "xtlogit")
def cmd_xtlogit(s: "Session", args: str) -> None:
    from ..estimation.panel2 import cmd_xtlogit as run
    run(s, args)


@command("stset", "stset")
def cmd_stset(s: "Session", args: str) -> None:
    from ..estimation.st import cmd_stset as run
    run(s, args)


@command("st", "st")
def cmd_st(s: "Session", args: str) -> None:
    from ..estimation.st import cmd_st as run
    run(s, args)


@command("streset", "streset")
def cmd_streset(s: "Session", args: str) -> None:
    from ..estimation.st import cmd_streset as run
    run(s, args)


@command("stdescribe", "stdes")
def cmd_stdescribe(s: "Session", args: str) -> None:
    from ..estimation.st import cmd_stdescribe as run
    run(s, args)


@command("stsum", "stsum")
def cmd_stsum(s: "Session", args: str) -> None:
    from ..estimation.st import cmd_stsum as run
    run(s, args)


@command("sts", "sts")
def cmd_sts(s: "Session", args: str) -> None:
    from ..estimation.km import cmd_sts as run
    run(s, args)


@command("stcox", "stcox")
def cmd_stcox(s: "Session", args: str) -> None:
    from ..estimation.stmodels import cmd_stcox as run
    run(s, args)


@command("streg", "streg")
def cmd_streg(s: "Session", args: str) -> None:
    from ..estimation.stmodels import cmd_streg as run
    run(s, args)


@command("mixed", "mixed")
def cmd_mixed(s: "Session", args: str) -> None:
    from ..estimation.mixed import cmd_mixed as run
    run(s, args)


@command("xtmixed", "xtmixed")
def cmd_xtmixed(s: "Session", args: str) -> None:
    from ..estimation.mixed import cmd_xtmixed as run
    run(s, args)


@command("melogit", "melogit")
def cmd_melogit(s: "Session", args: str) -> None:
    from ..estimation.meglm import cmd_melogit as run
    run(s, args)


@command("meprobit", "meprobit")
def cmd_meprobit(s: "Session", args: str) -> None:
    from ..estimation.meglm import cmd_meprobit as run
    run(s, args)


@command("mepoisson", "mepoisson")
def cmd_mepoisson(s: "Session", args: str) -> None:
    from ..estimation.meglm import cmd_mepoisson as run
    run(s, args)


@command("arima", "arima")
def cmd_arima(s: "Session", args: str) -> None:
    from ..estimation.tsmodels import cmd_arima as run
    run(s, args)


@command("var", "var")
def cmd_var(s: "Session", args: str) -> None:
    from ..estimation.tsmodels import cmd_var as run
    run(s, args)


@command("varsoc", "varsoc")
def cmd_varsoc(s: "Session", args: str) -> None:
    from ..estimation.tsmodels import cmd_varsoc as run
    run(s, args)


@command("vargranger", "vargranger")
def cmd_vargranger(s: "Session", args: str) -> None:
    from ..estimation.tsmodels import cmd_vargranger as run
    run(s, args)


@command("dfuller", "dfuller")
def cmd_dfuller(s: "Session", args: str) -> None:
    from ..estimation.tsmodels import cmd_dfuller as run
    run(s, args)


@command("corrgram", "corrgram")
def cmd_corrgram(s: "Session", args: str) -> None:
    from ..estimation.tsmodels import cmd_corrgram as run
    run(s, args)


@command("newey", "newey")
def cmd_newey(s: "Session", args: str) -> None:
    from ..estimation.tsmodels import cmd_newey as run
    run(s, args)


@command("mi", "mi")
def cmd_mi(s: "Session", args: str) -> None:
    from ..estimation.mi import cmd_mi as run
    run(s, args)


def _install_vce_bootstrap() -> None:
    from ..estimation.bootstrap import wrap_vce_bootstrap
    from .registry import REGISTRY
    for name in ("regress", "logit", "logistic", "probit", "poisson", "nbreg", "ologit", "oprobit",
                 "mlogit", "tobit", "glm", "areg"):
        if name in REGISTRY:
            wrap_vce_bootstrap(name)


_install_vce_bootstrap()

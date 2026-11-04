"""Estimation and postestimation commands (phase 5). The logic lives in
opendta/estimation; here only the names and abbreviations are registered."""

from __future__ import annotations

from .registry import lazy_commands

lazy_commands("tsops", "tsset", "xtset")
lazy_commands("regress", ("regress", "reg"))
lazy_commands("postest", "predict", "test", "testparm", "lincom", ("estimates", "est"))
lazy_commands("models", "logit", "logistic", "probit", "poisson", "cloglog")
lazy_commands("models2", "ologit", "oprobit", "mlogit", "nbreg", "tobit")
lazy_commands("glm", "glm")
lazy_commands("panel", "xtreg")
lazy_commands("xtgee", "xtgee")
lazy_commands("sem", "sem")
lazy_commands("gsem", "gsem")
lazy_commands("irt", "irt", "irtgraph")
lazy_commands("bayes", "bayesmh")
lazy_commands("bayes", "bayes", prefix=True)
lazy_commands("bayespost", "bayesstats", "bayestest", "bayesgraph")
lazy_commands("irf", "irf")
lazy_commands("panel", "areg", "ivregress")
lazy_commands("margins", "margins", "nlcom")
lazy_commands("svy", "svyset", "svydescribe")
lazy_commands("summaries", "mean", "proportion", "total", "ratio")
lazy_commands("svy", "svy", prefix=True)
lazy_commands("estat", "estat")
lazy_commands("anova", ("oneway", "onew"), "anova")
lazy_commands("bootstrap", "bootstrap", prefix=True)
lazy_commands("resample", ("jackknife", "jack"), "permute", ("simulate", "sim"), prefix=True)
lazy_commands("mlcmd", "ml")
lazy_commands("mlexp", "mlexp")
lazy_commands("mltools", "mleval", "mlsum", "mlvecsum", "mlmatsum")
lazy_commands("fpoly", "fp", "mfp", prefix=True)
lazy_commands("postest2", "testnl", "predictnl", "linktest", "contrast", "pwcompare", "pwmean", "suest")
lazy_commands("extra", "lrtest", "hausman", "xtsum", ("xtdescribe", "xtdes"), "qreg")
lazy_commands("panel2", "xtpoisson", "xtlogit")
lazy_commands("xtre", ("xtprobit", "xtprob"), "xtcloglog", "xtologit", "xtoprobit", "xttobit", "xtintreg",
                      "xtnbreg")
lazy_commands("xtlin", "xtgls", "xtpcse", "xtregar", "xtrc", "xthtaylor", "xtivreg")
lazy_commands("xtfrontier", "xtfrontier")
lazy_commands("xtroot", "xtunitroot")
lazy_commands("xtdpd", "xtabond", "xtdpdsys", "xtdpd")
lazy_commands("xtdesc", "xtline", "xttab", "xtdata")
lazy_commands("xtre", "quadchk")
lazy_commands("st", "stset", "st", "streset", ("stdescribe", "stdes"), "stsum")
lazy_commands("km", "sts")
lazy_commands("stcox", "stcox")
lazy_commands("streg", "streg")
lazy_commands("mixed", "mixed", "xtmixed")
lazy_commands("meglm", "melogit", "meprobit", "mepoisson", "mecloglog", "menbreg", "meologit", "meoprobit",
                       "meglm", "meqrlogit", "meqrpoisson")
lazy_commands("mestreg", "mestreg")
lazy_commands("tsmodels", "arima", "var", "varsoc", "vargranger", "dfuller", "corrgram", "newey")
lazy_commands("mi", "mi")


def _install_vce_bootstrap() -> None:
    from ..estimation.bootstrap import wrap_vce_bootstrap
    from .registry import REGISTRY
    for name in ("regress", "logit", "logistic", "probit", "cloglog", "poisson", "nbreg", "ologit", "oprobit",
                 "mlogit", "tobit", "glm", "areg"):
        if name in REGISTRY:
            wrap_vce_bootstrap(name)


_install_vce_bootstrap()

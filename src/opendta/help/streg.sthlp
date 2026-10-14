{smcl}
{title:Title}

{p 4 4 2}
{bf:streg} {hline 2} Parametric survival models

{title:Syntax}

{p 8 16 2}
{cmd:streg} [indepvars] [if] [in], {cmd:distribution(}{cmd:exponential}|{cmd:weibull}|{cmd:gompertz}|{cmd:lognormal}|{cmd:loglogistic}{cmd:)} [{cmd:time} {cmd:nohr} {cmd:tr} {cmd:vce(robust}|{cmd:cluster} var{cmd:)} {cmd:level(}#{cmd:)} {cmd:nolog} {cmd:noshow}]

{title:Description}

{pstd}
Fits parametric survival models by maximum likelihood, with censoring and delayed entry. The exponential and Weibull models are shown in the proportional-hazards metric (hazard ratios) unless {cmd:time} asks for the accelerated failure-time metric; Gompertz is proportional-hazards only; lognormal and loglogistic are accelerated failure-time models. The ancillary parameter ({cmd:/ln_p}, {cmd:/gamma}, {cmd:/ln_sig}, {cmd:/ln_gam}) and its transformations are shown below the coefficients. The generalized gamma distribution is not available.

{title:Postestimation}

{phang2}{cmd:predict} newvar [, {cmd:median time} {cmd:xb} {cmd:hr} {cmd:hazard} {cmd:surv} {cmd:csurv} {cmd:csnell} {cmd:mgale}]{p_end}

{smcl}
{title:Title}

{p 4 4 2}
{bf:stcox} {hline 2} Cox proportional hazards model

{title:Syntax}

{p 8 16 2}
{cmd:stcox} indepvars [if] [in] [, {cmd:nohr} {cmd:efron} {cmd:breslow} {cmd:strata(}varlist{cmd:)} {cmd:vce(robust}|{cmd:cluster} var{cmd:)} {cmd:basesurv(}newvar{cmd:)} {cmd:basechazard(}newvar{cmd:)} {cmd:basehc(}newvar{cmd:)} {cmd:level(}#{cmd:)} {cmd:nolog} {cmd:noshow}]

{title:Description}

{pstd}
Fits the Cox model by maximizing the partial likelihood, with the Breslow (default) or Efron approximation for tied failure times; delayed entry and strata enter through the risk sets. Results are shown as hazard ratios ({cmd:nohr} shows coefficients). With {cmd:vce(robust)} the variance uses the efficient score residuals, clustered on the {cmd:stset} id when there is one. The exact methods ({cmd:exactm}, {cmd:exactp}) are not available.

{title:Postestimation}

{phang2}{cmd:predict} newvar [, {cmd:hr} {cmd:xb} {cmd:stdp} {cmd:basesurv} {cmd:basechazard} {cmd:basehc} {cmd:csnell} {cmd:mgale} {cmd:partial}]{p_end}
{phang2}{cmd:estat phtest} [, {cmd:detail} {cmd:log} {cmd:km} {cmd:rank}]{space 3}(Grambsch-Therneau test of proportional hazards){p_end}

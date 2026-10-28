{smcl}
{title:Title}

{p 4 4 2}
{bf:melogit}, {bf:meprobit}, {bf:mepoisson} {hline 2} Multilevel mixed-effects logistic, probit and Poisson regression

{title:Syntax}

{p 8 16 2}
{cmd:melogit} depvar [fe_varlist] [if] [in] {cmd:||} {it:levelvar}{cmd::} [re_varlist] [{cmd:,} {cmd:covariance(}...{cmd:)} {cmd:noconstant}] [{cmd:,} {cmd:or} {cmd:intpoints(}#{cmd:)} {cmd:intmethod(mvaghermite}|{cmd:ghermite}|{cmd:laplace)} {cmd:level(}#{cmd:)} {cmd:nolog}]

{p 8 16 2}
{cmd:mepoisson} depvar [fe_varlist] ... [{cmd:,} {cmd:irr} {cmd:exposure(}varname{cmd:)} {cmd:offset(}varname{cmd:)} ...]

{title:Description}

{pstd}
Fits multilevel generalized linear mixed models with random intercepts and slopes. The likelihood of each group is integrated by mean-variance adaptive Gauss-Hermite quadrature (7 points per random effect by default). The fixed-effects model is fitted first and an LR test against it is reported. Nested levels are written one after the other, the highest first ({cmd:|| school: || class:}); the likelihood is integrated level by level with the points of each level adapted to its groups. Weights: {cmd:fweight}s, {cmd:pweight}s and {cmd:iweight}s at the observation level and the options {cmd:fweight(}varname{cmd:)} and {cmd:pweight(}varname{cmd:)} of a random-effects equation (constant within its groups). {cmd:vce(robust)} and {cmd:vce(cluster} clustvar{cmd:)} give the sandwich variance from the scores of the highest-level groups (the clusters must contain them); sampling weights imply {cmd:vce(robust)}, and no LR test is reported then.

{title:Postestimation}

{phang2}{cmd:predict} newvar [, {cmd:mu} {cmd:xb} {cmd:eta} {cmd:reffects} {cmd:fixedonly}]{p_end}
{phang2}{cmd:estat icc} (melogit, meprobit), {cmd:estat recovariance}{p_end}

{smcl}
{title:Title}

{p 4 4 2}
{bf:mestreg} {hline 2} Multilevel mixed-effects parametric survival models

{title:Syntax}

{p 8 16 2}
{cmd:mestreg} [fe_varlist] [if] [in] {cmd:||} {it:levelvar}{cmd::} [re_varlist] {cmd:,} {cmd:distribution(}exponential|weibull|lognormal|loglogistic|gamma{cmd:)} [{cmd:time} {cmd:nohr} {cmd:tratio} {cmd:covariance(}...{cmd:)} {cmd:noconstant} {cmd:intpoints(}#{cmd:)} {cmd:level(}#{cmd:)}]

{title:Description}

{pstd}
Fits parametric survival models, declared with {help stset}, with a Gaussian random effect per group. Exponential and Weibull are shown as hazard ratios (proportional-hazards metric; {cmd:nohr} for coefficients, {cmd:time} for the accelerated failure-time metric); lognormal, loglogistic and gamma are always in the time metric. The likelihood is integrated by mean-variance adaptive Gauss-Hermite quadrature. An LR test against the model without random effects is reported. Nested levels are written one after the other, the highest first ({cmd:|| school: || class:}); the likelihood is integrated level by level with the points of each level adapted to its groups. Weights: the weights of {help stset} at the observation level and the options {cmd:fweight(}varname{cmd:)} and {cmd:pweight(}varname{cmd:)} of a random-effects equation (constant within its groups). {cmd:vce(robust)} and {cmd:vce(cluster} clustvar{cmd:)} give the sandwich variance from the scores of the highest-level groups (the clusters must contain them); sampling weights imply {cmd:vce(robust)}, and no LR test is reported then.

{title:Postestimation}

{phang2}{cmd:predict} newvar [, {cmd:xb} {cmd:eta} {cmd:reffects} {cmd:median} {cmd:mean} {cmd:hazard}] (the last three given the posterior means of the random effects){p_end}

{title:Limitations}

{pstd}
{cmd:intmethod(laplace)} with more than one level is not available.

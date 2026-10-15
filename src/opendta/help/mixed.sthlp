{smcl}
{title:Title}

{p 4 4 2}
{bf:mixed} {hline 2} Multilevel mixed-effects linear regression

{title:Syntax}

{p 8 16 2}
{cmd:mixed} depvar [fe_varlist] [if] [in] {cmd:||} {it:levelvar}{cmd::} [re_varlist] [{cmd:,} {cmd:covariance(}{cmd:identity}|{cmd:independent}|{cmd:exchangeable}|{cmd:unstructured}{cmd:)} {cmd:noconstant}] [{cmd:||} ...] [{cmd:,} {cmd:mle} {cmd:reml} {cmd:variance} {cmd:stddeviations} {cmd:level(}#{cmd:)} {cmd:nolog} {cmd:nolrtest} {cmd:nogroup}]

{title:Description}

{pstd}
Fits linear models with random intercepts and random coefficients at one or more nested levels ({cmd:|| school: || class:}), by maximum likelihood (default) or restricted maximum likelihood ({cmd:reml}). The fixed effects are shown first; the random-effects parameters (variances, or standard deviations with {cmd:stddeviations}, and covariances or correlations) follow, with an LR test against the linear regression. {cmd:xtmixed} is a synonym. Crossed effects ({cmd:_all:}) and weights are not available.

{title:Postestimation}

{phang2}{cmd:predict} newvar [, {cmd:xb} {cmd:fitted} {cmd:reffects} {cmd:residuals} {cmd:rstandard}]{p_end}
{phang2}{cmd:estat icc}, {cmd:estat group}, {cmd:estat recovariance}{p_end}

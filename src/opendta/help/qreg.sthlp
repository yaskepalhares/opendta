{smcl}
{title:Title}

{p 4 4 2}
{bf:qreg}, {bf:bsqreg}, {bf:iqreg}, {bf:sqreg} {hline 2} Quantile regression

{title:Syntax}

{p 8 16 2}
{cmd:qreg} depvar indepvars [if] [in] [weight] [, quantile(#)]

{p 8 16 2}
{cmd:bsqreg} depvar indepvars [if] [in] [, quantile(#) reps(#) level(#)]

{p 8 16 2}
{cmd:iqreg} depvar indepvars [if] [in] [, quantiles(# #) reps(#) level(#) nodots]

{p 8 16 2}
{cmd:sqreg} depvar indepvars [if] [in] [, quantiles(# [# ...]) reps(#) level(#) nodots]

{title:Description}

{pstd}
Estimates a conditional quantile (default the median) by linear programming. Standard errors assume i.i.d. errors, with the sparsity estimated by the fitted method and the Hall-Sheather bandwidth.
{cmd:bsqreg} gives the same coefficients with bootstrap standard errors ({cmd:reps(20)} by
default); {cmd:iqreg} estimates the difference between two quantiles (.25 and .75 by
default) and {cmd:sqreg} several quantiles at once (equations {cmd:q25}, {cmd:q50}, ...), both
with the bootstrap covariance. OpenDTA's bootstrap draws differ from Stata's, so the standard
errors differ for the same seed.

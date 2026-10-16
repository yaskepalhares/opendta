{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:xtgee} {hline 2}}Population-averaged panel-data models by GEE{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:xtgee} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{it:weight}]
[{cmd:,} {opt f:amily(family)} {opt l:ink(link)} {opt c:orr(correlation)} {opt i(varname)} {opt t(varname)}
{opt vce(robust)} {cmd:eform} {cmd:nmp} {opt sc:ale(x2|dev|phi|#)} {opt exp:osure(varname)}
{opt off:set(varname)} {opt l:evel(#)} {cmd:nolog}]{p_end}

{p 8 16 2}{cmd:xtreg} ..., {cmd:pa}{space 4}{cmd:xtlogit} ..., {cmd:pa}{space 4}{cmd:xtpoisson} ..., {cmd:pa}{p_end}

{title:Description}

{pstd}
{cmd:xtgee} fits generalized linear models for panel data by generalized estimating
equations. The mean follows a GLM ({it:family} {cmd:gaussian}, {cmd:binomial}, {cmd:poisson},
{cmd:nbinomial}, {cmd:gamma} or {cmd:igaussian}, each with its canonical {it:link} unless
{cmd:link()} names {cmd:identity}, {cmd:log}, {cmd:logit}, {cmd:probit}, {cmd:cloglog} or
{cmd:power} {it:#}), and the observations of a panel are correlated through a working
correlation matrix: {cmd:independent}, {cmd:exchangeable} (the default), {cmd:ar} {it:#},
{cmd:stationary} {it:#}, {cmd:nonstationary} {it:#} or {cmd:unstructured}. The last four need
the time variable ({cmd:xtset} {it:panelvar timevar} or {cmd:t()}).
{p_end}

{pstd}
{cmd:xtreg, pa}, {cmd:xtlogit, pa} and {cmd:xtpoisson, pa} are {cmd:xtgee} with the Gaussian,
binomial (logit) and Poisson families. The standard errors are model based; {cmd:vce(robust)}
gives the sandwich estimator clustered on the panels. The working correlation is stored in
{cmd:e(R)}. {cmd:predict} offers {cmd:mu} (default), {cmd:xb} and {cmd:stdp}.
{p_end}

{title:Examples}

{phang2}{cmd:. xtset id year}{p_end}
{phang2}{cmd:. xtgee y x1 x2, family(binomial) link(logit) corr(exchangeable) vce(robust)}{p_end}
{phang2}{cmd:. xtreg y x1 x2, pa corr(ar 1)}{p_end}

{title:Also see}

{psee}{help xtreg}, {help xtlogit}, {help xtpoisson}, {help glm}{p_end}

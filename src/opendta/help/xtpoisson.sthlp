{smcl}
{title:Title}

{p 4 4 2}
{bf:xtpoisson} {hline 2} Fixed-effects and random-effects Poisson models

{title:Syntax}

{p 8 16 2}
{cmd:xtpoisson} depvar indepvars [if] [in] [, fe re pa irr]

{title:Description}

{pstd}
{cmd:fe} fits the conditional fixed-effects Poisson model (groups with all-zero outcomes are dropped). {cmd:re} (default) fits the gamma random-effects model and reports the LR test of alpha=0. {cmd:pa} fits the population-averaged model by GEE with an exchangeable correlation ({help xtgee} options such as {cmd:corr()} and {cmd:vce(robust)} are accepted).

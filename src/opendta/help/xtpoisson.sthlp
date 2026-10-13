{smcl}
{title:Title}

{p 4 4 2}
{bf:xtpoisson} {hline 2} Fixed-effects and random-effects Poisson models

{title:Syntax}

{p 8 16 2}
{cmd:xtpoisson} depvar indepvars [if] [in] [, fe re irr]

{title:Description}

{pstd}
{cmd:fe} fits the conditional fixed-effects Poisson model (groups with all-zero outcomes are dropped). {cmd:re} (default) fits the gamma random-effects model and reports the LR test of alpha=0.

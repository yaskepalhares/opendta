{smcl}
{title:Title}

{p 4 4 2}
{bf:qreg} {hline 2} Quantile regression

{title:Syntax}

{p 8 16 2}
{cmd:qreg} depvar indepvars [if] [in] [weight] [, quantile(#)]

{title:Description}

{pstd}
Estimates a conditional quantile (default the median) by linear programming. Standard errors assume i.i.d. errors, with the sparsity estimated by the fitted method and the Hall-Sheather bandwidth.

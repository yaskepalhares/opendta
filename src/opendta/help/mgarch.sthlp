{smcl}
{title:Title}

{p 4 4 2}
{bf:mgarch} {hline 2} Multivariate GARCH models

{title:Syntax}

{p 8 16 2}{cmd:mgarch} {c -(}{opt dvech}|{opt ccc}|{opt dcc}{c )-} {cmd:(}{it:depvars} [{cmd:=} {it:indepvars}] [{cmd:,} {opt noconstant}]{cmd:)} [...] [{it:if}] [{it:in}] [{cmd:,} {opt arch(numlist)} {opt garch(numlist)} {opt level(#)} {opt nolog}]{p_end}

{title:Description}

{pstd}
Fits diagonal vech, constant conditional correlation and dynamic conditional correlation models with Gaussian errors, by maximum likelihood (BHHH then BFGS) with the observed information covariance. The presample variances are the sample averages of the squared residuals (VERIFY). Not available: {cmd:mgarch vcc}, {cmd:distribution(t)}, {cmd:het()}, {cmd:constraints()}, {cmd:from()}, {cmd:vce(robust)}, {cmd:predict} and {cmd:estat} after it.
{p_end}

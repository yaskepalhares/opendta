{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:tobit} {hline 2}}Tobit regression{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:tobit} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{it:weight}] {cmd:,} {opt ll}[{cmd:(}{it:#}{cmd:)}] {opt ul}[{cmd:(}{it:#}{cmd:)}] [{it:options}]{p_end}

{title:Description}

{pstd}
Linear regression with a censored outcome. {opt ll} and {opt ul} give the censoring limits (without a number, the minimum or maximum of {it:depvar}). The output shows {cmd:/sigma} and the counts of censored observations.
{p_end}

{title:Options}

{phang}
{opt l:evel(#)} confidence level; default is {cmd:level(95)}.
{p_end}

{title:Examples}

{phang2}{cmd:. tobit wage educ exper, ll(0)}{p_end}

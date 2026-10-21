{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:xtfrontier} {hline 2}}Panel data: stochastic frontier models{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:xtfrontier} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {opt ti}|{opt tvd} {opt cost} {it:options}]{p_end}

{title:Description}

{pstd}
{cmd:xtfrontier} fits stochastic production (or cost) frontier models for panel data by maximum likelihood: the time-invariant inefficiency model ({cmd:ti}) and the time-varying decay model ({cmd:tvd}) of Battese and Coelli, with a truncated-normal inefficiency term.
{p_end}

{title:Options}

{phang}
{opt cost} fit a cost frontier.
{p_end}

{phang}
{opt ti} time-invariant inefficiency model; one of {opt ti} or {opt tvd} is required.
{p_end}

{phang}
{opt tvd} time-varying decay model, with parameter {cmd:eta}.
{p_end}

{phang}
{opt l:evel(#)} confidence level; default is {cmd:level(95)}.
{p_end}

{title:Examples}

{phang2}{cmd:. xtset id year}{p_end}
{phang2}{cmd:. xtfrontier ly lk ll}{p_end}
{phang2}{cmd:. xtfrontier ly lk ll, tvd}{p_end}

{title:Also see}

{pstd}
See also {help xtreg}.
{p_end}

{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:xtquad} {hline 2}}Panel data: random-effects models fit by quadrature{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:xtprobit} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {opt re}|{opt pa} {it:options}]{p_end}
{p 8 16 2}{cmd:xtcloglog} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {opt re}|{opt pa} {it:options}]{p_end}
{p 8 16 2}{cmd:xtologit} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:xtoprobit} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:xttobit} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {opt ll(#)} {opt ul(#)} {it:options}]{p_end}
{p 8 16 2}{cmd:xtintreg} {it:depvar1} {it:depvar2} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:xtnbreg} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {opt re}|{opt fe}|{opt pa} {opt irr} {it:options}]{p_end}
{p 8 16 2}{cmd:quadchk} [{it:#1} [{it:#2}]] [{cmd:,} {opt nooutput}]{p_end}

{title:Description}

{pstd}
These commands fit panel-data models with a Gaussian random effect, integrated by Gauss-Hermite quadrature: adaptive ({cmd:intmethod(mvaghermite)}, the default) or non-adaptive ({cmd:intmethod(ghermite)}), with 12 points by default. {cmd:xtnbreg} fits the Hausman, Hall and Griliches negative binomial models (random effects with a Beta distribution, or conditional fixed effects).
{p_end}

{pstd}
{cmd:quadchk} refits the last model with more and fewer quadrature points and compares the log likelihood and the coefficients, to check that the approximation is stable.
{p_end}

{title:Options}

{phang}
{opt intp:oints(#)} number of quadrature points; default is 12.
{p_end}

{phang}
{opt intm:ethod(intmethod)} {cmd:mvaghermite} or {cmd:ghermite}.
{p_end}

{phang}
{opt vce(vcetype)} {cmd:oim}, {cmd:robust} or {cmd:cluster} {it:clustvar}.
{p_end}

{phang}
{opt l:evel(#)} confidence level; default is {cmd:level(95)}.
{p_end}

{title:Examples}

{phang2}{cmd:. xtset id year}{p_end}
{phang2}{cmd:. xtprobit union age grade, re}{p_end}
{phang2}{cmd:. quadchk}{p_end}
{phang2}{cmd:. xtologit rating age, intpoints(8)}{p_end}

{title:Also see}

{pstd}
See also {help xtlogit}, {help xtpoisson}, {help xtgee} and {help xtreg}.
{p_end}

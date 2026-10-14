{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:regress} {hline 2}}Linear regression{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:regress} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:areg} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{it:weight}] {cmd:,} {opt a:bsorb(varname)} [{it:options}]{p_end}
{p 8 16 2}{cmd:ivregress} {cmd:2sls} {it:depvar} [{it:varlist1}] {cmd:(}{it:varlist2} {cmd:=} {it:instruments}{cmd:)} [{it:if}] [{it:in}] [{cmd:,} {it:options}]{p_end}

{title:Description}

{pstd}
{cmd:regress} fits a linear model by ordinary least squares. Factor variables ({cmd:i.}, {cmd:c.}, {cmd:#}, {cmd:##}) and time-series operators are allowed. {cmd:aweight}s, {cmd:fweight}s, {cmd:iweight}s and {cmd:pweight}s are accepted.
{p_end}

{pstd}
{cmd:areg} absorbs one categorical variable (fixed effects) without showing its coefficients. {cmd:ivregress 2sls} fits two-stage least squares.
{p_end}

{title:Options}

{phang}
{opt nocons:tant} suppress the constant term.
{p_end}

{phang}
{opt vce(vcetype)} {cmd:oim}, {cmd:robust}, {cmd:cluster} {it:clustvar} or {cmd:bootstrap}.
{p_end}

{phang}
{opt l:evel(#)} confidence level; default is {cmd:level(95)}.
{p_end}

{phang}
{opt beta} report standardized beta coefficients instead of confidence intervals.
{p_end}

{title:Examples}

{phang2}{cmd:. regress price mpg weight}{p_end}
{phang2}{cmd:. regress price mpg i.foreign, vce(robust)}{p_end}
{phang2}{cmd:. areg y x, absorb(firm)}{p_end}
{phang2}{cmd:. ivregress 2sls y x1 (x2 = z1 z2)}{p_end}

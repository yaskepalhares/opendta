{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:poisson} {hline 2}}Poisson and negative binomial regression{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:poisson} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:nbreg} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {it:options}]{p_end}

{title:Description}

{pstd}
Count-data models. {cmd:nbreg} adds the overdispersion parameter alpha ({cmd:/lnalpha} and {cmd:alpha} in the output) and the likelihood-ratio test of alpha=0.
{p_end}

{title:Options}

{phang}
{opt e:xposure(varname)} include ln({it:varname}) with coefficient 1.
{p_end}

{phang}
{opt off:set(varname)} include {it:varname} with coefficient 1.
{p_end}

{phang}
{opt irr} report incidence-rate ratios.
{p_end}

{phang}
{opt vce(vcetype)} {cmd:oim}, {cmd:robust}, {cmd:cluster} {it:clustvar} or {cmd:bootstrap}.
{p_end}

{phang}
{opt l:evel(#)} confidence level; default is {cmd:level(95)}.
{p_end}

{title:Remarks}

{pstd}
nbreg starts from a Poisson fit; the iteration log may differ from Stata's, the estimates do not.
{p_end}

{title:Examples}

{phang2}{cmd:. poisson deaths smokes i.agecat, exposure(pyears) irr}{p_end}
{phang2}{cmd:. nbreg daysabs math i.prog}{p_end}

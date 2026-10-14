{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:logit} {hline 2}}Logistic and probit regression for binary outcomes{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:logit} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt or} {it:options}]{p_end}
{p 8 16 2}{cmd:logistic} {it:depvar} {it:indepvars} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:probit} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {it:options}]{p_end}

{title:Description}

{pstd}
Maximum-likelihood models for 0/1 outcomes: any nonzero, nonmissing value of {it:depvar} is a success. {cmd:logistic} reports odds ratios; {cmd:logit, or} does the same. The iteration log, the header and the coefficient table follow the layout of Stata 14.
{p_end}

{title:Options}

{phang}
{opt or} report odds ratios (logit).
{p_end}

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
{opt nolog} do not show the iteration log.
{p_end}

{title:Remarks}

{pstd}
After estimation: {help predict}, {help margins}, {help test}, {help estat} ({cmd:estat classification}, {cmd:estat gof}, {cmd:estat ic}).
{p_end}

{title:Examples}

{phang2}{cmd:. logit low age smoke}{p_end}
{phang2}{cmd:. logistic low age i.race smoke}{p_end}
{phang2}{cmd:. probit foreign weight mpg, vce(robust)}{p_end}

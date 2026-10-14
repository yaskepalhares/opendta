{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:glm} {hline 2}}Generalized linear models{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:glm} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt f:amily(family)} {opt l:ink(link)} {it:options}]{p_end}

{title:Description}

{pstd}
Fits generalized linear models by maximum likelihood (Newton-Raphson with the observed information matrix, as Stata's default). Families: {cmd:gaussian}, {cmd:binomial} [{it:#}|{it:varname}], {cmd:poisson}, {cmd:nbinomial} [{it:#}], {cmd:gamma}, {cmd:igaussian}. Links: {cmd:identity}, {cmd:log}, {cmd:logit}, {cmd:probit}, {cmd:cloglog}, {cmd:power} {it:#}, {cmd:opower} {it:#}, {cmd:nbinomial}, {cmd:loglog}, {cmd:logc}. Without {opt link()}, the canonical link of the family is used.
{p_end}

{title:Options}

{phang}
{opt eform} report exponentiated coefficients.
{p_end}

{phang}
{opt e:xposure(varname)} {opt off:set(varname)} as in {help poisson}.
{p_end}

{phang}
{opt vce(vcetype)} {cmd:oim}, {cmd:robust}, {cmd:cluster} {it:clustvar} or {cmd:bootstrap}.
{p_end}

{phang}
{opt l:evel(#)} confidence level; default is {cmd:level(95)}.
{p_end}

{title:Examples}

{phang2}{cmd:. glm y x, family(binomial) link(log) eform}{p_end}
{phang2}{cmd:. glm deaths age, family(poisson) exposure(pyears)}{p_end}

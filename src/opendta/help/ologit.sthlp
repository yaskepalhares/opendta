{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:ologit} {hline 2}}Ordered and multinomial outcomes{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:ologit} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt or} {it:options}]{p_end}
{p 8 16 2}{cmd:oprobit} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:mlogit} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt b:aseoutcome(#)} {opt rrr} {it:options}]{p_end}

{title:Description}

{pstd}
{cmd:ologit} and {cmd:oprobit} fit ordered logit and probit models and report the cutpoints ({cmd:/cut1}, {cmd:/cut2}, ...). {cmd:mlogit} fits a multinomial logit; by default the most frequent outcome is the base.
{p_end}

{title:Options}

{phang}
{opt or} odds ratios (ologit).
{p_end}

{phang}
{opt rrr} relative-risk ratios (mlogit).
{p_end}

{phang}
{opt b:aseoutcome(#)} value of {it:depvar} used as the base outcome.
{p_end}

{phang}
{opt vce(vcetype)} {cmd:oim}, {cmd:robust}, {cmd:cluster} {it:clustvar} or {cmd:bootstrap}.
{p_end}

{phang}
{opt l:evel(#)} confidence level; default is {cmd:level(95)}.
{p_end}

{title:Examples}

{phang2}{cmd:. ologit rep78 price mpg}{p_end}
{phang2}{cmd:. mlogit insure age male, baseoutcome(1) rrr}{p_end}

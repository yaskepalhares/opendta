{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:mean} {hline 2}}Estimate means, proportions, ratios and totals{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:mean} {it:varlist} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt over(varlist)} {opt l:evel(#)} {opt vce()}]{p_end}
{p 8 16 2}{cmd:proportion} {it:varlist} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt over(varlist)} {opt l:evel(#)}]{p_end}
{p 8 16 2}{cmd:ratio} [{it:name}{cmd::}] {it:varname} {cmd:/} {it:varname} [{cmd:,} {opt over(varlist)}]{p_end}
{p 8 16 2}{cmd:total} {it:varlist} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt over(varlist)}]{p_end}

{title:Description}

{pstd}
Estimation commands for descriptive statistics: the results are stored in {cmd:e()} and accept {help test}, {help lincom} and {help nlcom}. With the {help svy} prefix they use the survey design.
{p_end}

{title:Examples}

{phang2}{cmd:. mean bmi, over(sex)}{p_end}
{phang2}{cmd:. proportion smoker}{p_end}
{phang2}{cmd:. svy: mean bmi, over(region)}{p_end}

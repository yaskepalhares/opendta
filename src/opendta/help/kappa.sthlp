{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:kappa} {hline 2}}Interrater agreement{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:kap} {it:varname1} {it:varname2} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt tab} {opt w:gt(w | w2 | matname)} {opt abs:olute}]{p_end}
{p 8 16 2}{cmd:kap} {it:varname1} {it:varname2} {it:varname3} [...] [{it:if}] [{it:in}] [{it:weight}]{p_end}
{p 8 16 2}{cmd:kappa} {it:varlist} [{it:if}] [{it:in}]{p_end}

{title:Description}

{pstd}
{cmd:kap} with two variables (one per rater) computes Cohen's kappa, the observed and expected agreement and the test that kappa is zero; {opt wgt(w)} weights disagreements linearly and {opt wgt(w2)} quadratically, by the order of the categories or, with {opt absolute}, by their values; {opt tab} shows the table of ratings. With three or more variables, each one a rater, {cmd:kap} computes Fleiss's kappa for each outcome and combined.
{p_end}

{pstd}
{cmd:kappa} takes the ratings in frequency form: one variable per outcome, holding how many raters assigned the subject to it.
{p_end}

{title:Examples}

{phang2}{cmd:. kap rater1 rater2}{p_end}
{phang2}{cmd:. kap rater1 rater2, wgt(w2) tab}{p_end}
{phang2}{cmd:. kap rater1 rater2 rater3}{p_end}
{phang2}{cmd:. kappa cat1 cat2 cat3}{p_end}

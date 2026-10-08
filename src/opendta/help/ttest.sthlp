{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:ttest} {hline 2}}Tests and confidence intervals for means and proportions{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:ttest} {it:varname} {cmd:==} {it:#} [{it:if}] [{it:in}] [{cmd:,} {opt l:evel(#)}]{p_end}
{p 8 16 2}{cmd:ttest} {it:varname1} {cmd:==} {it:varname2} [{cmd:,} {opt unp:aired} {opt une:qual} {opt w:elch}]{p_end}
{p 8 16 2}{cmd:ttest} {it:varname} [{cmd:,} {opt by(groupvar)} {opt une:qual} {opt w:elch}]{p_end}
{p 8 16 2}{cmd:prtest} {it:varname} {cmd:==} {it:#}  |  {cmd:prtest} {it:varname} [{cmd:,} {opt by(groupvar)}]{p_end}
{p 8 16 2}{cmd:ci} [{it:varlist}] [{cmd:,} {opt l:evel(#)} {opt b:inomial}]{p_end}

{title:Description}

{pstd}
{cmd:ttest} compares a mean with a value, two paired variables or the means
of two groups ({opt unequal} uses Satterthwaite's degrees of freedom,
{opt welch} Welch's). The output shows the three alternative hypotheses.
{cmd:prtest} does the same for proportions of 0/1 variables using the
normal approximation. {cmd:ci} gives confidence intervals for means, or
exact (Clopper-Pearson) intervals for proportions with {opt binomial}.
{p_end}

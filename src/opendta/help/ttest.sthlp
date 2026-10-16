{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:ttest} {hline 2}}Tests and confidence intervals for means and proportions{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:ttest} {it:varname} {cmd:==} {it:#} [{it:if}] [{it:in}] [{cmd:,} {opt l:evel(#)}]{p_end}
{p 8 16 2}{cmd:ttest} {it:varname1} {cmd:==} {it:varname2} [{cmd:,} {opt unp:aired} {opt une:qual} {opt w:elch}]{p_end}
{p 8 16 2}{cmd:ttest} {it:varname} [{cmd:,} {opt by(groupvar)} {opt une:qual} {opt w:elch}]{p_end}
{p 8 16 2}{cmd:prtest} {it:varname} {cmd:==} {it:#}  |  {cmd:prtest} {it:varname1} {cmd:==} {it:varname2}  |  {cmd:prtest} {it:varname} {cmd:,} {opt by(groupvar)}{p_end}
{p 8 16 2}{cmd:ci} [{cmd:means}] [{it:varlist}] [{cmd:,} {opt l:evel(#)} {opt poisson} {opt exp:osure(varname)}]{p_end}
{p 8 16 2}{cmd:ci proportions} [{it:varlist}] [{cmd:,} {opt l:evel(#)} {opt exact} {opt wald} {opt wilson} {opt agresti} {opt jeffreys}]{p_end}
{p 8 16 2}{cmd:ci variances} [{it:varlist}] [{cmd:,} {opt l:evel(#)} {opt sd} {opt bonett}]{p_end}

{title:Description}

{pstd}
{cmd:ttest} compares a mean with a value, two paired variables or the means
of two groups ({opt unequal} uses Satterthwaite's degrees of freedom,
{opt welch} Welch's). The output shows the three alternative hypotheses.
{cmd:prtest} does the same for proportions of 0/1 variables using the
normal approximation: one proportion against a value, two variables, or two
groups with {opt by()}. {cmd:ci} gives confidence intervals for means
(exact Poisson intervals with {opt poisson}), for proportions (exact
Clopper-Pearson by default; Wald, Wilson, Agresti-Coull or Jeffreys) and
for variances (chi-squared, or Bonett's interval, which does not assume
normality, with {opt bonett}; {opt sd} reports standard deviations).
{cmd:ci binomial} is the old form of {cmd:ci proportions}.
{p_end}

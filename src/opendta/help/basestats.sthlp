{smcl}
{title:Title}

{p 4 4 2}
{bf:ameans}, {bf:icc}, {bf:loneway}, {bf:esize}, {bf:esizei}, {bf:lnskew0}, {bf:ladder}, {bf:boxcox} {hline 2} Means, agreement, effect sizes and transformations

{title:Syntax}

{p 8 16 2}
{cmd:ameans} {it:varlist} [{cmd:,} {cmd:level(}#{cmd:)}]

{p 8 16 2}
{cmd:icc} {it:depvar target} [{it:rater}] [{cmd:,} {cmd:mixed} {cmd:oneway} {cmd:absolute} {cmd:consistency} {cmd:level(}#{cmd:)}]

{p 8 16 2}
{cmd:loneway} {it:response_var group_var} [{cmd:,} {cmd:level(}#{cmd:)}]

{p 8 16 2}
{cmd:esize twosample} {it:varname} {cmd:,} {cmd:by(}{it:groupvar}{cmd:)} [{cmd:cohensd} {cmd:hedgesg} {cmd:glassdelta} {cmd:pbcorr} {cmd:all}]   {cmd:esizei} {it:n1 m1 sd1 n2 m2 sd2}

{p 8 16 2}
{cmd:lnskew0} {it:newvar} {cmd:=} {it:exp} [{cmd:,} {cmd:level(}#{cmd:)}]   {cmd:ladder} {it:varname} [{cmd:,} {cmd:generate(}{it:newvar}{cmd:)} {cmd:noadjust}]   {cmd:gladder} {it:varname}

{p 8 16 2}
{cmd:boxcox} {it:depvar} [{it:indepvars}] [{cmd:,} {cmd:model(lhsonly)} {cmd:level(}#{cmd:)}]

{title:Description}

{pstd}
{cmd:ameans} reports arithmetic, geometric and harmonic means with confidence intervals. {cmd:icc} estimates intraclass correlations (McGraw and Wong) with F-based intervals. {cmd:loneway} does one-way ANOVA and the intraclass correlation. {cmd:esize} reports Cohen's d, Hedges's g, Glass's Delta and the point-biserial correlation with noncentral-t intervals. {cmd:lnskew0} finds k such that ln(x-k) has zero skewness. {cmd:ladder} tests each power of the ladder with D'Agostino's normality test. {cmd:boxcox} fits the left-hand-side Box-Cox model by maximum likelihood.

{title:Limitations}

{pstd}
{cmd:ameans} does not take weights; {cmd:icc, testvalue()}, {cmd:loneway, median} and {cmd:exact}, {cmd:esize, unequal} and {cmd:boxcox} models other than {cmd:lhsonly} give errors. {cmd:icc} needs a balanced table.

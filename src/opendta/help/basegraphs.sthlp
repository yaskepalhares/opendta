{smcl}
{title:Title}

{p 4 4 2}
{bf:quantile}, {bf:symplot}, {bf:serrbar}, {bf:dotplot}, {bf:sunflower}, {bf:lpoly}, {bf:cusum}, {bf:gladder} {hline 2} Distribution and smoothing graphs

{title:Syntax}

{p 8 16 2}
{cmd:quantile} {it:varname}   {cmd:symplot} {it:varname}   {cmd:gladder} {it:varname}

{p 8 16 2}
{cmd:serrbar} {it:mvar svar xvar} [{cmd:,} {cmd:scale(}#{cmd:)}]

{p 8 16 2}
{cmd:dotplot} {it:varlist} [{cmd:,} {cmd:over(}{it:varname}{cmd:)} {cmd:ny(}#{cmd:)}]

{p 8 16 2}
{cmd:sunflower} {it:yvar xvar} [{cmd:,} {cmd:binwidth(}#{cmd:)} {cmd:light(}#{cmd:)} {cmd:dark(}#{cmd:)} {cmd:flowersonly}]

{p 8 16 2}
{cmd:lpoly} {it:yvar xvar} [{cmd:,} {cmd:kernel(}{it:kernel}{cmd:)} {cmd:degree(}#{cmd:)} {cmd:bwidth(}#{cmd:)} {cmd:n(}#{cmd:)} {cmd:at(}{it:varname}{cmd:)} {cmd:generate(}{it:newvar}{cmd:)} {cmd:nograph} {cmd:noscatter}]

{p 8 16 2}
{cmd:cusum} {it:yvar xvar} [{cmd:,} {cmd:generate(}{it:newvar}{cmd:)}]

{title:Description}

{pstd}
{cmd:quantile} plots the ordered values against the fraction of the data and {cmd:symplot} the distances above the median against the distances below it. {cmd:serrbar} plots a mean with error bars. {cmd:dotplot} draws a dot plot. {cmd:sunflower} draws a density-distribution sunflower plot. {cmd:lpoly} draws a local polynomial smooth. {cmd:cusum} plots the cumulative sum of the standardized residuals of the regression of {it:yvar} on {it:xvar}. {cmd:gladder} shows the histograms of the ladder of powers in a 3x3 graph.

{title:Limitations}

{pstd}
The default bandwidth of {cmd:lpoly} is a reference-rule value, not Stata's plug-in estimator; {cmd:cusum} does not display the test statistics. The drawing details are marked VERIFY.

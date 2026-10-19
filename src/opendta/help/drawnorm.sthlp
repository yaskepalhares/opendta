{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:drawnorm} {hline 2}}Drawing multivariate normal data{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:drawnorm} {it:newvarlist} [{cmd:,} {opt clear} {opt double} {opt n(#)} {opt means(vector)} {opt sds(vector)} {opt corr(matrix)} | {opt cov(matrix)} {opt cstor:age(full|lower|upper)} {opt forcepsd} {opt seed(#)}]{p_end}
{p 8 16 2}{cmd:corr2data} {it:newvarlist} [{cmd:,} {it:the same options}]{p_end}

{title:Description}

{pstd}
{cmd:drawnorm} draws a sample from a multivariate normal distribution with the given means, standard deviations and correlations (or covariances). {cmd:corr2data} creates data whose sample means, standard deviations and correlations are exactly the ones given, useful to reproduce published analyses from summary statistics.
{p_end}

{title:Remarks}

{pstd}
The values follow the requested distribution but are not the same numbers Stata draws for the same seed.
{p_end}

{title:Examples}

{phang2}{cmd:. matrix C = (1, .5 \ .5, 1)}{p_end}
{phang2}{cmd:. drawnorm x y, n(500) corr(C) means(10 20) sds(2 4)}{p_end}
{phang2}{cmd:. corr2data a b, n(100) corr(C)}{p_end}

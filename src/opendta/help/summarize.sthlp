{smcl}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:summarize} {hline 2}}Summary statistics{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:summarize} [{it:varlist}] [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt d:etail} {opt mean:only} {opt f:ormat} {opt sep:arator(#)}]

{p 8 16 2}
{cmd:tabstat} {it:varlist} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt s:tatistics(list)} {opt by(varname)} {opt c:olumns(variables|statistics)} {opt not:otal} {opt case:wise}]

{p 8 16 2}
{cmd:centile} [{it:varlist}] [{cmd:,} {opt c:entile(numlist)}]{p_end}
{p 8 16 2}
{cmd:pctile} {it:newvar} {cmd:=} {it:exp} [{cmd:,} {opt nq(#)}]{p_end}
{p 8 16 2}
{cmd:xtile} {it:newvar} {cmd:=} {it:exp} [{cmd:,} {opt nq(#)}]{p_end}
{p 8 16 2}
{cmd:correlate} [{it:varlist}] [{cmd:,} {opt cov:ariance}]{p_end}
{p 8 16 2}
{cmd:pwcorr} [{it:varlist}] [{cmd:,} {opt sig} {opt obs}]

{title:Description}

{pstd}
{cmd:summarize} shows the number of observations, mean, standard deviation,
minimum and maximum. With {opt detail} it adds percentiles, the four
smallest and largest values, variance, skewness and kurtosis.
{opt meanonly} computes the results silently. All of them are stored in
{cmd:r()}.
{p_end}

{pstd}
{cmd:fweight}s, {cmd:aweight}s and {cmd:iweight}s are allowed.
Percentiles follow the weighted definition: with total weight W, the
{it:p}th percentile is the first ordered value whose cumulative weight
exceeds W{it:p}/100, or the average of two neighbours when it is equal.
{p_end}

{pstd}
{cmd:tabstat} builds compact tables of statistics, optionally by groups.
{cmd:centile} uses the interpolated definition (rank ({it:n}+1){it:p}/100).
{cmd:pctile} creates the cut points of quantile groups and {cmd:xtile}
the group number of each observation.
{p_end}

{title:Also see}

{pstd}
{help ttest}, {help tabulate}, {help egen}, {help collapse}
{p_end}

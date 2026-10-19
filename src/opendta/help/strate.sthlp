{smcl}
{title:Title}

{p 4 4 2}
{bf:strate} {hline 2} Rates, person-time, confidence intervals and tables for survival data

{title:Syntax}

{p 8 16 2}{cmd:strate} [{it:varlist}] [if] [in] [{cmd:,} {cmd:per(}#{cmd:)} {cmd:smr(}{it:varname}{cmd:)} {cmd:jackknife} {cmd:cluster(}{it:varname}{cmd:)} {cmd:output(}{it:filename}{cmd:)} {cmd:level(}#{cmd:)} {cmd:noshow}]{p_end}
{p 8 16 2}{cmd:stptime} [if] [in] [{cmd:,} {cmd:per(}#{cmd:)} {cmd:by(}{it:varlist}{cmd:)} {cmd:title(}{it:string}{cmd:)} {cmd:at(}{it:numlist}{cmd:)} {cmd:smr(}{it:varname}{cmd:)} {cmd:dd(}#{cmd:)} {cmd:jackknife} {cmd:cluster(}{it:varname}{cmd:)} {cmd:output(}{it:filename}{cmd:)} {cmd:level(}#{cmd:)} {cmd:noshow}]{p_end}
{p 8 16 2}{cmd:stci} [if] [in] [{cmd:,} {cmd:by(}{it:varlist}{cmd:)} {cmd:median} {cmd:p(}#{cmd:)} {cmd:rmean} {cmd:emean} {cmd:ccorr} {cmd:level(}#{cmd:)} {cmd:noshow}]{p_end}
{p 8 16 2}{cmd:stir} {it:exposurevar} [if] [in] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:stmh} {it:exposurevar} [{it:strata}] [if] [in] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:stmc} {it:exposurevar} [{it:strata}] [if] [in] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:ltable} {it:timevar} [{it:deathvar}] [{it:freqvar}] [if] [in] [{cmd:,} {cmd:by(}{it:varlist}{cmd:)} {cmd:intervals(}{it:numlist}{cmd:)} {cmd:survival} {cmd:failure} {cmd:hazard} {cmd:noadjust} {cmd:notable} {cmd:test} {cmd:level(}#{cmd:)}]{p_end}

{title:Description}

{pstd}
{cmd:strate} tabulates the incidence rate (failures per person-time) over the categories of {it:varlist}, with an exact Poisson confidence interval; {cmd:smr()} adds standardized mortality ratios against expected rates and {cmd:jackknife} or {cmd:cluster()} give jackknife limits. {cmd:output()} saves the table to a file. {cmd:stptime} reports person-time, failures and the rate, optionally by group and compared with reference rates. {cmd:stci} reports the confidence interval of the mean or of the median (or another percentile) survival time; {cmd:rmean} and {cmd:emean} give the restricted and the extended mean. {cmd:stir}, {cmd:stmh} and {cmd:stmc} report incidence-rate ratios: crude, Mantel-Haenszel stratified, and Mantel-Cox. {cmd:ltable} reports the life table for grouped data: interval, beginning total, deaths, lost, survival, standard error and confidence interval.

{pstd}
{cmd:graph}, {cmd:saving()} and {cmd:tvid()} of {cmd:ltable} and {cmd:graph} of {cmd:strate} and {cmd:stci} are not available yet and stop with error 198. Table layouts that the manual does not settle are marked VERIFY.

{title:Stored results}

{pstd}
{cmd:stptime} stores {cmd:r(ptime)}, {cmd:r(failures)}, {cmd:r(rate)}, {cmd:r(lb)} and {cmd:r(ub)}.

{title:See also}

{p 4 4 2}{help stset}, {help sts}, {help stsplit}, {help epitab}

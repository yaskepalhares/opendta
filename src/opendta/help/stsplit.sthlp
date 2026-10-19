{smcl}
{title:Title}

{p 4 4 2}
{bf:stsplit} {hline 2} Manipulate survival-time records: split, join, fill, generate, vary, base, snapspan

{title:Syntax}

{p 8 16 2}{cmd:stsplit} [{it:newvar}] [if] [in] [{cmd:,} {cmd:at(}{it:numlist}|{cmd:failures)} {cmd:every(}#{cmd:)} {cmd:after(}{it:spec}{cmd:)} {cmd:trim} {cmd:strata(}{it:varlist}{cmd:)} {cmd:riskset(}{it:newvar}{cmd:)} {cmd:nopreserve}]{p_end}
{p 8 16 2}{cmd:stjoin} [{cmd:,} {cmd:censored(}{it:varlist}{cmd:)}]{p_end}
{p 8 16 2}{cmd:stfill} {it:varlist} [if] [in] [{cmd:,} {cmd:baseline} {cmd:forward}]{p_end}
{p 8 16 2}{cmd:stvary} [{it:varlist}] [if] [in] [{cmd:,} {cmd:noshow}]{p_end}
{p 8 16 2}{cmd:stgen} {it:newvar} {cmd:=} {it:function}{cmd:(}{it:...}{cmd:)} [{cmd:,} {cmd:replace}]{p_end}
{p 8 16 2}{cmd:stbase} [{cmd:,} {cmd:at(}{it:...}{cmd:)} {cmd:gap(}#{cmd:)} {cmd:replace} {cmd:noshow} {cmd:nopreserve}]{p_end}
{p 8 16 2}{cmd:snapspan} {it:idvar} {it:timevar} {it:varlist} [{cmd:,} {cmd:generate(}{it:newvars}{cmd:)} {cmd:replace}]{p_end}

{title:Description}

{pstd}
These commands work on data that have been {cmd:stset}. {cmd:stsplit} splits each record into episodes at the given analysis times, at every {cmd:every()} units, or at the failure times ({cmd:at(failures)}), creating {it:newvar} with the start of each episode; {cmd:stjoin} joins consecutive episodes of a subject that agree on all the other variables; {cmd:stfill} fills in a time-varying variable from its first (or last, with {cmd:forward}) nonmissing value; {cmd:stvary} reports whether variables vary within subject; {cmd:stgen} generates variables from the history of each subject; {cmd:stbase} reshapes the data to one record per subject at a baseline; {cmd:snapspan} converts snapshot data to time-span records.

{pstd}
Options that are not implemented yet stop with error 198 ("not available yet in OpenDTA"). Details of the layout that the manual does not settle are marked VERIFY in the code and tested in compat/do/1901_survival.do.

{title:See also}

{p 4 4 2}{help stset}, {help sts}, {help stcox}, {help streg}, {help strate}

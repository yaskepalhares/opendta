{smcl}
{title:Title}

{p 4 4 2}
{bf:ct} {hline 2} Count-time data: ctset, ct, cttost, sttoct, sttocc

{title:Syntax}

{p 8 16 2}{cmd:ctset} {it:timevar} {it:nfailvar} {it:ncensvar} {it:nenterfvar} [{cmd:,} {cmd:by(}{it:varlist}{cmd:)} {cmd:noshow}]{p_end}
{p 8 16 2}{cmd:ct} [{cmd:,} {cmd:show}]{p_end}
{p 8 16 2}{cmd:cttost} [{cmd:,} {cmd:t0(}{it:varname}{cmd:)} {cmd:wvar(}{it:varname}{cmd:)} {cmd:clear} {cmd:nopreserve}]{p_end}
{p 8 16 2}{cmd:sttoct} [{cmd:,} {cmd:by(}{it:varlist}{cmd:)} {cmd:replace} {cmd:noshow}]{p_end}
{p 8 16 2}{cmd:sttocc} [{it:varlist}] [{cmd:,} {cmd:match(}{it:varlist}{cmd:)} {cmd:number(}#{cmd:)} {cmd:nodots} {cmd:generate(}{it:newvars}{cmd:)}]{p_end}

{title:Description}

{pstd}
{cmd:ctset} declares count-time data: for each time, the number of failures, censored and entering subjects. {cmd:ct} shows the declaration. {cmd:cttost} converts count-time data to survival-time data and {cmd:stset}s them; entry times of the subjects are assigned at random when only counts are known. {cmd:sttoct} converts survival-time data to count-time data. {cmd:sttocc} converts survival-time data to case-control data by sampling {cmd:number()} controls for each case from the subjects at risk when the case fails, matching on {cmd:match()}; it creates {cmd:_case}, {cmd:_set} and {cmd:_time}. Tied failure times are split at random, so set the seed for reproducible results.

{title:See also}

{p 4 4 2}{help stset}, {help strate}, {help stsplit}

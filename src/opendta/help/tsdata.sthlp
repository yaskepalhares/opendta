{smcl}
{title:Title}

{p 4 4 2}
{bf:tsappend, tsfill, tsreport, rolling} {hline 2} Time-series data tools

{title:Syntax}

{p 8 16 2}{cmd:tsappend} {cmd:,} {c -(}{opt add(#)} | {opt last(date)} {opt tsfmt(string)}{c )-} [{opt panel(#)}]{p_end}
{p 8 16 2}{cmd:tsfill} [{cmd:,} {opt full}]{p_end}
{p 8 16 2}{cmd:tsreport} [{it:varlist}] [{it:if}] [{it:in}] [{cmd:,} {opt detail} {opt casewise} {opt panel}]{p_end}
{p 8 16 2}{cmd:rolling} [{it:exp_list}] {cmd:,} {opt window(#)} [{opt recursive} {opt rrecursive} {opt step:size(#)} {opt start(time)} {opt end(time)} {opt keep(varname)} {opt clear} {opt saving(filename)} {opt nodots} {opt noisily} {opt reject(exp)} {opt force}]{cmd::} {it:command}{p_end}

{title:Description}

{pstd}
{cmd:tsappend} adds observations at the end of the data (and of each panel with {cmd:panel()}), {cmd:tsfill} inserts the observations that are missing from the time variable, {cmd:tsreport} reports the gaps and the missing values of the variables, and {cmd:rolling} runs a command on windows of the sample (fixed, recursive or reverse recursive) and collects the results into a new dataset. The data must be {cmd:tsset}.
{p_end}

{title:Examples}

{phang2}{cmd:. tsappend, add(4)}{p_end}
{phang2}{cmd:. tsfill}{p_end}
{phang2}{cmd:. rolling _b _se, window(60) clear: regress y x}{p_end}

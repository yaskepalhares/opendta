{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:compare} {hline 2}}Comparing variables and datasets{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:compare} {it:varname1} {it:varname2} [{it:if}] [{it:in}]{p_end}
{p 8 16 2}{cmd:cf} {it:varlist} {cmd:using} {it:filename} [{cmd:,} {opt all} {opt v:erbose}]{p_end}
{p 8 16 2}{cmd:datasignature} [{cmd:set}|{cmd:confirm}|{cmd:report}|{cmd:clear}] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:snapshot} {cmd:save} [{cmd:,} {opt l:abel("label")}] | {cmd:snapshot restore} {it:#} | {cmd:snapshot list} [{it:#}] | {cmd:snapshot label} {it:#} {it:"label"} | {cmd:snapshot erase} {it:numlist}{p_end}

{title:Description}

{pstd}
{cmd:compare} reports how often one variable is lower than, equal to and greater than another, with the differences. {cmd:cf} compares variables of the data in memory with those of a dataset on disk and exits with r(9) when they differ. {cmd:datasignature} computes a signature of the data (number of observations and variables and checksums of the values) and checks whether the data changed since it was set. {cmd:snapshot} keeps copies of the data in memory during the session and restores them.
{p_end}

{title:Remarks}

{pstd}
The checksums of {cmd:datasignature} are OpenDTA's own: the algorithm Stata uses is not documented, so signatures computed by the two programs differ.
{p_end}

{title:Examples}

{phang2}{cmd:. compare income income_rev}{p_end}
{phang2}{cmd:. cf _all using original, verbose}{p_end}
{phang2}{cmd:. snapshot save, label("before cleaning")}{p_end}

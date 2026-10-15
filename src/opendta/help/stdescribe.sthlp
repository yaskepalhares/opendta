{smcl}
{title:Title}

{p 4 4 2}
{bf:stdescribe} {hline 2} Describe survival-time data

{title:Syntax}

{p 8 16 2}
{cmd:stdescribe} [if] [in] [, {cmd:noshow}]

{title:Description}

{pstd}
Summarizes the {cmd:stset} data: number of subjects and records, entry and exit times, gaps, time at risk and failures, in total and per subject (mean, minimum, median and maximum). Saves r(N_sub), r(N_total), r(N_fail), r(risk) and others.

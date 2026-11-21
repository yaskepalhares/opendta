{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:memory} {hline 2}}Report on the memory used by the data and by the session{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:memory}{p_end}

{title:Description}

{pstd}
{cmd:memory} shows the bytes used by the data in memory, by the variable names, formats and labels,
and by matrices. It has the layout of Stata 14's report. OpenDTA does not preallocate memory,
so the {it:allocated} column equals {it:used}, and the items it does not track (ado-files, Mata,
stored results, overhead) show 0.

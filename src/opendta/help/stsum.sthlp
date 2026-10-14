{smcl}
{title:Title}

{p 4 4 2}
{bf:stsum} {hline 2} Summarize survival-time data

{title:Syntax}

{p 8 16 2}
{cmd:stsum} [if] [in] [, {cmd:by(}varname{cmd:)} {cmd:noshow}]

{title:Description}

{pstd}
Shows the time at risk, the incidence rate (failures divided by time at risk), the number of subjects and the 25th, 50th and 75th percentiles of survival time from the Kaplan-Meier estimate, in total and by groups.

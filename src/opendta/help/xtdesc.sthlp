{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:xtdesc} {hline 2}}Panel data: xtline, xttab and xtdata{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:xtline} {it:varlist} [{it:if}] [{it:in}] [{cmd:,} {opt overlay} {it:twoway_options}]{p_end}
{p 8 16 2}{cmd:xttab} {it:varname} [{it:if}] [{it:in}]{p_end}
{p 8 16 2}{cmd:xtdata} [{it:varlist}] [{it:if}] [{it:in}] {cmd:,} {opt fe}|{opt be} [{opt clear}]{p_end}

{title:Description}

{pstd}
{cmd:xtline} plots each panel's time series, one graph per panel or all on one graph with {cmd:overlay}. {cmd:xttab} tabulates a variable showing overall, between (panels) and within (percent of time) frequencies. {cmd:xtdata} replaces the data in memory by their fixed-effects ({cmd:fe}: deviations from the panel means plus the overall mean) or between ({cmd:be}: panel means) transformation; {cmd:re} and {cmd:fa} are not available yet.
{p_end}

{title:Options}

{phang}
{opt overlay} plot all panels on a single graph.
{p_end}

{phang}
{opt clear} allow {cmd:xtdata} to change data that have not been saved.
{p_end}

{title:Examples}

{phang2}{cmd:. xtset id year}{p_end}
{phang2}{cmd:. xtline wage, overlay}{p_end}
{phang2}{cmd:. xttab union}{p_end}
{phang2}{cmd:. xtdata wage age, fe clear}{p_end}

{title:Also see}

{pstd}
See also {help xtsum}, {help xtdescribe} and {help xtreg}.
{p_end}

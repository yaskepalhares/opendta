{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:bcal} {hline 2}}Business calendars{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:bcal check} [{it:varlist}] [{cmd:,} {opt rc0}]{p_end}
{p 8 16 2}{cmd:bcal describe} {it:calname}{p_end}
{p 8 16 2}{cmd:bcal load} {it:calname}{p_end}
{p 8 16 2}{cmd:bcal create} {it:filename} [{it:if}] [{it:in}] {cmd:,} {opt from(varname)} [{opt gen:erate(newvar)} {opt purp:ose(text)} {opt datef:ormat(ymd|ydm|myd|mdy|dym|dmy)} {opt range(fromdate todate)} {opt center:date(date)} {opt excl:udemissing(varlist)} {opt maxg:ap(#)} {opt replace} {opt pers:onal}]{p_end}

{title:Description}

{pstd}
A business calendar lists the days that count (for example, trading days), so that dates can be stored with the {cmd:%tb}{it:calname} format, where consecutive business days differ by 1. Calendars are {cmd:.stbcal} files in the working folder or the adopath. {cmd:bcal create} writes a calendar from the dates in a variable; {cmd:bcal describe} and {cmd:bcal load} read one; {cmd:bcal check} lists the {cmd:%tb} variables and their calendars. The functions {cmd:bofd(}{it:calname}{cmd:,} {it:td}{cmd:)} and {cmd:dofb(}{it:bd}{cmd:,} {it:calname}{cmd:)} convert between {cmd:%td} and business dates.
{p_end}

{title:Examples}

{phang2}{cmd:. bcal create trading, from(date) generate(bdate) replace}{p_end}
{phang2}{cmd:. format bdate %tbtrading}{p_end}
{phang2}{cmd:. generate d = dofb(bdate, "trading")}{p_end}

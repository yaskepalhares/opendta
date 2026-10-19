{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:putmata} {hline 2}}Copying variables into Mata{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:putmata} {it:putlist} [{it:if}] [{it:in}] [{cmd:,} {opt omit:missing} {opt v:iew} {opt replace}]{p_end}

{title:Description}

{pstd}
{cmd:putmata} copies variables into Mata: {it:varname} creates a Mata vector of the same name, {it:name}{cmd:=}{it:varname} one with another name, and {it:name}{cmd:=(}{it:varlist}{cmd:)} a matrix with one column per variable. {opt omitmissing} leaves out the observations with a missing value in any of them.
{p_end}

{title:Examples}

{phang2}{cmd:. putmata y X=(age weight)}{p_end}
{phang2}{cmd:. mata: invsym(X'X)}{p_end}

{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:display} {hline 2}}Display strings and values of expressions{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:display} [{it:display_directive} [...]]{p_end}

{title:Description}

{pstd}
Evaluates expressions and shows the results. Directives: {cmd:"}{it:text}{cmd:"}, {it:exp}, {it:%fmt} {it:exp}, {cmd:_skip(#)}, {cmd:_column(#)}, {cmd:_newline}[{cmd:(#)}], {cmd:_continue}, {cmd:_dup(#)}, {cmd:as text}|{cmd:result}|{cmd:error}|{cmd:input}, {cmd:_request(}{it:macname}{cmd:)} (reads a line typed by the user into the global {it:macname}, or into a local when the name starts with {cmd:_}).
{p_end}

{title:Examples}

{phang2}{cmd:. display 2 + 2}{p_end}
{phang2}{cmd:. display "mean = " %6.2f r(mean)}{p_end}
{phang2}{cmd:. display _n(2) as result c(N)}{p_end}

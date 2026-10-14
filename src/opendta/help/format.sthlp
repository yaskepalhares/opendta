{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:format} {hline 2}}Display formats, storage types and compress{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:format} {it:varlist} {it:%fmt}  |  {cmd:format} [{it:varlist}]{p_end}
{p 8 16 2}{cmd:recast} {it:type} {it:varlist} [{cmd:,} {opt force}]{p_end}
{p 8 16 2}{cmd:compress} [{it:varlist}]{p_end}

{title:Description}

{pstd}
{cmd:format} sets how values are displayed ({cmd:%9.2f}, {cmd:%10.0g}, {cmd:%td}, {cmd:%-20s}, {cmd:%9.0fc}...). {cmd:recast} changes the storage type; {cmd:compress} chooses the smallest type that holds the values exactly.
{p_end}

{title:Examples}

{phang2}{cmd:. format income %12.2fc}{p_end}
{phang2}{cmd:. format date %td}{p_end}
{phang2}{cmd:. compress}{p_end}

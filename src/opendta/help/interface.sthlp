{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:browse} {hline 2}}Windows of the interface{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:browse} [{it:varlist}] [{it:if}] [{it:in}] [{cmd:,} {opt nol:abel}]{p_end}
{p 8 16 2}{cmd:edit} [{it:varlist}] [{it:if}] [{it:in}]{p_end}
{p 8 16 2}{cmd:varmanage}  |  {cmd:doedit} [{it:filename}]{p_end}
{p 8 16 2}{cmd:help} [{it:topic}]  |  {cmd:search} {it:word}  |  {cmd:view} {it:filename}{p_end}
{p 8 16 2}{cmd:about}  |  {cmd:creturn list}  |  {cmd:clear} [{cmd:all}|{cmd:results}|{cmd:programs}|{cmd:mata}]{p_end}

{title:Description}

{pstd}
{cmd:browse} and {cmd:edit} open the Data Editor, read-only or editable (each edit becomes a command in Results and Review). {cmd:varmanage} opens the Variables Manager; {cmd:doedit} the Do-file Editor; {cmd:help} and {cmd:view} the Viewer.
{p_end}

{pstd}
The menus Statistics and Graphics open dialogs that build the command, show it before running, and offer OK, Submit and Copy to Command window. The Graph window has a simple Graph Editor (Ctrl+E) for titles, legend, scheme, size and log scales.
{p_end}

{title:Examples}

{phang2}{cmd:. browse if age > 60}{p_end}
{phang2}{cmd:. edit id age}{p_end}
{phang2}{cmd:. help regress}{p_end}

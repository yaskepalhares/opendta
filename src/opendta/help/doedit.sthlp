{smcl}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:doedit} {hline 2}}Open the Do-file Editor{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:doedit} [{it:filename}]{p_end}

{title:Description}

{pstd}
{cmd:doedit} opens the Do-file Editor, and the file if one is given (a
missing {cmd:.do} suffix is added). Without the interface (batch mode) it does
nothing.
{p_end}

{title:The editor}

{p2colset 9 30 32 2}{...}
{p2col :{bf:Do} (Ctrl+D)}execute the selection, or the whole file, with {cmd:do}{p_end}
{p2col :{bf:Run} (Ctrl+R)}the same with {cmd:run}, without echoing commands{p_end}
{p2col :{bf:Find} (Ctrl+F)}find and replace{p_end}
{p2colreset}{...}

{pstd}
A selection is extended to whole lines. As in Stata, the code is written to
a temporary file and executed with {cmd:do} (or {cmd:run}), so local macros and
{cmd:#delimit} behave as in a do-file. Commands, strings, comments and macros
are colored; each open file has its own tab, and unsaved changes are marked
with {cmd:*}.
{p_end}

{pstd}
The Graphics, Statistics and Data menus write a command template in the
Command window, ready to complete with variable names.
{p_end}

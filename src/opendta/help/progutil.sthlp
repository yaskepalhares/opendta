{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:pause} {hline 2}}Interrupt a program, wait, load a plugin, show a source{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:pause} {c -(}{cmd:on} | {cmd:off} | [{it:message}]{c )-}{break}
{cmd:sleep} {it:#}{break}
{cmd:viewsource} {it:filename}{break}
{cmd:program} {it:handle}{cmd:, plugin} [{cmd:using(}{it:filespec}{cmd:)}]{break}
{cmd:window} {c -(}{cmd:fopen} | {cmd:fsave}{c )-} {it:macroname} {cmd:"}{it:title}{cmd:"}
{cmd:"}{it:filter}{cmd:"} [{it:extension}]{break}
{cmd:window stopbox} {c -(}{cmd:stop} | {cmd:note} | {cmd:rusure}{c )-} {cmd:"}{it:line 1}{cmd:"} ...{break}
{cmd:window} {c -(}{cmd:push} | {cmd:manage} | {cmd:menu}{c )-} ...

{title:Description}

{pstd}
With {cmd:pause on}, a {cmd:pause} {it:message} inside a program shows the message and
suspends the program: you get the prompt {cmd:-> .} and may type any command (the local macros
of the program are visible). {cmd:end} or {cmd:q} resumes the program; {cmd:BREAK} stops it
with {cmd:r(1)}. With {cmd:pause off}, the default, {cmd:pause} does nothing.
{p_end}

{pstd}
{cmd:sleep} waits {it:#} milliseconds.
{p_end}

{pstd}
{cmd:viewsource} shows a file found on the adopath (give the extension, as in
{cmd:viewsource kappa.ado}). For a command that OpenDTA has built in,
{cmd:viewsource} {it:name}{cmd:.ado} shows OpenDTA's own source of that command, not Stata's. If there is nothing to show, the
error is {cmd:r(601)}.
{p_end}

{pstd}
{cmd:program} {it:handle}{cmd:, plugin} names a compiled Stata plugin. OpenDTA cannot load
it: it finds the file (otherwise {cmd:r(601)}) and then says so. {cmd:plugin call} does the
same.
{p_end}

{pstd}
{cmd:window fopen} and {cmd:fsave} put the file chosen in a dialog in the global
macro {it:macroname}; if the dialog is cancelled the macro is empty and the return code is
{cmd:601}. {cmd:window stopbox} shows a message box: {cmd:stop} ends with {cmd:r(1)}, {cmd:note} with 0,
{cmd:rusure} with 0 for Yes and 1 for No. {cmd:window push}, {cmd:manage} and {cmd:menu} act on the
windows and menus of the graphical interface. Without it (batch mode), {cmd:fopen}, {cmd:fsave} and
{cmd:rusure} behave as a cancelled dialog and {cmd:push}, {cmd:manage} and {cmd:menu} do nothing.
{p_end}

{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:cd} {hline 2}}Working directory and files{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:cd} [{cmd:"}{it:directory}{cmd:"}]  |  {cmd:pwd}{p_end}
{p 8 16 2}{cmd:dir} [{it:filespec}]  |  {cmd:ls}  |  {cmd:mkdir} {it:directory}{p_end}
{p 8 16 2}{cmd:erase} {it:filename}  |  {cmd:rm} {it:filename}{p_end}
{p 8 16 2}{cmd:type} {it:filename}  |  {cmd:copy} {it:from} {it:to} [{cmd:,} {opt replace}]{p_end}
{p 8 16 2}{cmd:which} {it:command}  |  {cmd:findfile} {it:filename}  |  {cmd:sysdir}  |  {cmd:adopath}  |  {cmd:discard}{p_end}

{title:Description}

{pstd}
{cmd:cd} changes the working directory (the status bar shows it). {cmd:which} tells whether a command is built in or an ado-file; {cmd:findfile}, {cmd:sysdir} and {cmd:adopath} show where ado-files are searched. OpenDTA never runs Stata's own ado-files.
{p_end}

{title:Examples}

{phang2}{cmd:. cd "~/projects/survey"}{p_end}
{phang2}{cmd:. pwd}{p_end}
{phang2}{cmd:. which regress}{p_end}

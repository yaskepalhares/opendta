{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:do} {hline 2}}Run do-files and control execution{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:do} {it:filename} [{it:arguments}] [{cmd:,} {opt nostop}]{p_end}
{p 8 16 2}{cmd:run} {it:filename} [{it:arguments}]  |  {cmd:include} {it:filename}{p_end}
{p 8 16 2}{cmd:capture} {it:command}  |  {cmd:quietly} {it:command}  |  {cmd:noisily} {it:command}{p_end}
{p 8 16 2}{cmd:exit} [{it:#}] [{cmd:,} {opt clear}]  |  {cmd:continue} [{cmd:,} {opt break}]{p_end}
{p 8 16 2}{cmd:error} {it:#}  |  {cmd:confirm} {cmd:variable}|{cmd:file}|{cmd:new variable}|{cmd:number}|... {it:args}{p_end}
{p 8 16 2}{cmd:version} {it:#} [{cmd::} {it:command}]{p_end}
{p 8 16 2}{cmd:timer} {cmd:on}|{cmd:off}|{cmd:list}|{cmd:clear} [{it:#}]{p_end}

{title:Description}

{pstd}
{cmd:do} runs a do-file showing the commands; {cmd:run} runs it silently; {cmd:include} runs it in the current context (local macros are shared). The Do-file Editor (Ctrl+9) runs the selection or the whole file. Loops and conditions are {cmd:foreach}, {cmd:forvalues}, {cmd:while} and {cmd:if}/{cmd:else}.
{p_end}

{title:Examples}

{phang2}{cmd:. do analysis.do}{p_end}
{phang2}{cmd:. capture confirm variable bmi}{p_end}
{phang2}{cmd:. if _rc display "no bmi"}{p_end}

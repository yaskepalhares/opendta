{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:local} {hline 2}}Macros, scalars and tokenize{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:local} {it:name} [{cmd:=}{it:exp} | {cmd:"}{it:text}{cmd:"} | {cmd::}{it:extended_fcn}]{p_end}
{p 8 16 2}{cmd:global} {it:name} [{cmd:=}{it:exp} | {cmd:"}{it:text}{cmd:"} | {cmd::}{it:extended_fcn}]{p_end}
{p 8 16 2}{cmd:macro list}  |  {cmd:macro drop} {it:names}{p_end}
{p 8 16 2}{cmd:scalar} [{cmd:define}] {it:name} {cmd:=} {it:exp}  |  {cmd:scalar list}  |  {cmd:scalar drop}{p_end}
{p 8 16 2}{cmd:tokenize} [{cmd:"}]{it:string}[{cmd:"}] [{cmd:,} {opt p:arse(pchars)}]{p_end}

{title:Description}

{pstd}
Local macros ({cmd:`name'}) live in the do-file or program that created them; global macros ({cmd:$name}) in the whole session. Extended functions include {cmd:word # of}, {cmd:list}, {cmd:subinstr}, {cmd:variable label}, {cmd:value label}, {cmd:dir}, {cmd:piece} and others. {cmd:tokenize} splits a string into {cmd:`1'}, {cmd:`2'}, ...
{p_end}

{title:Examples}

{phang2}{cmd:. local vars age bmi}{p_end}
{phang2}{cmd:. foreach v of local vars {}{p_end}
{phang2}{cmd:.     summarize `v'}{p_end}
{phang2}{cmd:. }}{p_end}
{phang2}{cmd:. scalar s = r(mean)}{p_end}
{phang2}{cmd:. tokenize "a b c"}{p_end}

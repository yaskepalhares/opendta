{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:set} {hline 2}}Settings{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:set} {it:setting} {it:value} [{cmd:,} {opt perm:anently}]{p_end}

{title:Description}

{pstd}
Settings accepted by OpenDTA include {cmd:obs}, {cmd:seed}, {cmd:rng}, {cmd:more}, {cmd:linesize}, {cmd:dp}, {cmd:type}, {cmd:varabbrev}, {cmd:trace}, {cmd:level}, {cmd:scheme}, {cmd:graphics} and the OpenDTA additions {help hints}, {help superscript}, {help numerics} and {help language}. {cmd:creturn list} shows the current values.
{p_end}

{title:Examples}

{phang2}{cmd:. set obs 100}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. set linesize 120}{p_end}
{phang2}{cmd:. set language pt}{p_end}

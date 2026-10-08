{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:preserve} {hline 2}}Keep a copy of the data and bring it back{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:preserve}

{p 8 16 2}
{cmd:restore} [{cmd:,} {cmd:not} {cmd:preserve}]

{title:Description}

{pstd}
{cmd:preserve} keeps a copy of the data in memory; {cmd:restore} brings it
back. Inside a program or do-file, preserved data come back automatically
when it ends, even after an error. {cmd:restore, not} cancels the copy;
{cmd:restore, preserve} brings the data back and keeps the copy.
{p_end}

{pstd}
Temporary names live in the same way: {cmd:tempvar}, {cmd:tempname} and
{cmd:tempfile} create names that disappear when the program ends.
{p_end}

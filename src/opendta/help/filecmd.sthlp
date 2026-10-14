{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:file} {hline 2}}Read and write text and binary files{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:file open} {it:handle} {cmd:using} {it:filename} {cmd:,} {cmd:read}|{cmd:write}|{cmd:append} [{opt text}|{opt binary}] [{opt replace}]{p_end}
{p 8 16 2}{cmd:file write} {it:handle} {it:directives}{p_end}
{p 8 16 2}{cmd:file read} {it:handle} {it:localmacro}{p_end}
{p 8 16 2}{cmd:file close} {it:handle}|{cmd:_all}  |  {cmd:file query}{p_end}

{title:Description}

{pstd}
Low-level file input and output for programs: {cmd:file write} takes the same directives as {help display}; after {cmd:file read}, {cmd:r(eof)} is 1 at the end of the file.
{p_end}

{title:Examples}

{phang2}{cmd:. tempname fh}{p_end}
{phang2}{cmd:. file open `fh' using "out.txt", write replace}{p_end}
{phang2}{cmd:. file write `fh' "n = " (_N) _n}{p_end}
{phang2}{cmd:. file close `fh'}{p_end}

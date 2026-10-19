{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:file} {hline 2}}Read and write text and binary files{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:file open} {it:handle} {cmd:using} {it:filename} {cmd:,} {cmd:read}|{cmd:write}|{cmd:read write} [{opt append}] [{opt text}|{opt binary}] [{opt replace}]{p_end}
{p 8 16 2}{cmd:file write} {it:handle} {it:directives}{p_end}
{p 8 16 2}{cmd:file read} {it:handle} {it:localmacro}{p_end}
{p 8 16 2}{cmd:file write} {it:handle} {it:%fmt} {cmd:(}{it:exp}{cmd:)} ...  |  {cmd:file read} {it:handle} {it:%fmt} {it:name}  (binary files){p_end}
{p 8 16 2}{cmd:file seek} {it:handle} {cmd:query}|{cmd:tof}|{cmd:eof}|{it:#}  |  {cmd:file set} {it:handle} {cmd:byteorder} {cmd:hilo}|{cmd:lohi}|{cmd:1}|{cmd:2}{p_end}
{p 8 16 2}{cmd:file close} {it:handle}|{cmd:_all}|{cmd:all}  |  {cmd:file query}{p_end}

{title:Description}

{pstd}
Low-level file input and output for programs: {cmd:file write} takes the same directives as {help display}; after {cmd:file read}, {cmd:r(eof)} is 1 at the end of the file.
{p_end}

{pstd}
Binary files use the formats {cmd:%1b}, {cmd:%2b} and {cmd:%4b} (integers with missing values), their signed and unsigned forms
({cmd:%1bs}, {cmd:%1bu}, ...), {cmd:%4z} and {cmd:%8z} (floating point) and {cmd:%}{it:#}{cmd:s} (a string of {it:#} bytes, cut or padded
with zeros). A number read goes to a scalar and a string to a local macro. {cmd:byteorder()} returns 1 (hilo) or 2 (lohi), the order of this
computer; {cmd:file set} changes the order used by a handle. {cmd:file seek} stores the position (the first byte is 0) in {cmd:r(loc)}.
{p_end}

{title:Examples}

{phang2}{cmd:. tempname fh}{p_end}
{phang2}{cmd:. file open `fh' using "out.txt", write replace}{p_end}
{phang2}{cmd:. file write `fh' "n = " (_N) _n}{p_end}
{phang2}{cmd:. file close `fh'}{p_end}

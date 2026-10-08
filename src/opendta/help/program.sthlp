{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:program} {hline 2}}Define programs{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:program} [{cmd:define}] {it:name} [{cmd:,} {cmd:rclass} | {cmd:eclass} | {cmd:sclass}]{break}
{space 4}{it:commands}{break}
{cmd:end}

{p 8 16 2}
{cmd:program drop} {it:name} | {cmd:_all}{space 4}{cmd:program dir}{space 4}{cmd:program list} {it:name}

{title:Description}

{pstd}
A program runs its commands with its own local macros: {cmd:`0'} holds
everything typed after the name and {cmd:`1'}, {cmd:`2'}, ... the words.
Use {help syntax} or {cmd:args} to read the arguments.
{p_end}

{pstd}
An {cmd:rclass} program stores results with {cmd:return scalar},
{cmd:return local} and {cmd:return matrix}; they appear in {cmd:r()} only
if the program ends without error. {cmd:eclass} programs use
{cmd:ereturn}. A program saved as {it:name}{cmd:.ado} in a folder of the
{cmd:adopath} is loaded automatically the first time it is used.
{p_end}

{title:Example}

        {cmd:program media, rclass}
        {cmd:    syntax varname [if] [in]}
        {cmd:    marksample touse}
        {cmd:    ...}
        {cmd:    return scalar mean = `m'}
        {cmd:end}

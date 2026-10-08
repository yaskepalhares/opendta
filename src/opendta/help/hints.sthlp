{smcl}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:set hints} {hline 2}}Explain errors{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:set hints} {c -(}{cmd:on}|{cmd:off}{c )-} [{cmd:,} {cmdab:perm:anently}]

{title:Description}

{pstd}
When an error stops a command, OpenDTA shows the usual error message and
return code {cmd:r(#)}. With {cmd:set hints on} (the default), one or two
lines in between explain what went wrong in the context of the command
and suggest a fix: similar variable names, similar file names in the
folder, the number of observations in memory, unbalanced quotes and so on.
{p_end}

{pstd}
The return codes never change, so do-files that test {cmd:_rc} behave the
same with hints on or off. Errors raised on purpose with {cmd:error #}
are not explained.
{p_end}

{pstd}
{cmd:permanently} also remembers the choice the next time OpenDTA opens.
The same switch is in the menu {bf:Help > Explain Errors}.
{p_end}

{title:Example}

        {cmd:. replace nota = 9 in 2}
        {err:Obs. nos. out of range}
          {txt:-> "in 2" asks for observations that do not exist; the dataset has 0 observations.}
        {err:r(198);}

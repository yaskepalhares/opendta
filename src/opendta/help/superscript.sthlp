{smcl}
{title:Title}

{p2colset 5 28 30 2}{...}
{p2col :{cmd:set superscript} {hline 2}}Superscript exponents in scientific notation{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:set superscript} {c -(}{cmd:on}|{cmd:off}{c )-} [{cmd:,} {cmdab:perm:anently}]

{title:Description}

{pstd}
Stata writes very large and very small numbers as {cmd:1.00000e+10} or
{cmd:2.50e-05}. With {cmd:set superscript on} (the default), OpenDTA writes
them the usual mathematical way, with the power of ten as a superscript:
{cmd:1.0000×10¹⁰} and {cmd:2.5×10⁻⁵}. Each number keeps the width of its
display format, so the mantissa may show one digit less.
{p_end}

{pstd}
Only what is shown on screen changes: {cmd:display}, {cmd:list},
{cmd:summarize}, tables, matrices, Mata output and the Data Editor. Text that
can be read back keeps the Stata form, so do-files work the same with the
setting on or off: {cmd:string()}, {cmd:tostring}, macros
({cmd:local x : display ...}), {cmd:file write}, {cmd:export delimited},
{cmd:outsheet}, and Mata's {cmd:strofreal()} and {cmd:sprintf()}.
Stored data never change.
{p_end}

{pstd}
{cmd:set superscript off} shows the Stata form ({cmd:e+10}), for example to
compare a log with one made by Stata.
{cmd:permanently} also remembers the choice the next time OpenDTA opens.
The same switch is in the menu {bf:Help > Superscript Exponents}.
{p_end}

{title:Example}

        {cmd:. display 1e10/3}
        {res:3.333×10⁹}

        {cmd:. set superscript off}
        {cmd:. display 1e10/3}
        {res:3.333e+09}

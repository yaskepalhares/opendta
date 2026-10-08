{smcl}
{title:Title}

{p2colset 5 16 18 2}{...}
{p2col :{cmd:use} {hline 2}}Load a .dta dataset{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:use} {it:filename} [{cmd:,} {cmd:clear} {cmdab:nol:abel}]

{p 8 16 2}
{cmd:use} [{varlist}] [{it:if}] [{it:in}] {cmd:using} {it:filename} [{cmd:,} {cmd:clear} {cmdab:nol:abel}]

{title:Description}

{pstd}
{cmd:use} replaces the data in memory with a dataset saved in {cmd:.dta}
format. When the file name has no extension, {cmd:.dta} is assumed.
Files written by versions 8 and later of the format are read directly;
older ones need the optional {cmd:pyreadstat} package.
{p_end}

{pstd}
The second form loads only some variables or observations; the {it:if}
condition may use variables that are not loaded.
{p_end}

{title:Options}

{phang}
{cmd:clear} allows the data in memory to be replaced even if they were
changed and not saved.

{phang}
{cmd:nolabel} loads the data without value labels.

{title:Examples}

        {cmd:. use pessoas}
        {cmd:. use nome idade using pessoas if idade >= 18, clear}

{title:Also see}

{pstd}
{help save}, {help import}, {help describe}
{p_end}

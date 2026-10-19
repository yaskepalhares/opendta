{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:describe} {hline 2}}Describe the dataset{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmdab:d:escribe} [{varlist}] [{cmd:,} {cmdab:s:hort} {cmdab:simple} {cmdab:f:ullnames} {cmdab:n:umbers} {cmd:varlist}]

{p 8 16 2}
{cmdab:d:escribe} [{varlist}] {cmd:using} {it:filename} [{cmd:,} {cmdab:s:hort} {cmdab:simple} {cmd:varlist}]

{p 8 16 2}
{cmdab:d:escribe} [{varlist}] {cmd:, replace} [{cmd:clear}]

{title:Description}

{pstd}
{cmd:describe} lists the variables with their storage type, display
format, value label and variable label, preceded by the number of
observations and variables, the dataset label and the date it was saved.
Variables with notes are marked with {cmd:*}. {cmd:numbers} shows the position of each variable,
{cmd:varlist} also stores {cmd:r(varlist)} and {cmd:r(sortlist)}, and {cmd:using} describes a {cmd:.dta} file without
loading it. With {cmd:replace} no report is shown: the data in memory become one observation per variable
({cmd:position}, {cmd:name}, {cmd:type}, {cmd:isnumeric}, {cmd:format}, {cmd:vallab}, {cmd:varlab}); {cmd:clear} lets it
replace data that were not saved.
{p_end}

{title:Also see}

{pstd}
{help list}, {help use}
{p_end}

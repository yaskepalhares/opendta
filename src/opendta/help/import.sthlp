{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:import} {hline 2}}Import and export text and Excel files{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:import delimited} [{cmd:using}] {it:filename} [{cmd:,} {cmd:clear}
{cmdab:delim:iters("}{it:char}{cmd:")} {cmdab:varn:ames(}{it:#}|{cmd:nonames)}
{cmd:case(}{it:preserve}|{it:lower}|{it:upper}{cmd:)} {cmd:asdouble}
{cmdab:rowr:ange(}[{it:start}][{cmd::}{it:end}]{cmd:)} {cmdab:colr:ange()}
{cmdab:stringc:ols(}{it:numlist}{cmd:)} {cmdab:numericc:ols(}{it:numlist}{cmd:)}]

{p 8 16 2}
{cmd:export delimited} [{varlist}] {cmd:using} {it:filename} [{it:if}] [{it:in}]
[{cmd:,} {cmd:replace} {cmdab:delim:iter()} {cmdab:novar:names} {cmdab:nol:abel} {cmd:quote}]

{p 8 16 2}
{cmd:import excel} [{cmd:using}] {it:filename} [{cmd:,} {cmd:clear} {cmd:sheet("}{it:name}{cmd:")}
{cmd:cellrange(}{it:A1:D20}{cmd:)} {cmd:firstrow} {cmd:allstring} {cmd:describe}]

{p 8 16 2}
{cmd:export excel} [{varlist}] {cmd:using} {it:filename} [{it:if}] [{it:in}]
[{cmd:,} {cmd:replace} {cmd:sheet("}{it:name}{cmd:"}[{cmd:, modify}|{cmd:replace}]{cmd:)}
{cmd:cell(}{it:B2}{cmd:)} {cmd:firstrow(}{it:variables}|{it:varlabels}{cmd:)} {cmdab:nol:abel}]

{title:Description}

{pstd}
{cmd:import delimited} reads comma- or tab-separated text, choosing the
smallest storage type that holds each column. {cmd:import excel} reads
{cmd:.xlsx} workbooks (and {cmd:.xls} with the optional {cmd:xlrd}
package); Excel dates become {cmd:%td} or {cmd:%tc} variables.
The older {cmd:insheet} and {cmd:outsheet} are also available.
{p_end}

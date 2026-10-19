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

{p 8 16 2}
{cmd:import sas} [{it:namelist}] [{it:if}] [{it:in}] {cmd:using} {it:filename} [{cmd:,} {cmd:clear} {cmd:case(}{it:preserve}|{it:lower}|{it:upper}{cmd:)}
{cmd:bcat(}{it:filename}{cmd:)} {cmd:encoding("}{it:enc}{cmd:")}]

{p 8 16 2}
{cmd:import spss} [{it:namelist}] [{it:if}] [{it:in}] {cmd:using} {it:filename} [{cmd:,} {cmd:clear} {cmd:case()} {cmd:zsav} {cmd:encoding()}]

{p 8 16 2}
{cmd:import sasxport} {it:filename} [{cmd:,} {cmd:clear} {cmd:novallabels}]

{p 8 16 2}
{cmd:export sasxport} [{it:varlist} {cmd:using}] {it:filename} [{it:if}] [{it:in}] [{cmd:,} {cmd:replace} {cmd:vallabfile(none)}]

{title:Description}

{pstd}
{cmd:import delimited} reads comma- or tab-separated text, choosing the
smallest storage type that holds each column. {cmd:import excel} reads
{cmd:.xlsx} workbooks (and {cmd:.xls} with the optional {cmd:xlrd}
package; {cmd:export excel} writes {cmd:.xls} with the optional {cmd:xlwt}
package); Excel dates become {cmd:%td} or {cmd:%tc} variables.
The older {cmd:insheet} and {cmd:outsheet} are also available.
{p_end}

{pstd}
{cmd:import sas}, {cmd:import spss} and {cmd:import sasxport} read through the optional {cmd:pyreadstat} package and show
{cmd:(}{it:k}{cmd: vars,} {it:n}{cmd: obs)}; {cmd:r(N)} and {cmd:r(k)} hold the counts. Names that Stata does not accept become
{cmd:v}{it:#} (the original name is kept in the characteristic {cmd:original_name}), whole numbers use the smallest integer type
that holds them, and dates become {cmd:%td} or {cmd:%tc}. SPSS value labels of numeric variables are named {cmd:labels}{it:#}; SAS
value labels come from {cmd:bcat()}. {cmd:export sasxport} writes version 5 files, whose variable names have up to 8 characters.
Value labels in {cmd:formats.xpf}, {cmd:import sasxport, describe} and {cmd:import haver} are not available.
{p_end}

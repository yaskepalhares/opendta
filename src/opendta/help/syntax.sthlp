{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:syntax} {hline 2}}Parse the arguments of a program{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:syntax} [{it:varlist}|{it:namelist}|{it:anything}] [{cmd:if}] [{cmd:in}]
[{cmd:using}] [{cmd:=}{it:exp}] [{it:weight}] [{cmd:,} {it:options}]

{title:Description}

{pstd}
{cmd:syntax} checks {cmd:`0'} against the description and fills the local
macros {cmd:varlist}, {cmd:if}, {cmd:in}, {cmd:using}, {cmd:exp},
{cmd:weight} and one macro per option. Parts in brackets are optional.
Capital letters in an option name give its shortest abbreviation:
{cmd:Detail} accepts {cmd:d}, {cmd:de}, ..., {cmd:detail}.
{p_end}

{pstd}
Option arguments: {cmd:Level(integer 95)}, {cmd:Title(string)},
{cmd:BY(varname)}, {cmd:Values(numlist)}, {cmd:noLOG} (macro {cmd:log}
holds {cmd:nolog}) and {cmd:*} for any other options, collected in
{cmd:`options'}. {cmd:gettoken}, {cmd:tokenize}, {cmd:marksample} and
{cmd:markout} complete the toolkit.
{p_end}

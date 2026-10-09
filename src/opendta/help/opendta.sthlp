{smcl}
{title:OpenDTA help}

{pstd}
OpenDTA runs analysis scripts (do-files), reads and writes {cmd:.dta}
datasets and keeps every change as a command. Type {cmd:help} followed by
a command name, or click a topic below.
{p_end}

{title:Data files}

{p2colset 5 24 26 2}{...}
{p2col :{help use}}load a .dta dataset{p_end}
{p2col :{help save}}save the data in memory (also {cmd:saveold}){p_end}
{p2col :{help import}}import delimited text and Excel files{p_end}
{p2col :{help infile}}read free-format and fixed-column text ({cmd:infile}, {cmd:infix}){p_end}
{p2colreset}{...}

{title:Working with data}

{p2colset 5 24 26 2}{...}
{p2col :{help describe}}describe the dataset{p_end}
{p2col :{help list}}list values{p_end}
{p2col :{help generate}}create and change variables ({cmd:generate}, {cmd:replace}){p_end}
{p2col :{help preserve}}save a copy of the data and bring it back{p_end}
{p2col :{help assert}}check that a condition holds{p_end}
{p2colreset}{...}

{title:Statistics and data management}

{p2colset 5 24 26 2}{...}
{p2col :{help summarize}}summary statistics ({cmd:tabstat}, {cmd:centile}, {cmd:correlate}...){p_end}
{p2col :{help ttest}}tests for means and proportions ({cmd:prtest}, {cmd:ci}){p_end}
{p2col :{help tabulate}}frequency tables ({cmd:tab1}, {cmd:tab2}, {cmd:table}){p_end}
{p2col :{help egen}}extensions to generate{p_end}
{p2col :{help collapse}}group statistics ({cmd:contract}, {cmd:expand}, {cmd:fillin}){p_end}
{p2col :{help merge}}combine datasets ({cmd:append}, {cmd:joinby}, {cmd:cross}){p_end}
{p2col :{help reshape}}wide and long forms{p_end}
{p2col :{help recode}}recode and convert ({cmd:encode}, {cmd:destring}, {cmd:split}...){p_end}
{p2col :{help duplicates}}duplicates, {cmd:isid}, {cmd:levelsof}{p_end}
{p2col :{help set_seed:set seed}}random numbers{p_end}
{p2colreset}{...}

{title:Programming}

{p2colset 5 24 26 2}{...}
{p2col :{help program}}define programs; {cmd:return}, {cmd:ereturn}{p_end}
{p2col :{help syntax}}parse a program's arguments{p_end}
{p2col :{help matrix}}matrices and matrix expressions{p_end}
{p2col :{help log}}record the session in a file{p_end}
{p2col :{help mata}}the matrix programming language{p_end}
{p2colreset}{...}

{title:OpenDTA additions}

{p2colset 5 24 26 2}{...}
{p2col :{help hints}}error explanations ({cmd:set hints}){p_end}
{p2col :{help superscript}}superscript exponents, 1.0×10¹⁰ ({cmd:set superscript}){p_end}
{p2colreset}{...}

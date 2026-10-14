{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:label} {hline 2}}Labels for datasets, variables and values{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:label data} [{cmd:"}{it:text}{cmd:"}]{p_end}
{p 8 16 2}{cmd:label variable} {it:varname} [{cmd:"}{it:text}{cmd:"}]{p_end}
{p 8 16 2}{cmd:label define} {it:lblname} {it:#} {cmd:"}{it:text}{cmd:"} [...] [{cmd:,} {opt add} {opt modify} {opt replace}]{p_end}
{p 8 16 2}{cmd:label values} {it:varlist} [{it:lblname}|{cmd:.}]{p_end}
{p 8 16 2}{cmd:label list} [{it:lblnames}]  |  {cmd:label dir}  |  {cmd:label drop} {it:lblnames}|{cmd:_all}{p_end}
{p 8 16 2}{cmd:notes} [{it:evarname}]{cmd::} {it:text}  |  {cmd:notes} [{cmd:list}]  |  {cmd:notes drop} ...{p_end}
{p 8 16 2}{cmd:char} {it:evarname}{cmd:[}{it:charname}{cmd:]} [{it:text}]  |  {cmd:char list}{p_end}

{title:Description}

{pstd}
Labels describe the data and are saved with the {cmd:.dta} file. Value labels attach text to numeric codes and are shown by {cmd:list}, {cmd:tabulate} and the Data Editor. {cmd:notes} and {cmd:char} keep notes and characteristics.
{p_end}

{title:Examples}

{phang2}{cmd:. label variable bmi "Body mass index"}{p_end}
{phang2}{cmd:. label define sexlbl 1 "Male" 2 "Female"}{p_end}
{phang2}{cmd:. label values sex sexlbl}{p_end}
{phang2}{cmd:. label list sexlbl}{p_end}

{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:drop} {hline 2}}Drop, keep, sort, rename and order{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:drop} {it:varlist}  |  {cmd:drop} {cmd:if} {it:exp}  |  {cmd:drop} {cmd:in} {it:range}{p_end}
{p 8 16 2}{cmd:keep} {it:varlist}  |  {cmd:keep} {cmd:if} {it:exp}  |  {cmd:keep} {cmd:in} {it:range}{p_end}
{p 8 16 2}{cmd:sort} {it:varlist} [{it:in}] [{cmd:,} {opt stable}]{p_end}
{p 8 16 2}{cmd:gsort} [{cmd:+}|{cmd:-}]{it:varname} [...] [{cmd:,} {opt gen:erate(newvar)} {opt m:first}]{p_end}
{p 8 16 2}{cmd:rename} {it:old} {it:new}  |  {cmd:rename} {cmd:(}{it:old1 old2}{cmd:)} {cmd:(}{it:new1 new2}{cmd:)}{p_end}
{p 8 16 2}{cmd:order} {it:varlist} [{cmd:,} {opt first} {opt last} {opt before(varname)} {opt after(varname)} {opt alpha:betic}]{p_end}
{p 8 16 2}{cmd:count} [{it:if}] [{it:in}]{p_end}

{title:Description}

{pstd}
Basic data management. {cmd:drop} and {cmd:keep} remove variables or observations; {cmd:sort} sorts in ascending order (missing values last; with {it:in}, only the observations in the range) and {cmd:gsort} also in descending order; {cmd:rename} and {cmd:order} change names and the order of the variables; {cmd:count} counts observations.
{p_end}

{title:Examples}

{phang2}{cmd:. drop if age < 18}{p_end}
{phang2}{cmd:. keep id age sex}{p_end}
{phang2}{cmd:. gsort -income}{p_end}
{phang2}{cmd:. rename (v1 v2) (age sex)}{p_end}
{phang2}{cmd:. order id, first}{p_end}
{phang2}{cmd:. count if missing(bmi)}{p_end}

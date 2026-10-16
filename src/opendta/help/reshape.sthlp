{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:reshape} {hline 2}}Convert data between wide and long forms{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:reshape long} {it:stubnames} {cmd:,} {opt i(varlist)} [{opt j(varname [values])} {opt string}]{p_end}
{p 8 16 2}{cmd:reshape wide} {it:stubnames} {cmd:,} {opt i(varlist)} {opt j(varname [values])} [{opt string}]{p_end}
{p 8 16 2}{cmd:reshape long}  |  {cmd:reshape wide}{p_end}
{p 8 16 2}{cmd:reshape i} {it:varlist}  |  {cmd:reshape j} {it:varname} [{it:values}] [{cmd:,} {opt string}]  |  {cmd:reshape xij} {it:stubnames}  |  {cmd:reshape xi} [{it:varlist}]{p_end}
{p 8 16 2}{cmd:reshape} [{cmd:query}]  |  {cmd:reshape clear}  |  {cmd:reshape error}{p_end}

{title:Description}

{pstd}
In wide form each {it:i} has one observation and variables such as
{cmd:inc80}, {cmd:inc81}; in long form each pair ({it:i}, {it:j}) has one
observation with {cmd:inc} and {cmd:year}. A {cmd:@} in a stub marks where
the {it:j} value goes ({cmd:inc@r} matches {cmd:inc80r}). Without stubs,
{cmd:reshape} repeats the last specification in the other direction.
{p_end}

{pstd}
The specification can also be given piece by piece ({cmd:reshape i},
{cmd:reshape j}, {cmd:reshape xij}) and then used with {cmd:reshape long}
or {cmd:reshape wide}; {cmd:reshape} or {cmd:reshape query} shows it and
{cmd:reshape clear} erases it. After a reshape that failed because {it:i}
(or {it:i} and {it:j}) do not identify the observations, or because a
variable is not constant within {it:i}, {cmd:reshape error} lists the
observations at fault.
{p_end}

{smcl}
{title:Title}

{p 4 4 2}
{bf:putexcel} {hline 2} Export results to an Excel file

{title:Syntax}

{p 8 16 2}
{cmd:putexcel} set filename [, sheet(name) replace modify]

{p 8 16 2}
{cmd:putexcel} cell = "text" | (exp) | matrix(name[, names rownames colnames]) ...

{title:Description}

{pstd}
Sets the target .xlsx file.

{pstd}
Writes text, the value of an expression or a matrix starting at {it:cell}. The Stata 14.0 form {cmd:putexcel} {it:cell=...} {cmd:using} {it:filename} is also accepted.

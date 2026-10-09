{smcl}
{title:Title}

{p2colset 5 23 25 2}{...}
{p2col :{cmd:duplicates} {hline 2}}Report, tag or drop duplicate observations{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:duplicates} {c -(}{cmd:report}|{cmd:examples}|{cmd:list}{c )-} [{it:varlist}] [{it:if}] [{it:in}]{p_end}
{p 8 16 2}{cmd:duplicates tag} [{it:varlist}] {cmd:,} {opt gen:erate(newvar)}{p_end}
{p 8 16 2}{cmd:duplicates drop} [{it:varlist}] [{cmd:,} {opt force}]{p_end}
{p 8 16 2}{cmd:isid} {it:varlist} [{cmd:,} {opt sort} {opt missok}]{p_end}
{p 8 16 2}{cmd:levelsof} {it:varname} [{cmd:,} {opt c:lean} {opt l:ocal(name)} {opt m:issing} {opt s:eparate(str)}]{p_end}

{title:Description}

{pstd}
Observations are duplicates when they agree on {it:varlist} (all variables
by default). {cmd:duplicates drop} keeps the first of each set; with a
varlist it needs {opt force}. {cmd:isid} checks that {it:varlist}
identifies the observations; {cmd:levelsof} lists the distinct values.
{p_end}

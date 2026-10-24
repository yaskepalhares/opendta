{smcl}
{title:Title}

{p2colset 5 16 18 2}{...}
{p2col :{cmd:list} {hline 2}}List values of variables{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmdab:l:ist} [{varlist}] [{it:if}] [{it:in}] [{cmd:,} {cmdab:nol:abel} {cmdab:noo:bs} {cmd:clean}
{cmdab:sep:arator(}{it:#}{cmd:)} {cmdab:sepby(}{it:varlist}{cmd:)} {cmdab:div:ider} {cmdab:ab:breviate(}{it:#}{cmd:)} {cmdab:noh:eader} {cmd:header(}{it:#}{cmd:)} {cmd:nodotz}
{cmd:mean}[{cmd:(}{it:varlist}{cmd:)}] {cmd:sum}[{cmd:(}{it:varlist}{cmd:)}] {cmd:N}[{cmd:(}{it:varlist}{cmd:)}] {cmdab:labvar(}{it:varname}{cmd:)} {cmd:nocompress} {cmd:string(}{it:#}{cmd:)} {cmdab:tab:le} {cmdab:disp:lay}]

{title:Description}

{pstd}
{cmd:list} shows the values of the variables in a table, using each
variable's display format. Value labels are shown instead of codes unless
{cmd:nolabel} is given. A separator line is drawn every 5 observations
({cmd:separator(}{it:#}{cmd:)}; 0 for none) or, with {cmd:sepby()}, whenever
the value of one of its variables changes. {cmd:header(}{it:#}{cmd:)} repeats the variable names every {it:#} observations and {cmd:nodotz} shows .z as a blank field. {cmd:mean}, {cmd:sum} and {cmd:N} add a row with the mean, the sum and the number of nonmissing values of the numeric variables (or of the ones named), labeled in the observation-number column or, with {cmd:labvar()}, in the column of that variable. {cmd:nocompress} is accepted and changes nothing. {cmd:string(}{it:#}{cmd:)} cuts string variables to {it:#} characters. {cmd:display} shows each observation in a box of its own, with the variables in groups of a line of names and a line of values; it is also used when the table is wider than the line, unless {cmd:table} is given. For browsing many observations, the Data Editor
({cmd:browse}) is often more convenient.
{p_end}

{title:Examples}

        {cmd:. list}
        {cmd:. list nome idade if idade > 30}
        {cmd:. list in 1/10, nolabel}

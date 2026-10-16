{smcl}
{title:Title}

{p2colset 5 16 18 2}{...}
{p2col :{cmd:list} {hline 2}}List values of variables{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmdab:l:ist} [{varlist}] [{it:if}] [{it:in}] [{cmd:,} {cmdab:nol:abel} {cmdab:noo:bs} {cmd:clean}
{cmdab:sep:arator(}{it:#}{cmd:)} {cmdab:sepby(}{it:varlist}{cmd:)} {cmdab:div:ider} {cmdab:ab:breviate(}{it:#}{cmd:)} {cmdab:noh:eader}]

{title:Description}

{pstd}
{cmd:list} shows the values of the variables in a table, using each
variable's display format. Value labels are shown instead of codes unless
{cmd:nolabel} is given. A separator line is drawn every 5 observations
({cmd:separator(}{it:#}{cmd:)}; 0 for none) or, with {cmd:sepby()}, whenever
the value of one of its variables changes. For browsing many observations, the Data Editor
({cmd:browse}) is often more convenient.
{p_end}

{title:Examples}

        {cmd:. list}
        {cmd:. list nome idade if idade > 30}
        {cmd:. list in 1/10, nolabel}

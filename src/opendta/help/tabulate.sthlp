{smcl}
{title:Title}

{p2colset 5 21 23 2}{...}
{p2col :{cmd:tabulate} {hline 2}}Frequency tables{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:tabulate} {it:varname} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt m:issing} {opt nol:abel} {opt sort} {opt gen:erate(stub)} {opt matcell(name)} {opt matrow(name)}]{p_end}
{p 8 16 2}{cmd:tabulate} {it:var1} {it:var2} [{cmd:,} {opt r:ow} {opt col:umn} {opt cell} {opt e:xpected} {opt nof:req} {opt ch:i2} {opt lr:chi2} {opt V} {opt exact} {opt gam:ma} {opt taub}]{p_end}
{p 8 16 2}{cmd:tab1} {it:varlist}  |  {cmd:tab2} {it:varlist}{p_end}
{p 8 16 2}{cmd:table} {it:rowvar} [{it:colvar}] [{cmd:,} {opt c:ontents(clist)} {opt row} {opt col} {opt f:ormat(%fmt)}]{p_end}

{title:Description}

{pstd}
{cmd:tabulate} with one variable lists the frequency, percent and
cumulative percent of each value. With two variables it builds a
contingency table; the options add percentages, expected frequencies and
tests of independence (Pearson and likelihood-ratio chi-squared,
Cramér's V, Goodman and Kruskal's gamma, Kendall's tau-b and Fisher's
exact test, computed by enumerating every table with the same margins).
{p_end}

{pstd}
{cmd:tab1} makes one-way tables of several variables, {cmd:tab2} all
two-way tables among them. {cmd:table} shows frequencies or other
statistics ({cmd:mean}, {cmd:sd}, {cmd:sum}, {cmd:min}, {cmd:max},
{cmd:median}, {cmd:p#}...) in a compact grid.
{p_end}

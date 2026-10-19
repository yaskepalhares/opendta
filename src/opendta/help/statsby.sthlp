{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:statsby} {hline 2}}Collecting statistics for groups{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:statsby} [{it:exp_list}] [{cmd:,} {opt by(varlist)} {opt clear} {opt sa:ving(filename)} {opt to:tal} {opt sub:sets} {opt nodots} {opt n:oisily} {opt nolegend} {opt v:erbose}] {cmd::} {it:command}{p_end}

{title:Description}

{pstd}
{cmd:statsby} runs {it:command} for each group of {opt by()} and replaces the data with one observation per group holding the results named in {it:exp_list} ({it:name}{cmd:=}{it:exp}, for example {cmd:mean=r(mean)}); {cmd:_b} and {cmd:_se} collect the coefficients and standard errors after estimation commands. {opt total} adds the whole sample; {opt subsets} all combinations of the groups.
{p_end}

{title:Examples}

{phang2}{cmd:. statsby mean=r(mean) sd=r(sd), by(region) clear: summarize income}{p_end}
{phang2}{cmd:. statsby _b, by(sex) clear: regress bmi age}{p_end}

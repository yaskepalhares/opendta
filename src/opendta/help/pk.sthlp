{smcl}
{title:Title}

{p 4 4 2}
{bf:pk} {hline 2} Pharmacokinetic (biopharmaceutical) data

{title:Syntax}

{p 8 16 2}{cmd:pkexamine} {it:time} {it:concentration} [{it:if}] [{it:in}] [{cmd:,} {opt fit(#)} {opt t:rapezoid}]{p_end}
{p 8 16 2}{cmd:pksumm} {it:id} {it:time} {it:concentration} [{it:if}] [{it:in}] [{cmd:,} {opt fit(#)} {opt t:rapezoid} {opt stat(measure)}]{p_end}
{p 8 16 2}{cmd:pkcollapse} {it:time} {it:concentration} [{it:if}]{cmd:,} {opt id(idvar)} [{opt stat(measures)} {opt keep(varlist)} {opt fit(#)} {opt t:rapezoid}]{p_end}
{p 8 16 2}{cmd:pkshape} {it:id} {it:sequence} {it:period1} {it:period2} [...]{cmd:,} {opt o:rder(string)} [{opt out:come(newvar)} {opt tre:atment(newvar)} {opt carry:over(newvar)} {opt seq:uence(newvar)} {opt per:iod(newvar)}]{p_end}
{p 8 16 2}{cmd:pkequiv} {it:outcome} {it:treatment} {it:period} {it:sequence} {it:id} [{it:if}] [{it:in}] [{cmd:,} {opt comp:are(# #)} {opt lim:it(#)} {opt level(#)}]{p_end}
{p 8 16 2}{cmd:pkcross} {it:outcome} [{it:if}] [{it:in}] [{cmd:,} {opt seq:uence(varname)} {opt tre:atment(varname)} {opt per:iod(varname)} {opt id(varname)} {opt carry:over(varname|none)}]{p_end}

{title:Description}

{pstd}
{cmd:pkexamine} computes the maximum concentration and its time, the area under the curve (cubic
spline, or the trapezoidal rule), the elimination rate and half-life from the last {cmd:fit()}
points (3 by default) and the area extrapolated to infinity in three ways. {cmd:pksumm} computes
the measures for each subject and summarizes them; {cmd:pkcollapse} keeps one observation per
subject with the measures. {cmd:pkshape} turns crossover data with one variable per period into
one observation per subject and period, with the treatment and carryover from {cmd:order()}
({cmd:order(ab ba)}). {cmd:pkcross} fits the crossover ANOVA and {cmd:pkequiv} the classic
confidence interval for average bioequivalence.
{p_end}

{pstd}
Not available yet: the graphs, and {cmd:fieller}, {cmd:symmetric}, {cmd:anderson} and {cmd:tost}
in {cmd:pkequiv}, {cmd:model()} in {cmd:pkcross}.
{p_end}

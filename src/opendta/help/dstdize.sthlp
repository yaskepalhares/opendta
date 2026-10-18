{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:dstdize} {hline 2}}Direct and indirect standardization{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:dstdize} {it:charvar} {it:popvar} {it:stratavars} [{it:if}] [{it:in}] {cmd:,} {opt by(groupvars)} [{cmd:using} {it:filename} | {opt base(# | string)}] [{opt l:evel(#)} {opt nores}]{p_end}
{p 8 16 2}{cmd:istdize} {it:casevar} {it:popvar} {it:stratavars} [{it:if}] [{it:in}] {cmd:using} {it:filename} {cmd:,} {opt popvars(casevar popvar)} | {opt rate(ratevar #)} [{opt by(groupvars)} {opt l:evel(#)}]{p_end}

{title:Description}

{pstd}
{cmd:dstdize} computes directly standardized rates: each group's stratum-specific rates of {it:charvar} per {it:popvar} applied to a standard population, which is the populations of all groups combined, the group chosen with {opt base()}, or the one read from a dataset with {cmd:using}. It reports the crude and adjusted rates with a normal confidence interval, and leaves the matrices {cmd:r(Nobs)}, {cmd:r(crude)}, {cmd:r(adj)}, {cmd:r(lb)}, {cmd:r(ub)} and {cmd:r(se)}.
{p_end}

{pstd}
{cmd:istdize} computes indirectly standardized rates: the standard stratum rates (from the dataset named with {cmd:using}, as cases and population with {opt popvars()} or as rates with {opt rate()}) applied to each group's population give the expected cases; the standardized mortality ratio (SMR) is observed over expected, with an exact Poisson interval, and the adjusted rate is the standard crude rate times the SMR.
{p_end}

{title:Examples}

{phang2}{cmd:. dstdize deaths pop agecat, by(state)}{p_end}
{phang2}{cmd:. istdize deaths pop agecat using stdpop, popvars(stddeaths stdpop) by(state)}{p_end}

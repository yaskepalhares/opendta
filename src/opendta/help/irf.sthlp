{smcl}
{title:Title}

{p2colset 5 16 18 2}{...}
{p2col :{cmd:irf} {hline 2}}Impulse-response functions and variance decompositions{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:irf set} {it:filename} [{cmd:,} {opt replace}]{p_end}
{p 8 16 2}{cmd:irf create} {it:irfname} [{cmd:,} {opt set(filename)} {opt replace} {opt step(#)} {opt order(varlist)} {opt nose}]{p_end}
{p 8 16 2}{cmd:irf table} [{it:stats}] [{cmd:,} {opt irf(irfnames)} {opt imp:ulse(varlist)} {opt res:ponse(varlist)} {opt noci} {opt std:error} {opt l:evel(#)}]{p_end}
{p 8 16 2}{cmd:irf graph} [{it:stats}] [{cmd:,} {opt irf(irfnames)} {opt imp:ulse(varlist)} {opt res:ponse(varlist)} {opt noci} {opt l:evel(#)} {it:graph options}]{p_end}
{p 8 16 2}{cmd:irf describe}  |  {cmd:irf drop} {it:irfnames}  |  {cmd:irf rename} {it:old} {it:new}{p_end}

{p 4 4 2}{it:stats}: {cmd:irf}, {cmd:oirf}, {cmd:cirf}, {cmd:coirf}, {cmd:fevd}{p_end}

{title:Description}

{pstd}
After {help var}, {cmd:irf create} computes the impulse-response functions up to
{cmd:step()} periods (default 8): simple ({cmd:irf}), orthogonalized with the Cholesky factor
of {cmd:e(Sigma)} in the order of the variables or of {cmd:order()} ({cmd:oirf}), their
cumulative sums ({cmd:cirf}, {cmd:coirf}) and the forecast-error variance decomposition
({cmd:fevd}), with asymptotic standard errors. The results are saved in the active
{cmd:.irf} file, chosen with {cmd:irf set} or {cmd:set()}; one file holds several sets of
results, each with its {it:irfname}.
{p_end}

{pstd}
{cmd:irf table} shows the statistics with confidence intervals, one column per
statistic and pair impulse/response; {cmd:irf graph} draws them, one panel per pair.
{p_end}

{title:Examples}

{phang2}{cmd:. var dln_inv dln_inc dln_consump, lags(1/2)}{p_end}
{phang2}{cmd:. irf create order1, set(myirf) step(10)}{p_end}
{phang2}{cmd:. irf table oirf fevd, impulse(dln_inc) response(dln_consump)}{p_end}
{phang2}{cmd:. irf graph oirf, impulse(dln_inc)}{p_end}

{title:OpenDTA notes}

{pstd}
The {cmd:.irf} file is a dataset in OpenDTA's own layout. Bootstrap standard errors
({cmd:bs}), {cmd:irf} after {cmd:svar} and {cmd:vec}, and {cmd:irf ctable}, {cmd:cgraph},
{cmd:ograph} and {cmd:add} are not available yet.
{p_end}

{smcl}
{title:Title}

{p 4 4 2}
{bf:stcrreg} {hline 2} Competing-risks regression (Fine and Gray)

{title:Syntax}

{p 8 16 2}{cmd:stcrreg} [{it:indepvars}] [if] [in] {cmd:,} {cmd:compete(}{it:varname} {cmd:==} {it:numlist}{cmd:)} [{cmd:noshr} {cmd:noshow} {cmd:noheader} {cmd:nolog} {cmd:level(}#{cmd:)} {cmd:vce(}{it:vcetype}{cmd:)}]{p_end}

{title:Description}

{pstd}
{cmd:stcrreg} fits the Fine-Gray regression of the subdistribution hazard on {cmd:stset} data, where {cmd:compete()} names the events that compete with the failure event. Coefficients are reported as subhazard ratios (SHR) unless {cmd:noshr} is given. The variance is the robust sandwich estimate, clustered on the subject id by default; {cmd:vce(cluster} {it:varname}{cmd:)} and {cmd:vce(robust)} are accepted.

{pstd}
{cmd:tvc()}, {cmd:texp()}, {cmd:offset()}, {cmd:constraints()}, {cmd:collinear} and {cmd:noadjust} are not available yet and stop with error 198. After {cmd:stcrreg}, {cmd:predict} gives {opt xb} and {opt shr}, and {helpb stcurve} plots the cumulative incidence function.

{title:See also}

{p 4 4 2}{help stcox}, {help stcurve}, {help stset}

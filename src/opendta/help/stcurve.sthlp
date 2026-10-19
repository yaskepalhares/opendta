{smcl}
{title:Title}

{p 4 4 2}
{bf:stcurve} {hline 2} Plot the survivor, hazard, cumulative hazard or cumulative incidence function after stcox, streg or stcrreg

{title:Syntax}

{p 8 16 2}{cmd:stcurve} {cmd:,} {c -(}{cmd:survival}|{cmd:hazard}|{cmd:cumhaz}|{cmd:cif}{c )-} [{cmd:at(}{it:varname=#} ...{cmd:)} {cmd:at1(}...{cmd:)} ... {cmd:at10(}...{cmd:)} {cmd:range(}# #{cmd:)} {cmd:alpha1} {cmd:noshow} {it:graph_options}]{p_end}

{title:Description}

{pstd}
{cmd:stcurve} draws the requested function at the covariate values given in {cmd:at()}; covariates not given are held at their means in the estimation sample. {cmd:cif} is available only after {cmd:stcrreg}, {cmd:survival} and {cmd:hazard} only after {cmd:streg} (and {cmd:survival} and {cmd:cumhaz} after {cmd:stcox}). After {cmd:streg} with {cmd:frailty()} or {cmd:shared()} the curves are conditional on a frailty of one ({cmd:alpha1}); the unconditional curves are not available yet.

{pstd}
{cmd:outfile()}, {cmd:width()}, {cmd:kernel()}, {cmd:noboundary}, {cmd:fixedonly}, {cmd:unconditional}, {cmd:marginal} and {cmd:addplot()} are not available yet and stop with error 198.

{title:See also}

{p 4 4 2}{help stcox}, {help streg}, {help stcrreg}, {help sts}

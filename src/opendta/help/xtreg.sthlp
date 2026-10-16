{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:xtreg} {hline 2}}Panel data: declaration and linear models{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:xtset} {it:panelvar} [{it:timevar}] [{cmd:,} {opt delta(#)}]{p_end}
{p 8 16 2}{cmd:tsset} {it:timevar} [{cmd:,} {opt delta(#)}]{p_end}
{p 8 16 2}{cmd:xtreg} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {cmd:re}|{cmd:fe}|{cmd:be}|{cmd:mle}|{cmd:pa} {it:options}]{p_end}

{title:Description}

{pstd}
{cmd:xtset} declares the panel and time variables used by the {cmd:xt} commands and by the time-series operators {cmd:L.}, {cmd:F.}, {cmd:D.} and {cmd:S.}; {cmd:tsset} declares a single time series. Without arguments, both show the current settings.
{p_end}

{pstd}
{cmd:xtreg} fits random-effects (GLS, the default), fixed-effects (within), between-effects (regression on the panel means; {cmd:wls} weights each panel by its size), maximum-likelihood random-effects and population-averaged ({cmd:pa}, GEE with an exchangeable correlation; see {help xtgee}) models. After {cmd:fe}, {cmd:re} and {cmd:be}, {cmd:predict} offers {cmd:xb}, {cmd:stdp}, {cmd:ue}, {cmd:u}, {cmd:e} and {cmd:xbu}. See also {help xtgee}, {help xtsum}, {help xtdescribe}, {help xtlogit} and {help xtpoisson}.
{p_end}

{title:Options}

{phang}
{opt re} {opt fe} {opt be} {opt mle} {opt pa} the estimator.
{p_end}

{phang}
{opt vce(vcetype)} {cmd:oim}, {cmd:robust}, {cmd:cluster} {it:clustvar} or {cmd:bootstrap}.
{p_end}

{phang}
{opt l:evel(#)} confidence level; default is {cmd:level(95)}.
{p_end}

{title:Examples}

{phang2}{cmd:. xtset id year}{p_end}
{phang2}{cmd:. xtreg ln_wage age tenure, fe}{p_end}
{phang2}{cmd:. xtreg ln_wage age tenure, re vce(robust)}{p_end}

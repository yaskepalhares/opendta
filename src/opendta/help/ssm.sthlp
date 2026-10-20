{smcl}
{title:Title}

{p 4 4 2}
{bf:ucm, sspace, dfactor} {hline 2} State-space models

{title:Syntax}

{p 8 16 2}{cmd:ucm} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {opt model(}{it:none|ntrend|dconstant|llevel|rwalk|dtrend|lldtrend|rwdrift|lltrend|strend|rtrend}{cmd:)} {opt seasonal(#)} {opt cycle(#)} {opt level(#)} {opt nolog}]{p_end}
{p 8 16 2}{cmd:sspace} {cmd:(}{it:state eq}{cmd:)} ... {cmd:(}{it:observation eq}{cmd:)} [{it:if}] [{it:in}] [{cmd:,} {opt diffuse} {opt level(#)} {opt nolog}]{p_end}
{p 8 16 2}{cmd:dfactor} {cmd:(}{it:depvars} {cmd:=} [{it:exog}] [{cmd:,} {opt noconstant}]{cmd:)} {cmd:(}{it:factor} {cmd:=} {cmd:,} {opt ar(numlist)}{cmd:)} [{it:if}] [{it:in}] [{cmd:,} {opt level(#)} {opt nolog}]{p_end}

{title:Description}

{pstd}
All three fit linear Gaussian state-space models by maximum likelihood with the Kalman filter; nonstationary states are diffuse (the initial variance is 10^7 and the observations that depend on diffuse states are left out of the likelihood, VERIFY). The variance parameters are estimated on the log scale and their standard errors come from the delta method. {cmd:ucm} fits a trend (default {cmd:rwalk}), an optional seasonal component and up to three stochastic cycles of order 1 to 3. {cmd:sspace} takes state equations ({it:z L.z e.z, state}) and observation equations ({it:y z e.y}); coefficients can be fixed ({it:z@1}) or shared ({it:L.z@a}). {cmd:dfactor} fits one dynamic factor with unit innovation variance, loadings in each observed equation and diagonal observation errors.
{p_end}

{pstd}
Not available: {cmd:constraints()}, {cmd:from()} and {cmd:vce()}; in {cmd:ucm}, the {cmd:frequency()} suboption of {cmd:cycle()}; in {cmd:sspace}, {cmd:covstate()} and {cmd:covobserved()} other than {cmd:diagonal} and states that depend on lags higher than 1; in {cmd:dfactor}, more than one factor, {cmd:ma()}, error autoregressions and {cmd:covariance()}; {cmd:predict} and {cmd:estat} after them.
{p_end}

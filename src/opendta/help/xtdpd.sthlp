{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:xtdpd} {hline 2}}Panel data: dynamic panel-data estimation (Arellano-Bond and system GMM){p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:xtabond} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {opt lags(#)} {opt maxldep(#)} {opt maxlags(#)} {opt pre(varlist)} {opt endogenous(varlist)} {opt inst(varlist)} {opt twostep} {opt vce(robust)} {opt noconstant}]{p_end}
{p 8 16 2}{cmd:xtdpdsys} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:xtdpd} {it:depvar} {it:indepvars} [{it:if}] [{it:in}] {cmd:,} {opt dgmmiv(varlist)} [{opt lgmmiv(varlist)} {opt div(varlist)} {opt liv(varlist)} {opt twostep} {opt vce(robust)} {opt noconstant}]{p_end}
{p 8 16 2}{cmd:estat abond} [{cmd:,} {opt artests(#)}]{p_end}
{p 8 16 2}{cmd:estat sargan}{p_end}

{title:Description}

{pstd}
{cmd:xtabond} fits the Arellano-Bond difference-GMM estimator for dynamic panel-data models; {cmd:xtdpdsys} adds the level equation (Arellano-Bover/Blundell-Bond system GMM); {cmd:xtdpd} is the general command in which the model and the instruments are specified directly. With {cmd:twostep} and {cmd:vce(robust)} the standard errors have the Windmeijer correction.
{p_end}

{pstd}
{cmd:estat abond} tests for serial correlation in the first-differenced errors and {cmd:estat sargan} for the overidentifying restrictions (not valid after {cmd:vce(robust)}).
{p_end}

{title:Options}

{phang}
{opt lags(#)} lags of the dependent variable; default is 1.
{p_end}

{phang}
{opt maxldep(#)} maximum lag of the dependent variable used as instrument.
{p_end}

{phang}
{opt pre(varlist)} {opt endogenous(varlist)} predetermined and endogenous regressors.
{p_end}

{phang}
{opt twostep} two-step estimator; the default is one-step.
{p_end}

{title:Examples}

{phang2}{cmd:. xtset id year}{p_end}
{phang2}{cmd:. xtabond n w k, lags(2) vce(robust)}{p_end}
{phang2}{cmd:. xtdpdsys n w k, twostep}{p_end}
{phang2}{cmd:. estat abond}{p_end}

{title:Also see}

{pstd}
See also {help xtreg} and {help xtivreg}.
{p_end}

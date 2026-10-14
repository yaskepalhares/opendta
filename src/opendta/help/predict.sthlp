{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:predict} {hline 2}}Predictions after estimation{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:predict} [{it:type}] {it:newvar} [{it:if}] [{it:in}] [{cmd:,} {it:statistic}]{p_end}

{title:Description}

{pstd}
Creates a variable with predictions, residuals or related statistics from the last estimates. Without {cmd:if} or {cmd:in}, predictions are made for every observation, including those outside the estimation sample (use {cmd:if e(sample)} to restrict).
{p_end}

{title:Options}

{phang}
{opt xb} linear prediction (default after regress).
{p_end}

{phang}
{opt pr} probability of a positive outcome (default after logit, probit; with {opt outcome(#)} after ologit, oprobit, mlogit).
{p_end}

{phang}
{opt n} number of events (default after poisson, nbreg).
{p_end}

{phang}
{opt mu} mean of the response (glm).
{p_end}

{phang}
{opt r:esiduals} residuals.
{p_end}

{phang}
{opt stdp} {opt stdf} {opt stdr} standard errors of the prediction, of the forecast and of the residual.
{p_end}

{phang}
{opt hat} {opt rsta:ndard} {opt rstu:dent} {opt cooksd} leverage and influence statistics after regress.
{p_end}

{phang}
{opt sc:ore} equation-level scores.
{p_end}

{title:Examples}

{phang2}{cmd:. regress price mpg weight}{p_end}
{phang2}{cmd:. predict yhat}{p_end}
{phang2}{cmd:. predict r, residuals}{p_end}
{phang2}{cmd:. logit low age smoke}{p_end}
{phang2}{cmd:. predict p if e(sample)}{p_end}

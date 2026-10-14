{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:margins} {hline 2}}Marginal means, predictive margins and marginal effects{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:margins} [{it:marginlist}] [{it:if}] [{it:in}] [{cmd:,} {it:options}]{p_end}

{title:Description}

{pstd}
Computes predictive margins (averages of predictions with the factor variables in {it:marginlist} set to each level) and average marginal effects, with delta-method standard errors. {help marginsplot} draws the last results.
{p_end}

{title:Options}

{phang}
{opt dydx(varlist)} marginal effects (derivatives) of the variables.
{p_end}

{phang}
{opt eyex()} {opt dyex()} {opt eydx()} elasticities and semi-elasticities.
{p_end}

{phang}
{opt at(...)} fix covariates at values, e.g. {cmd:at(age=(20 40 60))}.
{p_end}

{phang}
{opt atmeans} set the other covariates at their means.
{p_end}

{phang}
{opt predict(...)} the prediction to average, e.g. {cmd:predict(xb)}.
{p_end}

{phang}
{opt over(varlist)} margins computed within groups.
{p_end}

{phang}
{opt l:evel(#)} confidence level; default is {cmd:level(95)}.
{p_end}

{phang}
{opt post} post the margins as estimation results.
{p_end}

{title:Examples}

{phang2}{cmd:. logit low age i.smoke}{p_end}
{phang2}{cmd:. margins smoke}{p_end}
{phang2}{cmd:. margins, dydx(age)}{p_end}
{phang2}{cmd:. margins smoke, at(age=(20 30))}{p_end}
{phang2}{cmd:. marginsplot}{p_end}

{smcl}
{title:Title}

{p2colset 5 16 18 2}{...}
{p2col :{cmd:sem} {hline 2}}Structural equation models{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:sem} {it:paths} [{it:if}] [{it:in}] [{cmd:,} {opt stand:ardized} {opt cov:ariance(e.y1*e.y2 ...)}
{opt var:iance(name@#)} {opt lat:ent(names)} {opt nocapslatent} {opt meth:od(ml|mlmv|adf)}
{opt gr:oup(varname)} {opt ginv:ariant(classes)} {opt l:evel(#)} {cmd:nolog}]{p_end}

{p 8 16 2}{it:paths}: {cmd:(}{it:y} {cmd:<-} {it:x1 x2}{cmd:)}, {cmd:(}{it:x1 x2} {cmd:->} {it:y}{cmd:)},
{cmd:(}{it:F} {cmd:->} {it:a b c}{cmd:)}, {cmd:(}{it:y} {cmd:<-} {it:x}{cmd:@1)} ...{p_end}

{p 8 16 2}{cmd:estat gof} [{cmd:,} {opt stats(all|chi2|rmsea|ic|indices|residuals)}]{p_end}
{p 8 16 2}{cmd:estat} {cmd:eqgof}|{cmd:eqtest}|{cmd:framework}|{cmd:ggof}|{cmd:ginvariant}|{cmd:mindices}|{cmd:residuals}|{cmd:scoretests}|{cmd:stable}|{cmd:stdize:}|{cmd:summarize}|{cmd:teffects}{p_end}
{p 8 16 2}{cmd:predict} {it:newvar}{cmd:,} {opt xb(depvar)} | {opt latent(name)}{p_end}

{title:Description}

{pstd}
{cmd:sem} fits linear structural equation models by maximum likelihood: path models with
observed variables, confirmatory factor analysis and models that combine both. Names that
start with a capital letter are latent variables ({cmd:nocapslatent} and {cmd:latent()}
change that). {cmd:method(mlmv)} uses every observation under missing at random;
{cmd:method(adf)} is asymptotic distribution free; {cmd:group()} fits the model in each group with
the parameter classes of {cmd:ginvariant()} constrained equal across groups; {cmd:name@label}
constrains paths with the same label to be equal. Each latent variable is anchored by its first indicator (loading fixed at 1)
and has mean 0. Endogenous variables get an intercept and an error variance; exogenous
variables, observed or latent, have free variances and covariances. {cmd:@}{it:#} fixes a
path at a value, and {cmd:covariance()} adds error covariances.
{p_end}

{pstd}
The output lists the endogenous and exogenous variables, the coefficients by equation
(measurement and structural parts), means, variances and covariances, and the
likelihood-ratio test of the model against the saturated model. {cmd:standardized} reports
the standardized solution. {cmd:estat gof} shows the fit statistics (with {cmd:stats(all)}:
RMSEA, AIC, BIC, CFI, TLI, SRMR and the coefficient of determination). {cmd:predict} gives
the linear prediction of an observed endogenous variable or factor scores (regression
method).
{p_end}

{title:Examples}

{phang2}{cmd:. sem (y <- x1 x2)}{p_end}
{phang2}{cmd:. sem (Anxiety -> a1 a2 a3 a4)}{p_end}
{phang2}{cmd:. sem (F -> x1 x2 x3) (G -> y1 y2 y3) (G <- F), standardized}{p_end}
{phang2}{cmd:. estat gof, stats(all)}{p_end}

{title:OpenDTA notes}

{pstd}
Layouts of the {cmd:estat} subcommands that the manual does not show follow its examples and
may differ in detail from Stata. Nonlinear models are fitted by {help gsem}.
{p_end}

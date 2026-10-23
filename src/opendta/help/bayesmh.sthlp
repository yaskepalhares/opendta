{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:bayesmh} {hline 2}}Bayesian models by Metropolis-Hastings sampling{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:bayesmh} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}]{cmd:,}
{opt likel:ihood(model)} {opt prior(params, dist)} [{opt block(params)} {opt init:ial()}
{opt burn:in(#)} {opt mcmcs:ize(#)} {opt thinn:ing(#)} {opt rseed(#)} {opt clev:el(#)} {opt hpd}
{opt nomodelsummary} {opt dots}]{p_end}

{p 8 16 2}{cmd:bayes} [{cmd:,} {it:options}]{cmd::} {cmd:regress}|{cmd:logit}|{cmd:probit}|{cmd:poisson} {it:depvar} {it:indepvars} ...{p_end}

{p 8 16 2}{cmd:bayesstats summary}|{cmd:ess}|{cmd:ic} [{it:params}] [{cmd:,} {opt hpd} {opt clev:el(#)}]{p_end}
{p 8 16 2}{cmd:bayestest interval} {it:param}{cmd:,} {opt lower(#)} {opt upper(#)}{p_end}
{p 8 16 2}{cmd:bayestest model} {it:name1} {it:name2} ...{p_end}
{p 8 16 2}{cmd:bayesgraph diagnostics}|{cmd:trace}|{cmd:histogram}|{cmd:kdensity}|{cmd:ac} {it:params}{p_end}

{title:Description}

{pstd}
{cmd:bayesmh} samples the posterior distribution of a model with an adaptive random-walk
Metropolis-Hastings sampler. Parameters are written in braces: {cmd:{c -(}}{it:depvar}{cmd::}{it:varlist}{cmd:{c )-}}
is the linear form of the equation, {cmd:{c -(}}{it:depvar}{cmd::}{cmd:{c )-}} all of its coefficients and
{cmd:{c -(}var{c )-}} a scalar parameter. Likelihoods: {cmd:normal()}, {cmd:lognormal()},
{cmd:t()}, {cmd:exponential}, {cmd:poisson}, {cmd:logit}, {cmd:probit} and {cmd:binomial()}.
Priors: {cmd:normal()}, {cmd:lognormal()}, {cmd:uniform()}, {cmd:igamma()}, {cmd:gamma()},
{cmd:beta()}, {cmd:exponential()}, {cmd:chi2()}, {cmd:flat} and {cmd:jeffreys}.
{p_end}

{pstd}
The table reports the posterior mean, standard deviation, Monte Carlo standard error, median
and equal-tailed (or {cmd:hpd}) credible interval, with the acceptance rate, the efficiency
and the log marginal likelihood (Laplace-Metropolis). The {cmd:bayes} prefix fits the listed
estimation commands with default priors. {cmd:bayesstats}, {cmd:bayestest} and {cmd:bayesgraph}
summarize, test and plot the draws.
{p_end}

{title:Examples}

{phang2}{cmd:. bayesmh y x, likelihood(normal({c -(}var{c )-})) prior({c -(}y:{c )-}, normal(0,100)) prior({c -(}var{c )-}, igamma(.01,.01)) rseed(17)}{p_end}
{phang2}{cmd:. bayes, rseed(17): logit d x1 x2}{p_end}
{phang2}{cmd:. bayesstats summary {c -(}y:x{c )-}}{p_end}
{phang2}{cmd:. bayestest interval {c -(}y:x{c )-}, lower(0)}{p_end}

{title:OpenDTA notes}

{pstd}
The draws come from OpenDTA's own random numbers and differ from Stata's; the summaries agree
up to Monte Carlo error. The {cmd:bayes} prefix accepts {cmd:regress}, {cmd:logit},
{cmd:probit} and {cmd:poisson} only.
{p_end}

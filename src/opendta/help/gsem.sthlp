{smcl}
{title:Title}

{p2colset 5 17 19 2}{...}
{p2col :{cmd:gsem} {hline 2}}Generalized structural equation models{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:gsem} {it:paths} [{it:if}] [{it:in}] [{cmd:,} {opt cov:ariance(L1*L2 ...)}
{opt var:iance(name@#)} {opt lat:ent(names)} {opt nocapslatent} {opt intm:ethod(mvaghermite|mcaghermite|ghermite)}
{opt intp:oints(#)} {opt l:evel(#)} {opt iter:ate(#)} {cmd:nolog}]{p_end}

{p 8 16 2}{it:paths}: {cmd:(}{it:y} {cmd:<-} {it:x1 x2 L}{cmd:,} {opt fam:ily(name)} {opt li:nk(name)}{cmd:)},
{cmd:(}{it:L} {cmd:->} {it:y1 y2 y3}{cmd:, logit)}, {cmd:(}{it:y} {cmd:<-} {it:x} {cmd:M1[}{it:id}{cmd:])} ...{p_end}

{p 8 16 2}{cmd:predict} {it:newvar}{cmd:,} {opt mu(depvar)} | {opt xb(depvar)} | {opt lat:ent(name)}{p_end}

{title:Description}

{pstd}
{cmd:gsem} fits structural equation models whose responses need not be normal. Each response
has a family ({cmd:gaussian}, {cmd:bernoulli}, {cmd:binomial}, {cmd:poisson}, {cmd:nbinomial},
{cmd:gamma}, {cmd:ordinal}) and a link; the shorthands {cmd:regress}, {cmd:logit},
{cmd:probit}, {cmd:cloglog}, {cmd:poisson}, {cmd:nbreg}, {cmd:ologit}, {cmd:oprobit} and
{cmd:gamma} set both. Names that start with a capital letter are latent variables; {cmd:M[}{it:id}{cmd:]}
is a latent variable at the level of the clusters of {it:id} (a random effect). The first path
of each latent variable is fixed at 1; latent variables are uncorrelated unless
{cmd:covariance()} says otherwise.
{p_end}

{pstd}
The latent variables are integrated out by Gauss-Hermite quadrature, adaptive by default, and
the likelihood is maximized by Newton steps. The table lists the coefficients by response, the
ancillary parameters (cutpoints, variances, dispersion) and the variances and covariances of the
latent variables. {cmd:predict} gives the mean or the linear prediction of a response, or the
empirical Bayes means of a latent variable.
{p_end}

{title:Examples}

{phang2}{cmd:. gsem (y <- x1 x2, logit)}{p_end}
{phang2}{cmd:. gsem (L -> y1 y2 y3 y4, logit)}{p_end}
{phang2}{cmd:. gsem (y <- x M1[school], poisson)}{p_end}
{phang2}{cmd:. predict mu1, mu(y1)}{p_end}

{title:OpenDTA notes}

{pstd}
A model cannot mix observation-level latent variables with cluster-level ones. The iteration
log and the layout of the ancillary parameters follow the manual's examples and may differ in
detail from Stata.
{p_end}

{smcl}
{title:Title}

{p 4 4 2}
{bf:hotelling}, {bf:mvreg}, {bf:manova}, {bf:canon}, {bf:mvtest} {hline 2} Multivariate tests and models

{title:Syntax}

{p 8 16 2}
{cmd:hotelling} varlist [if] [in] [weight] [{cmd:,} {cmd:by(}varname{cmd:)} {cmd:notable}]

{p 8 16 2}
{cmd:mvreg} depvars {cmd:=} indepvars [if] [in] [weight] [{cmd:,} {cmd:noconstant} {cmd:corr} {cmd:level(}#{cmd:)}]

{p 8 16 2}
{cmd:manova} depvars {cmd:=} termlist [if] [in] [weight] [{cmd:,} {cmd:sequential} {cmd:continuous(}varlist{cmd:)} {cmd:noconstant} {cmd:partial} {cmd:dropemptycells}]

{p 8 16 2}
{cmd:canon} ({it:varlist1}) ({it:varlist2}) [if] [in] [weight] [{cmd:,} {cmd:lc(}#{cmd:)} {cmd:noconstant} {cmd:stdcoef} {cmd:test(}numlist{cmd:)} {cmd:level(}#{cmd:)}]

{p 8 16 2}
{cmd:mvtest} {cmd:means} varlist [{cmd:,} {cmd:by(}varname{cmd:)}]   {cmd:mvtest covariances} | {cmd:correlations} varlist [{cmd:,} {cmd:by(}varname{cmd:)}]   {cmd:mvtest normality} varlist [{cmd:,} {cmd:stats(}mardia hz{cmd:)}]

{title:Description}

{pstd}
{cmd:hotelling} tests that the means of varlist are zero (one sample) or equal across the two groups of {cmd:by()}, with Hotelling's T-squared and its F transformation. {cmd:mvreg} fits a multivariate regression (seemingly unrelated equations with equal regressors) and posts {cmd:e(b)}, {cmd:e(V)} and {cmd:e(Sigma)}. {cmd:manova} reports Wilks' lambda, Pillai's trace, Lawley-Hotelling's trace and Roy's largest root with their F approximations (flags {cmd:e} exact, {cmd:a} approximate, {cmd:u} upper bound); after it, {cmd:manovatest} term [term ...] [{cmd:/} term ...] [{cmd:,} {cmd:ytransform(}matname{cmd:)}] tests the listed terms jointly, against the residual or the terms after the slash ({cmd:test()} and {cmd:showorder} are not implemented yet). {cmd:canon} reports canonical correlations and the coefficients of the canonical variates. {cmd:mvtest} tests means, covariance matrices (Box's M and likelihood-ratio tests of diagonal or spherical structure), correlation matrices and multivariate normality (Mardia's skewness and kurtosis, Henze-Zirkler).

{pstd}
Not implemented yet: {cmd:manova} with {cmd:/} error terms, {cmd:mvtest means} with {cmd:pairwise} and {cmd:mvtest correlations, pairwise}, {cmd:mvtest covariances} with {cmd:compound} and {cmd:block}, {cmd:mvtest normality} with {cmd:univariate}, {cmd:bivariate}, {cmd:by()} and Doornik-Hansen.

{title:Postestimation}

{phang2}{cmd:predict} newvar [{cmd:,} {cmd:xb} {cmd:residuals} {cmd:equation(}eqno{cmd:)}] (after {cmd:mvreg}){p_end}

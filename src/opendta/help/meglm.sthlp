{smcl}
{title:Title}

{p 4 4 2}
{bf:meglm}, {bf:mecloglog}, {bf:menbreg}, {bf:meologit}, {bf:meoprobit}, {bf:meqrlogit}, {bf:meqrpoisson} {hline 2} Multilevel mixed-effects generalized linear models

{title:Syntax}

{p 8 16 2}
{cmd:meglm} depvar [fe_varlist] [if] [in] {cmd:||} {it:levelvar}{cmd::} [re_varlist] [{cmd:,} {cmd:family(}{it:family}{cmd:)} {cmd:link(}{it:link}{cmd:)} {cmd:covariance(}...{cmd:)} {cmd:noconstant} {it:options}]

{p 8 16 2}
{cmd:mecloglog}, {cmd:menbreg}, {cmd:meologit}, {cmd:meoprobit} depvar [fe_varlist] ... {cmd:||} {it:levelvar}{cmd::} ... [{cmd:,} {it:options}]

{p 8 16 2}
{cmd:meqrlogit}, {cmd:meqrpoisson} depvar [fe_varlist] ... {cmd:||} {it:levelvar}{cmd::} ... [{cmd:,} {it:options}]

{title:Description}

{pstd}
Fit generalized linear mixed models with one level of random effects (intercepts and slopes). {cmd:meglm} accepts {cmd:family()} gaussian, bernoulli, binomial, gamma, nbinomial (mean or constant dispersion), ordinal and poisson, with the usual links. {cmd:mecloglog} is the binary complementary log-log model, {cmd:menbreg} the negative binomial (option {cmd:dispersion(mean|constant)}), {cmd:meologit} and {cmd:meoprobit} the ordered models. {cmd:meqrlogit} and {cmd:meqrpoisson} are the QR-decomposition commands; they fit the same models as {cmd:melogit} and {cmd:mepoisson}.

{pstd}
The likelihood of each group is integrated by mean-variance adaptive Gauss-Hermite quadrature ({cmd:intmethod(mvaghermite)}, 7 points by default), non-adaptive Gauss-Hermite or the Laplace approximation ({cmd:intmethod(laplace)}). The fixed-effects model is fitted first and an LR test against it is reported.

{title:Options}

{phang}{cmd:intpoints(}#{cmd:)}, {cmd:intmethod(}...{cmd:)}, {cmd:laplace}: integration. {cmd:variance}, {cmd:stddeviations}: parameterization of the random-effects table. {cmd:level(}#{cmd:)}, {cmd:nolog}, {cmd:or}, {cmd:irr}, {cmd:offset(}varname{cmd:)}, {cmd:exposure(}varname{cmd:)}: as in {help melogit}.{p_end}

{title:Postestimation}

{phang2}{cmd:predict} newvar [, {cmd:mu} {cmd:xb} {cmd:eta} {cmd:reffects} {cmd:fixedonly}]{p_end}
{phang2}{cmd:estat icc}, {cmd:estat recovariance}, {cmd:estat group}{p_end}

{title:Limitations}

{pstd}
Only one level of random effects; weights and {cmd:vce(robust)} are not available.

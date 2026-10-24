{smcl}
{title:Title}

{p 4 4 2}
{bf:intreg}, {bf:truncreg}, {bf:heckman}, {bf:heckprobit}, {bf:heckoprobit}, {bf:biprobit} {hline 2} Interval, truncated, selection and bivariate probit models

{title:Syntax}

{p 8 16 2}
{cmd:intreg} depvar1 depvar2 [indepvars] [if] [in] [weight] [{cmd:,} {cmd:het(}varlist{cmd:)} ...]

{p 8 16 2}
{cmd:truncreg} depvar [indepvars] [if] [in] [weight] [{cmd:,} {cmd:ll(}#{cmd:)} {cmd:ul(}#{cmd:)} ...]

{p 8 16 2}
{cmd:heckman} depvar [indepvars] [if] [in] [weight] {cmd:,} {cmd:select(}[depvar_s {cmd:=}] varlist_s{cmd:)} [{cmd:twostep} {cmd:mills(}newvar{cmd:)}|{cmd:mle}] [{cmd:offset(}varname{cmd:)}] ...

{p 8 16 2}
{cmd:heckprobit} depvar [indepvars] [if] [in] {cmd:,} {cmd:select(}...{cmd:)} [...]

{p 8 16 2}
{cmd:heckoprobit} depvar [indepvars] [if] [in] {cmd:,} {cmd:select(}...{cmd:)} [...]

{p 8 16 2}
{cmd:biprobit} depvar1 depvar2 [indepvars] [if] [in] [weight] [{cmd:,} {cmd:offset1(}varname{cmd:)} {cmd:offset2(}varname{cmd:)} ...] ({cmd:partial} is not available yet){p_end}
{p 8 16 2}
{cmd:biprobit} {cmd:(}depvar1 {cmd:=} indepvars1{cmd:)} {cmd:(}depvar2 {cmd:=} indepvars2{cmd:)} [if] [in] [weight] [{cmd:,} ...]

{title:Description}

{pstd}
Models with censored, truncated or selected outcomes fitted by maximum likelihood. The bivariate normal distribution function used by {cmd:biprobit}, {cmd:heckprobit} and {cmd:heckoprobit} is computed with Plackett's identity and Gauss-Legendre quadrature, with exact derivatives. {cmd:heckman} offers the two-step estimator with Greene's corrected variance. The correlation is estimated as {cmd:/athrho} (atanh rho), {cmd:/lnsigma} and {cmd:/lambda} are reported through the delta method.

{title:Postestimation}

{phang2}{cmd:predict} newvar [{cmd:,} {cmd:xb} {cmd:stdp} {cmd:pmargin} {cmd:p11} {cmd:ycond} {cmd:yexpected} {cmd:ystar(}#,#{cmd:)} ...]{p_end}
{phang2}{cmd:test}, {cmd:lincom}, {cmd:margins}, {cmd:estat ic}{p_end}

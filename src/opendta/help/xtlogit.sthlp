{smcl}
{title:Title}

{p 4 4 2}
{bf:xtlogit} {hline 2} Fixed-effects and random-effects logit models

{title:Syntax}

{p 8 16 2}
{cmd:xtlogit} depvar indepvars [if] [in] [, fe re or]

{title:Description}

{pstd}
{cmd:fe} fits the conditional (fixed-effects) logit model; groups with all-positive or all-negative outcomes are dropped. {cmd:re} (default) integrates a normal random effect by adaptive Gauss-Hermite quadrature with 12 points and reports sigma_u, rho and the LR test of rho=0.

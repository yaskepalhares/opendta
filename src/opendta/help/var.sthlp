{smcl}
{title:Title}

{p 4 4 2}
{bf:var} {hline 2} Vector autoregressive models

{title:Syntax}

{p 8 16 2}
{cmd:var} depvarlist [if] [in] [{cmd:,} {cmd:lags(}numlist{cmd:)} {cmd:exog(}varlist{cmd:)} {cmd:noconstant} {cmd:dfk} {cmd:small} {cmd:level(}#{cmd:)}]

{p 8 16 2}
{cmd:varsoc} depvarlist [if] [in] [{cmd:,} {cmd:maxlag(}#{cmd:)} {cmd:exog(}varlist{cmd:)} {cmd:noconstant}]

{p 8 16 2}
{cmd:vargranger}

{title:Description}

{pstd}
{cmd:var} fits a VAR by OLS equation by equation (lags 1 and 2 by default) and reports the log likelihood, FPE, AIC, HQIC, SBIC, the determinant of the ML covariance of the residuals and per-equation statistics. {cmd:varsoc} reports the lag-order selection statistics over a common sample, marking with * the order chosen by each criterion. {cmd:vargranger} reports Granger causality Wald tests after {cmd:var}.

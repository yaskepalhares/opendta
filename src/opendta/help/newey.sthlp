{smcl}
{title:Title}

{p 4 4 2}
{bf:newey} {hline 2} Regression with Newey-West standard errors

{title:Syntax}

{p 8 16 2}
{cmd:newey} depvar [indepvars] [if] [in] [weight]{cmd:,} {cmd:lag(}#{cmd:)} [{cmd:force}]

{title:Description}

{pstd}
OLS coefficients with standard errors robust to heteroskedasticity and autocorrelation up to lag #, with Bartlett weights. The data must be {cmd:tsset}.

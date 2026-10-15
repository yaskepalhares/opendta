{smcl}
{title:Title}

{p 4 4 2}
{bf:dfuller} {hline 2} Augmented Dickey-Fuller unit-root test

{title:Syntax}

{p 8 16 2}
{cmd:dfuller} varname [if] [in] [{cmd:,} {cmd:noconstant} {cmd:trend} {cmd:drift} {cmd:regress} {cmd:lags(}#{cmd:)}]

{title:Description}

{pstd}
Regresses D.y on L.y (and lags of D.y, a constant, a trend) and tests whether the coefficient on L.y is zero. Critical values are interpolated in the number of observations from Fuller (1976); the approximate p-value is MacKinnon's (1994). With {cmd:drift}, the statistic has a t distribution.

{smcl}
{title:Title}

{p 4 4 2}
{bf:corrgram} {hline 2} Tabulate the correlogram

{title:Syntax}

{p 8 16 2}
{cmd:corrgram} varname [if] [in] [{cmd:,} {cmd:lags(}#{cmd:)} {cmd:noplot} {cmd:yw}]

{title:Description}

{pstd}
Lists the autocorrelations, the partial autocorrelations (from regressions on the lags, or Yule-Walker with {cmd:yw}), the Ljung-Box Q statistics and their p-values, with character plots of the two correlations.

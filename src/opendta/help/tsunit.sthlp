{smcl}
{title:Title}

{p 4 4 2}
{bf:dfgls, pperron} {hline 2} Unit-root tests

{title:Syntax}

{p 8 16 2}{cmd:dfgls} {it:varname} [{it:if}] [{it:in}] [{cmd:,} {opt maxlag(#)} {opt notrend} {opt ers}]{p_end}
{p 8 16 2}{cmd:pperron} {it:varname} [{it:if}] [{it:in}] [{cmd:,} {opt noconstant} {opt trend} {opt regress} {opt lags(#)}]{p_end}

{title:Description}

{pstd}
{cmd:dfgls} performs the Elliott-Rothenberg-Stock modified Dickey-Fuller test on the GLS-detrended series for each lag up to {cmd:maxlag()} and reports the lag chosen by the Ng-Perron sequential rule, the Schwarz criterion and the modified Akaike criterion. {cmd:pperron} performs the Phillips-Perron test with Newey-West variance. The 5% and 10% critical values of {cmd:dfgls} are interpolated from the tables of Elliott et al. and Cheung and Lai; the table for the mean-only case is the asymptotic one (marked VERIFY in the source).
{p_end}

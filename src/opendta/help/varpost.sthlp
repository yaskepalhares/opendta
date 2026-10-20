{smcl}
{title:Title}

{p 4 4 2}
{bf:varlmar, varnorm, varstable, varwle, varbasic, fcast} {hline 2} Postestimation for VAR, SVAR and VEC

{title:Syntax}

{p 8 16 2}{cmd:varlmar} [{cmd:,} {opt mlag(#)} {opt estimates(name)}]{p_end}
{p 8 16 2}{cmd:varnorm} [{cmd:,} {opt jbera} {opt skewness} {opt kurtosis} {opt cholesky} {opt estimates(name)}]{p_end}
{p 8 16 2}{cmd:varstable} [{cmd:,} {opt graph} {opt amat(name)} {opt estimates(name)} {it:graph options}]{p_end}
{p 8 16 2}{cmd:varwle} [{cmd:,} {opt estimates(name)}]{p_end}
{p 8 16 2}{cmd:varbasic} {it:varlist} [{it:if}] [{it:in}] [{cmd:,} {opt lags(numlist)} {opt irf} {opt fevd} {opt nograph} {opt step(#)}]{p_end}
{p 8 16 2}{cmd:fcast} {c -(}{opt compute} {it:prefix} {cmd:,} {opt step(#)} {opt dynamic(time)} {opt nose}{c )-} | {cmd:fcast graph} {it:varlist} [{cmd:,} {opt observed}]{p_end}

{title:Description}

{pstd}
{cmd:varlmar} is the Lagrange-multiplier test for residual autocorrelation, {cmd:varnorm} the Jarque-Bera, skewness and kurtosis tests on the orthogonalized residuals, {cmd:varstable} the eigenvalue stability condition (with the unit-circle graph) and {cmd:varwle} the Wald tests of exclusion of each lag. {cmd:varbasic} fits a VAR and graphs the orthogonalized impulse responses ({cmd:irf} the simple ones, {cmd:fevd} the variance decomposition). {cmd:fcast compute} creates the forecast variables after {cmd:var}, {cmd:svar} or {cmd:vec} and extends the data when the forecast goes beyond the last period; {cmd:fcast graph} plots them.
{p_end}

{pstd}
Not available: {cmd:fcast} with {cmd:bs}, {cmd:bsp} and {cmd:exog} variables, and {cmd:fcast} after a VAR with non-consecutive lags.
{p_end}

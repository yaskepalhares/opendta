{smcl}
{title:Title}

{p 4 4 2}
{bf:arch, arfima} {hline 2} ARCH-family and fractionally integrated models

{title:Syntax}

{p 8 16 2}{cmd:arch} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {opt arch(numlist)} {opt garch(numlist)} {opt saarch()} {opt tarch()} {opt aarch()} {opt narch()} {opt narchk()} {opt abarch()} {opt atarch()} {opt sdgarch()} {opt earch()} {opt egarch()} {opt parch()} {opt tparch()} {opt aparch()} {opt nparch()} {opt nparchk()} {opt pgarch()} {opt archm} {opt archmlags()} {opt archmexp()} {opt ar()} {opt ma()} {opt arima(#p,#d,#q)} {opt distribution(gaussian|t|ged [#])} {opt het(varlist)} {opt arch0()} {opt condobs()} {opt noconstant} {opt vce(oim|opg|robust)} {opt level(#)} {opt nolog}]{p_end}
{p 8 16 2}{cmd:arfima} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {opt ar(numlist)} {opt ma(numlist)} {opt noconstant} {opt level(#)} {opt nolog}]{p_end}

{title:Description}

{pstd}
{cmd:arch} fits regressions with ARCH, GARCH and the other conditional-variance terms, ARCH-in-mean, ARMA disturbances and Gaussian, Student t or generalized-error distributions, by maximum likelihood with BHHH steps followed by BFGS ({cmd:technique(bhhh 5 bfgs 10)}); the default variance is OPG. {cmd:narchk()} and {cmd:nparchk()} use one shift parameter per lag. {cmd:constraints()} and {cmd:from()} are not available. {cmd:arfima} fits ARFIMA models by exact maximum likelihood with the observed information matrix; {cmd:mpl}, {cmd:smemory} and robust variances are not available.
{p_end}

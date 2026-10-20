{smcl}
{title:Title}

{p 4 4 2}
{bf:prais} {hline 2} Prais-Winsten and Cochrane-Orcutt regression

{title:Syntax}

{p 8 16 2}{cmd:prais} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {opt corc} {opt rho:type(}{it:regress}|{it:freg}|{it:tscorr}|{it:dw}|{it:theil}|{it:nagar}{cmd:)} {opt ssesearch} {opt two:step} {opt nodw} {opt noconstant} {opt level(#)}]{p_end}

{title:Description}

{pstd}
Fits a linear regression with AR(1) disturbances by iterated feasible GLS. The default is Prais-Winsten (keeps the first observation); {cmd:corc} drops it. {cmd:rhotype()} chooses the estimator of the autocorrelation, {cmd:ssesearch} searches for rho by minimizing the sum of squared residuals and {cmd:twostep} stops after the first transformation. The output ends with the original and the transformed Durbin-Watson statistics. {cmd:vce(robust)} and {cmd:vce(cluster)} are not available.
{p_end}

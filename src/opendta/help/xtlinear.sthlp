{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:xtlinear} {hline 2}}Panel data: FGLS, PCSE, AR(1), random-coefficients, Hausman-Taylor and IV models{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:xtgls} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {opt panels(iid|hetero|correlated)} {opt corr(independent|ar1|psar1)} {opt igls} {it:options}]{p_end}
{p 8 16 2}{cmd:xtpcse} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {opt corr(independent|ar1|psar1)} {opt hetonly}|{opt independent} {opt pairwise} {it:options}]{p_end}
{p 8 16 2}{cmd:xtregar} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}] [{cmd:,} {opt fe}|{opt re} {opt rhotype(type)} {opt twostep} {it:options}]{p_end}
{p 8 16 2}{cmd:xtrc} {it:depvar} [{it:indepvars}] [{it:if}] [{it:in}]{p_end}
{p 8 16 2}{cmd:xthtaylor} {it:depvar} {it:indepvars} [{it:if}] [{it:in}] {cmd:,} {opt endog(varlist)} [{opt constant(varlist)} {opt varying(varlist)} {opt amacurdy}]{p_end}
{p 8 16 2}{cmd:xtivreg} {it:depvar} [{it:varlist1}] {cmd:(}{it:varlist2} {cmd:=} {it:varlist_iv}{cmd:)} [{it:if}] [{it:in}] [{cmd:,} {opt re}|{opt be}|{opt fe}|{opt fd}|{opt ec2sls} {opt nosa} {it:options}]{p_end}

{title:Description}

{pstd}
{cmd:xtgls} fits panel-data linear models by feasible generalized least squares with heteroskedastic or cross-sectionally correlated errors and AR(1) autocorrelation; {cmd:xtpcse} keeps the pooled OLS (or Prais-Winsten) coefficients and computes panel-corrected standard errors. {cmd:xtregar} fits fixed- and random-effects models with AR(1) disturbances, {cmd:xtrc} the Swamy random-coefficients model, {cmd:xthtaylor} the Hausman-Taylor estimator and {cmd:xtivreg} panel instrumental-variables estimators (G2SLS and EC2SLS random effects, between, within and first differences).
{p_end}

{title:Options}

{phang}
{opt i(varname)} panel variable, if different from the one set by {cmd:xtset}.
{p_end}

{phang}
{opt l:evel(#)} confidence level; default is {cmd:level(95)}.
{p_end}

{phang}
{opt rhotype(type)} how the autocorrelation is estimated: {cmd:regress}, {cmd:freg}, {cmd:tscorr}, {cmd:dw}, {cmd:theil} or {cmd:nagar}.
{p_end}

{title:Examples}

{phang2}{cmd:. xtset id year}{p_end}
{phang2}{cmd:. xtgls y x1 x2, panels(hetero) corr(ar1)}{p_end}
{phang2}{cmd:. xtpcse y x1 x2, corr(ar1)}{p_end}
{phang2}{cmd:. xtivreg y x1 (x2 = z), fe}{p_end}
{phang2}{cmd:. xthtaylor y x1 x2 w, endog(x2)}{p_end}

{title:Also see}

{pstd}
See also {help xtreg}, {help xtgee} and {help ivregress}.
{p_end}

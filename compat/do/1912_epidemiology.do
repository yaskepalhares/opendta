* Epidemiology, second part: rocreg with the parametric ROC model fit by maximum likelihood
* (probit ml), with ROC and control covariates, and the cluster() option (cluster-robust standard
* errors after ml, cluster resampling in the bootstrap). Every command whose behavior the manual
* does not settle runs under capture noisily, so that one Stata run shows the whole case.
clear
set seed 1912
set obs 300
gen d = runiform() < .4
gen x = round(runiform(), .01)
gen z = runiform() < .5
gen id = ceil(_n/3)
gen y = round(10 + 2*z + invnormal(runiform())*(1 + .3*d) + d*(1.2 + .8*x), .01)
gen y2 = round(y + invnormal(runiform()), .01)

* probit ml without covariates: the estimates are those of two normal samples
capture noisily rocreg d y, probit ml
capture noisily ereturn list
capture noisily matrix list e(b)
capture noisily matrix list e(V)
* ROC and control covariates
capture noisily rocreg d y2, probit ml roccov(x) ctrlcov(z)
capture noisily ereturn list
capture noisily rocreg d y2, probit ml roccov(x) ctrlcov(z) ctrlmodel(linear)
capture noisily rocreg d y2, probit ml ctrlcov(z) ctrlmodel(strata)
capture noisily rocreg d y y2, probit ml
* cluster-robust standard errors
capture noisily rocreg d y2, probit ml roccov(x) ctrlcov(z) cluster(id)
capture noisily ereturn list
capture noisily rocreg d y, ml
* cluster() with the bootstrap: the draws depend on Stata's resampling order (VERIFY), so the
* layout and e() are what this part checks
capture noisily rocreg d y, cluster(id) breps(50) bseed(5) nodots
capture noisily ereturn list
capture noisily rocreg d y, probit cluster(id) breps(50) bseed(5) nodots
capture noisily rocreg d y, cluster(id) nobstrata breps(50) bseed(5) nodots

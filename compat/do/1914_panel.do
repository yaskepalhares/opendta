* Panel data, second part: xtivreg with first and regress; xtdata with re; xtrc with betas;
* xtunitroot fisher with drift, demean and pperron, and hadri with kernel(); xtabond and xtdpd
* with artests() and hascons. Every command whose behavior the manual does not settle runs under
* capture noisily, so that one Stata run shows the whole case.
clear
set seed 1914
set obs 300
gen id = ceil(_n/10)
gen t = mod(_n-1,10)+1
xtset id t
gen z = round(invnormal(runiform()),.01)
gen x2 = round(invnormal(runiform()),.01)
gen u = round(invnormal(runiform()),.01)
gen x1 = round(z + .5*u + invnormal(runiform()),.01)
gen y = 1 + x1 + x2 + u + id/10

* xtivreg: first and regress
capture noisily xtivreg y x2 (x1 = z), fe first
capture noisily xtivreg y x2 (x1 = z), re first
capture noisily xtivreg y x2 (x1 = z), be first
capture noisily xtivreg y x2 (x1 = z), fd first
capture noisily xtivreg y x2 (x1 = z), fe regress
capture noisily ereturn list
capture noisily xtivreg y x2 (x1 = z), re regress
capture noisily xtreg y x1 x2, re
capture noisily xtivreg y x2 (x1 = z), first regress

* xtrc, betas
preserve
capture noisily keep if id <= 5
capture noisily xtrc y x1 x2, betas
capture noisily ereturn list
capture noisily matrix list e(beta_ps)
restore

* xtdata, re
capture noisily xtreg y x1 x2, re
local r = e(sigma_u)/e(sigma_e)
preserve
capture noisily xtdata y x1 x2, re ratio(`r') clear
capture noisily describe
capture noisily summarize
capture noisily regress y x1 x2 constant, noconstant
restore
capture noisily xtdata y x1 x2, re clear
capture noisily xtdata y x1 x2, fa clear

* xtunitroot
clear
set obs 400
gen id = ceil(_n/20)
gen t = mod(_n-1,20)+1
xtset id t
gen e = round(invnormal(runiform()), .01)
gen y = e
sort id t, stable
by id: replace y = .5*y[_n-1] + e if _n > 1
capture noisily xtunitroot fisher y, dfuller lags(1) drift
capture noisily return list
capture noisily xtunitroot fisher y, dfuller lags(1) demean
capture noisily xtunitroot fisher y, dfuller lags(1) trend demean
capture noisily xtunitroot fisher y, pperron lags(2)
capture noisily return list
capture noisily xtunitroot fisher y, pperron lags(2) trend
capture noisily xtunitroot fisher y, pperron lags(2) demean
capture noisily xtunitroot fisher y, dfuller drift trend
capture noisily xtunitroot hadri y, kernel(bartlett 3)
capture noisily return list
capture noisily xtunitroot hadri y, kernel(parzen 3) robust
capture noisily xtunitroot hadri y, kernel(quadraticspectral 2) trend
capture noisily xtunitroot hadri y, kernel(bartlett)

* xtabond and xtdpd: artests() and hascons
clear
set obs 300
gen id = ceil(_n/10)
gen t = mod(_n-1,10)+1
xtset id t
gen x = round(invnormal(runiform()),.01)
gen y = round(invnormal(runiform()),.01)
sort id t, stable
by id: replace y = .4*y[_n-1] + x + y if _n > 1
gen one = 1
capture noisily xtabond y x, artests(3) vce(robust)
capture noisily ereturn list
capture noisily estat abond
capture noisily xtdpd y L.y x one, dgmmiv(y) lgmmiv(y) div(x) liv(one) hascons
capture noisily xtdpd y L.y x, dgmmiv(y) div(x) artests(3) vce(robust)
capture noisily estat abond

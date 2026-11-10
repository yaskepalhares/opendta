* estat dwatson after regress ([R] regress postestimation time series): layout, r(), a gap in
* the time variable, and the errors without tsset, with panel data and after another model.
clear
set obs 30
gen t = _n
gen x = mod(_n*7, 11)
gen y = 2 + 0.5*x + sin(_n)
quietly regress y x
capture noisily estat dwatson
display _rc
tsset t
estat dwatson
return list
drop if t == 10
quietly regress y x
estat dwatson
quietly regress y x if t > 5
estat dwat
gen id = 1 + (t > 15)
xtset id t
quietly regress y x
capture noisily estat dwatson
display _rc
tsset t
gen d = y > 4
quietly logit d x
capture noisily estat dwatson
display _rc

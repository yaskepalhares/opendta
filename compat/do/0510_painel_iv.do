* Fase 5b: xtreg, areg e ivregress.
clear
estimates clear
ereturn clear
set seed 1010
set obs 60
gen id = ceil(_n/6)
gen t = mod(_n-1, 6) + 1
gen a = id/3
gen x1 = round(runiform()*10, .1) + a
gen x2 = round(runiform()*5, .01)
gen z = x1 + round(runiform()*2, .01)
gen y = 1 + .5*x1 - x2 + 2*a + round(runiform()*2, .01)
xtset id t
xtreg y x1 x2, fe
ereturn list
xtreg y x1 x2, re
xtreg y x1 x2, fe vce(robust)
areg y x1 x2, absorb(id)
ivregress 2sls y x2 (x1 = z)
ivregress 2sls y x2 (x1 = z), small
ereturn list

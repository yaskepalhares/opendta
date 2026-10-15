* Robust and cluster variance with frequency weights: each observation counts w times, so the
* variance is that of the expanded data (OpenDTA's vce.sandwich). VERIFY against Stata 14.
clear
set seed 514
set obs 120
gen x = runiform()
gen g = mod(_n, 10)
gen t = ceil(-ln(runiform())/exp(-1 + x)*10)
gen d = runiform() < .8
gen fw = 1 + mod(_n, 3)
gen k = mod(_n, 3)
gen y = runiform() < invlogit(x)
logit y x [fweight=fw], vce(cluster g) nolog
matrix list e(V), format(%12.0g)
ologit k x [fweight=fw], vce(cluster g) nolog
matrix list e(V), format(%12.0g)
ologit k x [fweight=fw], vce(robust) nolog
matrix list e(V), format(%12.0g)
stset t [fweight=fw], failure(d)
stcox x, vce(cluster g) nolog
matrix list e(V), format(%12.0g)
streg x, distribution(weibull) vce(robust) nolog
matrix list e(V), format(%12.0g)
streg x, distribution(exponential) vce(cluster g) nolog
matrix list e(V), format(%12.0g)

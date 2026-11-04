* Models, part 8 (block 20f): teffects ra, ipw, ipwra and aipw with multivalued treatments. No
* reference until the Stata run. Every command whose behavior the manual does not settle runs
* under capture noisily, so that one Stata run shows the whole case.
clear
set obs 1200
set seed 1925
gen x = round(invnormal(runiform()), .01)
gen z = round(invnormal(runiform()), .01)
gen double u = runiform()
gen double p1 = exp(.3*x + .2*z)
gen double p2 = exp(-.2 + .5*x - .3*z)
gen t = 0
replace t = 1 if u < p1/(1+p1+p2)
replace t = 2 if u >= p1/(1+p1+p2) & u < (p1+p2)/(1+p1+p2)
gen y = round(1 + .5*x + t + .3*t*x + invnormal(runiform()), .01)
capture noisily teffects ra (y x) (t)
capture noisily ereturn list
capture noisily teffects ra (y x) (t), pomeans
capture noisily teffects ipw (y) (t x z)
capture noisily teffects ipwra (y x) (t x z)
capture noisily teffects aipw (y x) (t x z)
capture noisily teffects ipwra (y x) (t x z), atet
capture noisily teffects ipwra (y x) (t x z), atet tlevel(1)
capture noisily teffects ra (y x) (t), control(1)
capture noisily teffects ipw (y) (t x z, probit)

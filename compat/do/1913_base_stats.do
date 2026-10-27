* Base statistics, second part: ml with the lf0, lf1, lf2, gf0 methods, constraints() and svy;
* boxcox with model(rhsonly), model(lambda), model(theta) and notrans(); icc with testvalue() and
* unbalanced one-way data; loneway with mean, median and exact; esize with unequal and welch;
* pwmean with unequal; ameans with weights; lv with generate and tail(); stem with width().
* Every command whose behavior the manual does not settle runs under capture noisily, so that one
* Stata run shows the whole case.
clear
set seed 1913
set obs 200
gen x1 = round(invnormal(runiform()), .01)
gen x2 = round(invnormal(runiform()), .01)
gen y = runiform() < invlogit(.3 + .8*x1 - .5*x2)
gen st = mod(_n, 4)
gen psu = ceil(_n/5)
gen wt = 1 + mod(_n, 3)

* --- ml evaluators of the lf0, lf1, lf2 and gf0 types ---
program lf0logit
    args todo b lnfj
    tempvar xb
    mleval `xb' = `b'
    quietly replace `lnfj' = ln(invlogit(`xb')) if $ML_y1 == 1
    quietly replace `lnfj' = ln(invlogit(-`xb')) if $ML_y1 == 0
end
program lf1logit
    args todo b lnfj g1
    tempvar xb
    mleval `xb' = `b'
    quietly replace `lnfj' = ln(invlogit(`xb')) if $ML_y1 == 1
    quietly replace `lnfj' = ln(invlogit(-`xb')) if $ML_y1 == 0
    if `todo' == 0 exit
    quietly replace `g1' = $ML_y1 - invlogit(`xb')
end
program lf2logit
    args todo b lnfj g1 H
    tempvar xb
    mleval `xb' = `b'
    quietly replace `lnfj' = ln(invlogit(`xb')) if $ML_y1 == 1
    quietly replace `lnfj' = ln(invlogit(-`xb')) if $ML_y1 == 0
    if `todo' == 0 exit
    quietly replace `g1' = $ML_y1 - invlogit(`xb')
    if `todo' == 1 exit
    tempvar w
    quietly gen double `w' = invlogit(`xb')*invlogit(-`xb')
    mlmatsum `lnfj' `H' = -`w'
end
capture noisily ml model lf0 lf0logit (y = x1 x2)
capture noisily ml maximize
capture noisily ereturn list
capture noisily ml model lf1 lf1logit (y = x1 x2)
capture noisily ml maximize
capture noisily ml model lf1 lf1logit (y = x1 x2), vce(robust)
capture noisily ml maximize
capture noisily ml model lf2 lf2logit (y = x1 x2)
capture noisily ml maximize
capture noisily ml model gf0 lf0logit (y = x1 x2)
capture noisily ml maximize
* constraints()
constraint define 1 x1 = -x2
capture noisily ml model lf2 lf2logit (y = x1 x2), constraints(1)
capture noisily ml maximize
capture noisily matrix list e(Cns)
capture noisily ereturn list
capture noisily ml model d0 lf0logit (y = x1 x2), constraints(1)
* svy
svyset psu [pw=wt], strata(st)
capture noisily ml model lf1 lf1logit (y = x1 x2), svy
capture noisily ml maximize
capture noisily ereturn list
capture noisily svy: logit y x1 x2

* --- boxcox ---
clear
set obs 60
gen x = 1 + mod(_n*7, 13)
gen z = mod(_n, 3)
gen y = round(exp(.5 + .1*x + .2*z + mod(_n*11,7)/20), .01)
capture noisily boxcox y x z, model(rhsonly) notrans(z)
capture noisily ereturn list
capture noisily boxcox y x, model(lambda)
capture noisily boxcox y x z, model(lambda) notrans(z) nolog
capture noisily boxcox y x z, model(theta) notrans(z)
capture noisily ereturn list
capture noisily boxcox y x, notrans(z) nolog
capture noisily boxcox y z, model(rhsonly) nolog

* --- icc and loneway ---
clear
set obs 40
gen t = ceil(_n/4)
gen j = mod(_n-1,4)+1
gen y = t + mod(_n*7,5)*.5 + j*.3
capture noisily icc y t j, testvalue(.5)
capture noisily return list
capture noisily icc y t j, mixed consistency testvalue(.3)
capture noisily icc y t j, mixed absolute testvalue(.3)
capture noisily icc y t, testvalue(.2)
capture noisily return list
capture noisily loneway y t, exact
capture noisily return list
capture noisily loneway y t, mean
capture noisily loneway y t, median
drop if _n == 3 | _n == 17
capture noisily icc y t
capture noisily return list
capture noisily icc y t, testvalue(.4)
capture noisily icc y t j
capture noisily loneway y t, exact

* --- esize, pwmean, ameans ---
capture noisily esizei 20 5 2 25 4 3, unequal
capture noisily return list
capture noisily esizei 20 5 2 25 4 3, welch all
capture noisily esize twosample y if j < 3, by(j) unequal all
capture noisily pwmean y, over(j) unequal
capture noisily pwmean y, over(j) unequal effects mcompare(bonferroni)
capture noisily pwmean y, over(j) unequal cieffects mcompare(scheffe)
gen w = 1 + mod(_n, 3)
capture noisily ameans y [fw=w]
capture noisily return list
capture noisily ameans y [aw=w]

* --- lv and stem ---
clear
set obs 60
gen y = round(exp(.5 + .1*mod(_n*7, 13) + mod(_n*11,7)/20), .01)
gen x = 1 + mod(_n*7, 13)
capture noisily lv y, tail(8)
capture noisily lv y, tail(4) generate
capture noisily list _mid _psigma in 1/5
capture noisily lv y, generate
capture noisily lv x, tail(3)
capture noisily stem x, width(5)
capture noisily stem x, width(2)
capture noisily stem x, width(10)
capture noisily stem y, lines(2)

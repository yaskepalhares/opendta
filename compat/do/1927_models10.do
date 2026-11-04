* Models, part 10 (block 20h): eteffects with probit, fractional and exponential outcome models
* (the linear model is checked in case 1908). No reference until the Stata run. Every command
* whose behavior the manual does not settle runs under capture noisily, so that one Stata run
* shows the whole case.
clear
set obs 2000
set seed 1927
gen x = round(invnormal(runiform()), .01)
gen z = round(invnormal(runiform()), .01)
gen u = round(invnormal(runiform()), .01)
gen t = (.2 + .7*z + .5*x + u) > 0
gen yl = round(1 + x + 1.5*t + .6*u + .7*invnormal(runiform()), .01)
gen yb = (-.2 + .5*x + .8*t + .5*u + invnormal(runiform()) > 0)
gen yf = invlogit(-.2 + .5*x + .8*t + .5*u + invnormal(runiform()))
gen yc = round(exp(.2 + .3*x + .5*t + .3*u)*runiform()*2)
capture noisily eteffects (yl x) (t x z), atet
capture noisily eteffects (yb x, probit) (t x z)
capture noisily ereturn list
capture noisily eteffects (yb x, probit) (t x z), pomeans
capture noisily eteffects (yf x, fractional) (t x z)
capture noisily eteffects (yc x, exponential) (t x z)
capture noisily eteffects (yc x, exponential) (t x z), atet

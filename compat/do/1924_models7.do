* Models, part 7 (block 20e): slogit with dimension() above 1; exlogistic with binomial(),
* group(), condvars() and estconstant; expoisson with exposure(), offset() and condvars(). No
* reference until the Stata run. Every command whose behavior the manual does not settle runs
* under capture noisily, so that one Stata run shows the whole case.
clear
set obs 1500
set seed 1924
gen x1 = round(invnormal(runiform()), .01)
gen x2 = round(invnormal(runiform()), .01)
gen x3 = round(invnormal(runiform()), .01)
gen double u1 = runiform()
gen double e1 = .5*x1 + .3*x2
gen double e2 = .4*x2 - .6*x3
gen double s1 = exp(.5 - e1)
gen double s2 = exp(.2 - e2)
gen double s3 = exp(-.1 - .5*e1 - .5*e2)
gen double tot = 1 + s1 + s2 + s3
gen y = 4
replace y = 1 if u1 < s1/tot
replace y = 2 if u1 >= s1/tot & u1 < (s1+s2)/tot
replace y = 3 if u1 >= (s1+s2)/tot & u1 < (s1+s2+s3)/tot
capture noisily slogit y x1 x2 x3, dim(2)
capture noisily ereturn list
capture noisily matrix list e(b)
capture noisily slogit y x1 x2 x3, dim(3)
capture noisily slogit y x1 x2 x3, dim(2) baseoutcome(1)

clear
set obs 40
set seed 1924
gen g = ceil(_n/8)
gen x1 = round(runiform())
gen x2 = round(runiform()*2)
gen n = 1 + round(runiform()*3)
gen y = (runiform() < .3 + .3*x1)
gen yb = round(n*(.2 + .4*x1*runiform()))
gen t = 1 + round(runiform()*2)
gen c = round(t*exp(.2*x1)*runiform()*2)
capture noisily exlogistic y x1, condvars(x2)
capture noisily exlogistic y x1 x2, estconstant
capture noisily exlogistic yb x1 x2, binomial(n)
capture noisily ereturn list
capture noisily exlogistic y x1 x2, group(g)
capture noisily expoisson c x1 x2, exposure(t)
capture noisily ereturn list
gen lt = ln(t)
capture noisily expoisson c x1 x2, offset(lt)
capture noisily expoisson c x1, condvars(x2) offset(x2)

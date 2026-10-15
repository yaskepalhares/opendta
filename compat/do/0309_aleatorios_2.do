* Phase 3d: data to reproduce rnormal() and Stata's other
* distributions. The runiform() draws with the same seed show the available uniforms;
* by comparing, one can find out how many each function consumes and how.
clear
set obs 8
set seed 123
gen double u = runiform()
set seed 123
gen double z = rnormal()
set seed 123
gen double k = runiformint(1, 6)
set seed 123
gen double b = rbinomial(10, .5)
set seed 123
gen double p = rpoisson(3)
set seed 123
gen double c = rchi2(2)
set seed 123
gen double e = rexponential(1)
set seed 123
gen double t = rt(5)
format u-t %20.17f
list u z k, noobs
list b p c, noobs
list e t, noobs
set seed 123
display %20.17f rnormal() _n %20.17f runiform()
set seed 123
display %20.17f rnormal() _n %20.17f rnormal() _n %20.17f runiform()

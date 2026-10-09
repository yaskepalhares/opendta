* Fase 5b: poisson (offset, exposure, irr, vce) e predict.
clear
estimates clear
ereturn clear
set seed 707
set obs 90
gen x1 = round(runiform()*10, .1)
gen x2 = round(runiform()*5, .01)
gen t = 1 + mod(_n, 4)
gen c = floor(-ln(runiform())*exp(.1 + .1*x1))
poisson c x1 x2
ereturn list
poisson c x1 x2, irr
poisson c x1 x2, exposure(t)
poisson c x1 x2, offset(x2)
poisson c x1 x2, vce(robust)
predict n
predict ir, ir
predict xb, xb
summarize n ir xb
estimates store p1
poisson c x1
estimates store p2
estimates table p1 p2, stats(N ll)

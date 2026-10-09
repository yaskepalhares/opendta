* Fase 5a: predict, test, testparm, lincom e estimates.
clear
estimates clear
ereturn clear
set seed 31
set obs 40
gen g = 1 + mod(_n, 3)
gen x1 = round(runiform()*10, .1)
gen x2 = round(runiform()*5, .01)
gen y = 1 + x1 - 2*x2 + g + round(runiform()*4, .01)
replace x1 = . in 3
capture noisily predict p0
regress y x1 x2 i.g
predict yhat
predict double xb2, xb
predict r, residuals
predict sp, stdp
predict sf, stdf
predict hv, hat
predict rs, rstandard
predict rt, rstudent
predict cd, cooksd
predict rin if g == 1, residuals
summarize yhat xb2 r sp sf hv rs rt cd rin
list y yhat r sp hv cd in 1/5
test x1
test x1 x2
test x1 = x2
test x1 + x2 = -1
test (x1 = 0) (x2 = -2)
test 2.g 3.g
test 2.g = 3.g
test x1, notest
test x2, accumulate
test _b[x1] = 1
return list
testparm x1 x2
testparm i.g, equal
lincom x1 + x2
lincom x1 - 2*x2
lincom 3.g - 2.g
lincom x1, level(90)
return list
estimates store completo
regress y x1
estimates store simples
estimates dir
estimates table completo simples
estimates table completo simples, se
estimates table completo simples, stats(N r2) b(%9.3f)
estimates table completo simples, star
estimates restore completo
display e(cmd) " " e(N)
estimates replay simples
estimates drop simples
estimates dir
capture noisily estimates restore simples
capture noisily test x9

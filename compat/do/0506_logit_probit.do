* Fase 5b: logit, logistic e probit (log de iterações, tabelas, pós-estimação).
clear
estimates clear
ereturn clear
set seed 611
set obs 80
gen x1 = round(runiform()*10, .1)
gen x2 = round(runiform()*5, .01)
gen g = 1 + mod(_n, 3)
gen u = runiform()
gen y = u < invlogit(-2 + .3*x1 - .2*x2 + .5*(g==2))
gen w = 1 + mod(_n, 2)
logit y x1 x2
ereturn list
logit
logit y x1 x2 i.g
logit y x1 x2, or
logit y x1 x2, nolog level(90)
logit y x1 x2, vce(robust)
logit y x1 x2, vce(cluster g)
logit y x1 x2 [fweight=w]
logit y x1 x2 [pweight=w]
predict p
predict xb, xb
predict sp, stdp
summarize p xb sp
test x1 x2
lincom x1 - x2
logistic y x1 x2
logistic
logistic y x1 i.g, coef
probit y x1 x2
probit y x1 x2 i.g, vce(robust)
predict pp
summarize pp
gen y0 = 0
capture noisily logit y0 x1

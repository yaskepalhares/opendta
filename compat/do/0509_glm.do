* Fase 5b: glm (famílias binomial, poisson e gaussian; eform).
clear
estimates clear
ereturn clear
set seed 909
set obs 120
gen x1 = round(runiform()*10, .1)
gen x2 = round(runiform()*5, .01)
gen c = floor(-ln(runiform())*exp(.1 + .1*x1) * (0.5 + runiform()))
gen y = runiform() < invlogit(-2 + .3*x1 - .2*x2)
gen yr = runiform() < .04*exp(.15*x1)
gen yc = 2 + x1 - x2 + round(runiform()*3, .01)
glm y x1 x2, family(binomial)
ereturn list
glm yr x1 x2, family(binomial) link(log) eform
glm c x1 x2, family(poisson)
glm c x1 x2, family(poisson) vce(robust)
glm c x1 x2, family(poisson) eform
glm yc x1 x2
ereturn list
predict mu
summarize mu

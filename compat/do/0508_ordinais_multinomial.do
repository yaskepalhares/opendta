* Fase 5b: ologit, oprobit e mlogit.
clear
estimates clear
ereturn clear
set seed 811
set obs 120
gen x1 = round(runiform()*10, .1)
gen x2 = round(runiform()*5, .01)
gen u = runiform()
gen r = 1 + (u < invlogit(-2+.3*x1)) + (u < invlogit(-4+.3*x1)) + (runiform()<.2)
label define rl 1 "baixo" 2 "medio" 3 "alto" 4 "muito alto"
label values r rl
ologit r x1 x2
ereturn list
ologit r x1 x2, vce(robust)
oprobit r x1 x2
mlogit r x1 x2
ereturn list
mlogit r x1, baseoutcome(1)
mlogit r x1, rrr
test x1

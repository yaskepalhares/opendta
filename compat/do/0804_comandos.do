* Fase 8: lowess, lrtest, hausman, xtsum, xtdescribe, xtpoisson, xtlogit,
* qreg e putexcel
clear
set seed 804
set obs 200
gen id = ceil(_n/10)
gen t = mod(_n-1, 10) + 1
xtset id t
gen x = round(runiform()*4, .01)
gen u = (id - 10)/10
gen yc = floor(exp(.2 + .3*x + u) * runiform()*2)
gen yb = runiform() < invlogit(-1 + .6*x + u)
gen y = 1 + 2*x + u + round(runiform(), .01)

lowess y x, generate(lw) nograph
summarize lw
lowess y x, bwidth(.3) mean generate(lw2) nograph
summarize lw2
xtsum y x
xtdescribe
xtpoisson yc x, fe
xtpoisson yc x, re
xtlogit yb x, fe
xtlogit yb x, re
quietly xtreg y x, fe
estimates store fe
quietly xtreg y x, re
estimates store re
hausman fe re
return list
quietly logit yb x
estimates store small
quietly logit yb x u
estimates store big
lrtest big small
return list
qreg y x
qreg y x, quantile(.25)
ereturn list
quietly regress y x
matrix b = e(b)
putexcel set odta_res.xlsx, replace
putexcel A1 = "Resultados" B2 = (e(r2)) A3 = matrix(b), names
erase odta_res.xlsx

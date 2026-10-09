* Fase 5b: margins e nlcom.
clear
estimates clear
ereturn clear
set seed 611
set obs 80
gen x1 = round(runiform()*10, .1)
gen x2 = round(runiform()*5, .01)
gen g = 1 + mod(_n, 3)
gen y = runiform() < invlogit(-2 + .3*x1 - .2*x2 + .5*(g==2))
gen yc = 1 + x1 - x2 + g + round(runiform(), .01)
logit y x1 x2 i.g
margins
margins g
margins, dydx(x1)
margins, dydx(*)
margins, at(x1=(2 4 6))
margins, atmeans
margins g, atmeans
return list
regress yc x1 x2 i.g
margins g
margins, dydx(x1 x2)
margins, at(x1=(2(2)6))
nlcom _b[x1]/_b[x2]
nlcom (r1: _b[x1]/_b[x2]) (r2: exp(_b[x1]))

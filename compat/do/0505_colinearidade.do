* Fase 5a: qual variável o regress omite por colinearidade, pesos e opções.
clear
estimates clear
ereturn clear
set seed 505
set obs 50
gen x1 = round(runiform()*10, .1)
gen x2 = round(runiform()*5, .01)
gen x3 = 2*x1
gen x4 = x1 + x2
gen d1 = _n <= 10
gen d2 = _n > 10 & _n <= 30
gen d3 = _n > 30
gen g = 1 + d2 + 2*d3
gen y = 1 + x1 - x2 + 2*d2 + round(runiform()*3, .01)
regress y x3 x1 x2
regress y x1 x3 x2
regress y x1 x2 x4
regress y x4 x1 x2
regress y d1 d2 d3
regress y d3 d2 d1
regress y d1 d2 d3, noconstant
regress y i.g d2
regress y x1 if d1
gen um = 1
regress y x1 um
regress y x1 x2 [aweight=g], vce(robust)
regress y x1 x2 [iweight=g]
regress y x1 x2, vce(hc2)
regress y x1 x2, vce(hc3)
regress y x1 x2, vce(cluster g)
regress y x1 x2, level(99) noheader
regress, level(80)

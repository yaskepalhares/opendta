* Fase 5a: tsset, xtset e operadores de séries temporais.
clear
estimates clear
set obs 15
gen t = _n
gen y = 10 + 2*t + mod(t*7, 5)
gen x = mod(t*3, 7)
capture noisily gen ly = L.y
tsset t
tsset
return list
gen ly = L.y
gen l2y = L2.y
gen fy = F.y
gen dy = D.y
gen d2y = D2.y
gen sy = S.y
gen lfy = LF.y
list t y ly l2y fy dy d2y sy lfy in 1/5
regress y L.y x
regress y L(0/2).x
regress D.y L.x
regress y L.y L.x LD.x
tsset, clear
clear
set obs 16
gen t = _n
replace t = t + 2 if _n > 8
gen y = mod(t*5, 9)
tsset t
gen ly = L.y
list t y ly in 6/11
clear
set obs 12
gen id = ceil(_n/4)
gen t = mod(_n-1,4)+1
gen y = id*3 + t + mod(_n,3)
xtset id t
xtset
return list
gen ly = L.y
list in 1/6
regress y L.y i.id
drop in 5
xtset id t
clear
set obs 6
gen d = mdy(1, 1, 2000) + _n - 1
format d %td
gen y = _n^2
tsset d
gen dy = D.y
list
clear
set obs 4
gen t = _n
replace t = 2 in 3
capture noisily tsset t

* Fase 5a: variáveis fatoriais e interações em regress.
clear
estimates clear
set seed 77
set obs 48
gen g = 1 + mod(_n, 3)
gen h = mod(_n, 2)
gen k = 1 + mod(_n, 4)
gen x = round(runiform()*10, .1)
gen y = 2 + .5*x + g + 2*h - h*g + round(runiform()*3, .01)
label define glab 1 "baixo" 2 "medio" 3 "alto"
label values g glab
regress y i.g
regress y i.k x
regress y ib2.k x
regress y ib(last).k
regress y ibn.g, noconstant
regress y i.g##i.h
testparm i.g
testparm i.g#i.h
regress y g#h
regress y i.h##c.x
regress y i.g#c.x
regress y c.x##c.x
regress y i.(g h) x
regress y i(1 3).k
regress y i.k
testparm i.k
matrix list e(b)
regress y i.g, noconstant
gen gneg = g - 2
capture noisily regress y i.gneg
gen gfrac = g/2
capture noisily regress y i.gfrac

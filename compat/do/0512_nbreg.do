* Fase 5b: nbreg (Poisson, modelo só com a constante, modelo completo).
clear
estimates clear
ereturn clear
set seed 909
set obs 120
gen x1 = round(runiform()*10, .1)
gen x2 = round(runiform()*5, .01)
gen c = floor(-ln(runiform())*exp(.1 + .1*x1) * (0.5 + runiform()))
nbreg c x1 x2
ereturn list
nbreg c x1 x2, irr
nbreg c x1 x2, vce(robust) nolog

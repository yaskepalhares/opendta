* Fase 6: svyset, svydescribe, mean/proportion/total/ratio e svy:.
clear
estimates clear
ereturn clear
set seed 600
set obs 200
gen estrato = ceil(_n/50)
gen upa = ceil(_n/10)
gen peso = 50 + mod(_n, 7)*10
gen sexo = mod(_n, 2)
label define sexol 0 "masc" 1 "fem"
label values sexo sexol
gen idade = 18 + floor(runiform()*60)
gen ativo = runiform() < invlogit(-1 + .02*idade - .5*sexo)
gen imc = 22 + .05*idade + 2*sexo + round(runiform()*4, .01)
mean imc
mean imc idade, over(sexo)
mean imc [pweight=peso]
proportion ativo
total ativo
ratio imc/idade
capture noisily svy: mean imc
svyset upa [pweight=peso], strata(estrato)
svyset
svydescribe
svy: mean imc
ereturn list
svy: mean imc, over(sexo)
svy: proportion ativo
svy: total ativo
svy: ratio imc/idade
svy, subpop(sexo): mean imc
svy: regress imc idade i.sexo
svy: logistic ativo idade i.sexo
svy: logit ativo idade
svy: tabulate ativo
svy: tabulate ativo sexo
svyset _n [pweight=peso]
svy: mean imc
gen npop = 20
svyset upa [pweight=peso], strata(estrato) fpc(npop)
svy: mean imc
gen tx = .25
svyset upa [pweight=peso], strata(estrato) fpc(tx)
svy: mean imc

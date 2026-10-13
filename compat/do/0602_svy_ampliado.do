* Fase 6b: svy com vários estágios, réplicas (jackknife, BRR, bootstrap, SDR), estat effects e mais modelos
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
gen nivel = 1 + (imc > 26) + (imc > 28)
gen cont = floor(idade/10)
gen upa2 = mod(_n, 2)
gen npop = 20
gen n2 = 10
svyset upa [pweight=peso], strata(estrato)
svy: mean imc
estat effects
svy: mean imc idade, over(sexo)
estat effects
svy: logit ativo idade
estat effects
svy, vce(jackknife): mean imc
svy, vce(jackknife) mse: mean imc
svy, vce(jackknife): regress imc idade i.sexo
svy: ologit nivel idade i.sexo
svy: oprobit nivel idade
svy: mlogit nivel idade
svy: mlogit nivel idade, rrr
svy: nbreg cont imc
svy: glm ativo idade, family(binomial) eform
svy, subpop(sexo): ologit nivel idade
svyset upa2 [pweight=peso], strata(estrato) vce(brr)
svy: mean imc
svy: proportion ativo
svyset upa [pweight=peso], strata(estrato) fpc(npop) || _n, fpc(n2)
svy: mean imc
forvalues r = 1/10 {
    gen bw`r' = peso*(1 + .1*sin(_n*`r'))
}
svyset [pweight=peso], bsrweight(bw*) vce(bootstrap)
svy: mean imc
svyset [pweight=peso], sdrweight(bw*) vce(sdr)
svy: mean imc
svyset [pweight=peso], jkrweight(bw*) vce(jackknife)
svy: mean imc

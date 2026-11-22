* Probe (survey): factor notation and equation names in post-estimation after svy: poisson/logit/regress/glm.
* Each block runs under capture noisily so that one run shows every command.
version 14
set more off
clear
set seed 123
set obs 600
gen psu = ceil(_n/20)
gen strat = mod(psu,3)+1
gen w = 1+round(runiform()*3,.1)
gen ano = 2006 + floor(runiform()*10)
gen sexo = runiform()<.5
gen idade = 18+floor(runiform()*60)
gen faixa = 1+floor(runiform()*3)
gen pc = 1+floor(runiform()*10)+sexo
gen yb = runiform() < .3+.1*sexo
svyset psu [pw=w], strata(strat)

display "=== svy: poisson, naming forms"
capture noisily svy: poisson pc ano i.sexo ib2.faixa c.idade
foreach e in "ano" "c.ano" "1.sexo" "3.faixa" "_b[ano]" "[pc]ano" "[pc]_b[ano]" "[pc]c.idade" "ib2.faixa" "ano + [pc]c.idade" {
    display "--- lincom `e'"
    capture noisily lincom `e'
}
display "--- lincom eform / irr then return list"
capture noisily lincom ano, eform
return list
capture noisily lincom ano, irr
return list
foreach e in "ano" "[pc]ano" "c.ano" "_b[ano] = 0" "[pc]_b[ano] = 0" "3.faixa = 1.faixa" "1.sexo = c.idade" {
    display "--- test `e'"
    capture noisily test `e'
}
capture noisily test [pc]: ano
capture noisily testparm i.sexo
capture noisily testparm c.ano
capture noisily testparm i.faixa
capture noisily testparm i(1 3).faixa
capture noisily nlcom exp(_b[ano])
capture noisily nlcom [pc]_b[ano]
capture noisily nlcom (r: exp(_b[ano]))
foreach e in "_b[ano]" "_se[ano]" "_b[pc:ano]" "[pc]_b[ano]" "[pc]ano" "_b[1.sexo]" "_b[3.faixa]" "_b[2b.faixa]" "_b[c.idade]" {
    display "--- display `e'"
    capture noisily display `e'
}
display "--- margins/contrast/estat"
capture noisily margins sexo
capture noisily margins r.sexo
capture noisily margins ar.faixa
capture noisily margins faixa, at(ano=(2006 2014))
capture noisily margins, dydx(c.ano)
capture noisily margins, eydx(ano)
capture noisily margins, expression(exp(xb()))
capture noisily contrast faixa
capture noisily pwcompare faixa
capture noisily predictnl p = exp(xb()), se(sep)
foreach s in vce effects svyset sd strata size cv summarize {
    display "--- estat `s'"
    capture noisily estat `s'
}
display "=== svy with subpop and interaction"
capture noisily svy, subpop(if sexo==1): poisson pc ano c.idade
capture noisily lincom ano, irr
return list
capture noisily svy: poisson pc ano i.sexo##c.idade, irr
capture noisily lincom 1.sexo + 30*1.sexo#c.idade
capture noisily margins sexo#c.idade
display "=== svy: logit, or"
capture noisily svy: logit yb ano i.sexo i.faixa, or
capture noisily lincom ano, or
return list
capture noisily lincom [yb]ano, or
capture noisily lincom 2.faixa - 3.faixa, or
capture noisily test [yb]ano
display "=== svy: regress"
capture noisily svy: regress pc ano i.sexo
capture noisily lincom ano
capture noisily lincom [pc]ano
display "=== glm and svy: glm"
capture noisily glm pc ano i.sexo, family(poisson) link(log) vce(robust) eform
return list
capture noisily lincom ano, eform
capture noisily svy: glm pc ano i.sexo, family(poisson) link(log)
capture noisily lincom [pc]ano
capture noisily lincom ano, irr
ereturn list
display "=== replicate-weight prefixes"
capture noisily svy bootstrap: poisson pc ano
capture noisily svy jackknife: poisson pc ano

* Fase 8: oneway, anova, tabulate summarize() e estat
clear
set obs 74
gen mpg = 12 + mod(_n*7, 30) + (mod(_n, 4) == 0)*mod(_n, 9)
gen price = 3000 + mod(_n*131, 9000) + mod(_n, 7)^3
gen foreign = _n > 52
label define orig 0 "Domestic" 1 "Foreign"
label values foreign orig
gen rep = mod(_n, 5) + 1
gen peso = 2000 + 20*mpg + mod(_n*17, 300)
label variable peso "Peso do carro"

oneway mpg rep
return list
oneway mpg rep, tabulate
oneway price foreign, tabulate means
tabulate foreign, summarize(mpg)
tabulate rep, summarize(price) means
anova mpg rep foreign
ereturn list
anova mpg rep##foreign
anova price rep c.peso
regress price mpg peso
estat ic
return list
estat vce
estat vce, correlation
estat summarize
estat hettest
return list
estat hettest, rhs
estat hettest, iid
estat ovtest
return list
logit foreign mpg peso
estat classification
return list
estat gof
return list
estat gof, group(5)
estat ic

* Phase 8: nonparametric tests, normality, sdtest and immediate commands
clear
set obs 74
gen mpg = 12 + mod(_n*7, 30) + (mod(_n, 4) == 0)*mod(_n, 9)
gen price = 3000 + mod(_n*131, 9000) + mod(_n, 7)^3
gen foreign = _n > 52
label define orig 0 "Domestic" 1 "Foreign"
label values foreign orig
gen rep = mod(_n, 5) + 1
gen peso = 2000 + 20*mpg + mod(_n*17, 300)

ranksum mpg, by(foreign)
return list
ranksum price, by(foreign) porder
signrank mpg = 25
return list
signrank mpg = peso
kwallis mpg, by(rep)
return list
spearman mpg price
return list
spearman mpg price peso
swilk mpg price peso
return list
sktest mpg price peso
sktest mpg, noadjust
sdtest mpg, by(foreign)
return list
sdtest mpg == 8
sdtest mpg == price
* Stata 14.0 syntax (cii means/proportions only exists from 14.1 on)
cii 166 19509 4379
return list
cii 100 30
cii 100 30, wilson
cii 100 30, agresti
cii 100 30, jeffreys
cii 100 30, wald
cii 50 7, poisson
ttesti 24 62.6 15.8 75
ttesti 20 20 5 32 15 4
ttesti 20 20 5 32 15 4, unequal
prtesti 50 .52 .7
tabi 30 18 \ 38 14, chi2 exact
tabi 10 20 30 \ 5 15 25, row

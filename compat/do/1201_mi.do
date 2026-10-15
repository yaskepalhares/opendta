* Phase 10: multiple imputation (mi set, register, impute, estimate, passive, xeq)
clear
set seed 1201
set obs 200
gen age = 20 + floor(50*runiform())
gen female = runiform() < .5
gen bmi = 22 + .08*age + 1.5*female + rnormal()*3
gen smoke = runiform() < invlogit(-1 + .02*age)
gen y = 5 + .3*bmi - .5*female + rnormal()*2
replace bmi = . if runiform() < .2
replace smoke = . if runiform() < .15
mi set mlong
mi register imputed bmi smoke
mi describe
mi impute regress bmi = age female y, add(5) rseed(1234)
mi describe
mi estimate: regress y bmi female
mi impute chained (regress) bmi (logit) smoke = age female y, replace rseed(99)
mi estimate: regress y bmi female smoke
mi estimate: logit smoke bmi age
mi passive: generate bmi2 = bmi^2
mi xeq 0 1: summarize bmi
mi extract 2, clear
describe, short
count

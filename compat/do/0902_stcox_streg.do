* Phase 10: stcox, estat phtest, predict after stcox and streg
clear
set seed 902
set obs 120
gen id = _n
gen drug = 1 + (runiform() < .5)
gen age = 40 + floor(20*runiform())
gen sex = runiform() < .5
gen double tt = ceil(-ln(runiform())/(0.03*exp(0.7*(drug==1) + 0.04*(age-50))))
gen died = runiform() < .8
replace tt = 30 if tt > 30
replace died = 0 if tt == 30
stset tt, failure(died)

stcox drug age
ereturn list
stcox drug age, nohr
stcox drug age, efron
stcox drug age i.sex, vce(robust)
stcox drug age, strata(sex)
stcox drug age
estat phtest
estat phtest, detail
estat phtest, log detail
predict hr
predict xb, xb
predict s0, basesurv
predict H0, basechazard
predict cs, csnell
predict mg, mgale
summarize hr xb s0 H0 cs mg

streg drug age, distribution(exponential)
streg drug age, distribution(exponential) time
streg drug age, distribution(weibull)
ereturn list
streg drug age, distribution(weibull) time
streg drug age, distribution(weibull) nohr
streg drug age, distribution(gompertz)
streg drug age, distribution(lognormal)
streg drug age, distribution(loglogistic)
streg drug age, distribution(weibull) vce(robust)
predict med, median time
predict cs2, csnell
summarize med cs2

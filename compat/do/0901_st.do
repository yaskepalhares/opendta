* Phase 10: survival-time data: stset, stdescribe, stsum, sts list,
* sts test and sts generate
clear
set seed 901
set obs 60
gen id = _n
gen drug = 1 + (runiform() < .5)
gen age = 40 + floor(20*runiform())
gen double tt = ceil(-ln(runiform())/(0.02*exp(0.6*(drug==1) + 0.03*(age-50))))
gen died = runiform() < .75
replace tt = 40 if tt > 40
replace died = 0 if tt == 40
replace tt = . in 60

stset tt, failure(died)
stdescribe
stsum
stsum, by(drug)
sts list in 1/20
sts list, by(drug) at(5 10 20 30)
sts list, failure at(10 20)
sts test drug
sts test drug, wilcoxon
sts test drug, strata(age)
sts generate s = s se = se(s) lb = lb(s) ub = ub(s) na = na
summarize s se lb ub na
st
stset tt, failure(died==1) scale(7)
stsum
stset tt, failure(died) origin(time 2) enter(time 5)
stdescribe
stset, clear

clear
input id t0 t1 died x
1 0 5 0 1
1 5 9 1 0
2 0 4 0 1
2 4 12 0 1
3 2 6 1 0
4 0 3 0 0
4 3 3 1 0
5 0 . 1 1
6 1 8 1 1
7 0 10 1 0
7 10 14 1 0
end
stset t1, id(id) failure(died)
stset t1, id(id) failure(died) time0(t0)
stdescribe
sts list
stset t1, id(id) failure(died) time0(t0) exit(time .)
stsum

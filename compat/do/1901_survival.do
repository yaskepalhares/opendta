* Survival analysis, the remaining commands: strate, stptime, stci, stir, stmh, stmc, ltable,
* stsplit, stjoin, stfill, stgen, stvary, stbase, snapspan, ctset, ct, cttost, sttoct, sttocc,
* stcrreg, stcurve, streg with distribution(ggamma) and frailty()/shared(), stcox with the exactm
* and exactp ties. Every command whose output the manual does not settle runs under capture
* noisily, so that one Stata run shows the whole case. The data are built without random numbers.
clear
set obs 60
gen id = _n
gen g = ceil(_n/4)
gen x = mod(_n*7,11)/5 - 1
gen byte grp = 1 + mod(_n,3)
gen byte trt = mod(_n,2)
gen double t = 1 + mod(_n*13,23) + mod(_n,4)/2
gen byte ev = mod(_n,5) != 0
replace ev = 2 if ev == 1 & mod(_n,7) == 0
gen double pop = 10 + mod(_n,5)
stset t, failure(ev==1) id(id)

* rates, person-time, intervals
capture noisily strate
capture noisily strate grp
capture noisily strate grp, per(100)
capture noisily strate grp trt, jackknife
capture noisily strate grp, output(odta_rates, replace)
capture noisily describe using odta_rates
capture noisily stptime
capture noisily stptime, by(grp) per(100) title(Person-years)
capture noisily stptime, by(grp) at(5 10 20)
capture noisily display r(rate)
capture noisily stci
capture noisily stci, by(grp)
capture noisily stci, p(25) by(grp)
capture noisily stci, rmean
capture noisily stci, emean
capture noisily stir trt
capture noisily stmh trt grp
capture noisily stmc trt grp
capture noisily ltable t ev
capture noisily ltable t ev, intervals(5)
capture noisily ltable t ev, by(grp) intervals(5) survival
capture noisily ltable t ev, intervals(5) failure
capture noisily ltable t ev, intervals(5) hazard
capture noisily ltable t ev, intervals(5) noadjust
capture noisily ltable t ev, by(grp) intervals(5) test notable

* manipulation of survival-time records
preserve
capture noisily stsplit per, at(5 10 20)
capture noisily list id _t0 _t _d per in 1/8
capture noisily stjoin
capture noisily count
restore
preserve
capture noisily stsplit, every(7)
capture noisily count
restore
preserve
capture noisily stsplit per, at(failures)
capture noisily count
restore
capture noisily stvary x trt
capture noisily stfill x, baseline
capture noisily stgen nfail = nfailures()
capture noisily stgen everf = ever(ev==1)
capture noisily summarize nfail everf
preserve
capture noisily stbase, replace
restore
preserve
capture noisily snapspan id t x, generate(t0) replace
restore

* count-time data
clear
input byte time byte nfail byte ncens byte nenter
0 0 0 11
1 1 0 0
2 2 1 0
3 1 1 0
4 2 0 0
5 0 3 0
end
capture noisily ctset time nfail ncens nenter
capture noisily ct
capture noisily cttost
capture noisily summarize _t0 _t _d
capture noisily stptime
capture noisily sttoct nf nc ne
capture noisily list

* case-control sampling (distinct failure times, every control is sampled)
clear
set obs 12
gen id = _n
gen double t = _n
gen byte d = inlist(_n,2,5,8)
gen byte sex = mod(_n,2)
stset t, failure(d) id(id)
preserve
capture noisily sttocc, number(20) nodots
capture noisily describe
capture noisily sort _set _case
capture noisily list _set _case _time id, sepby(_set)
restore
preserve
capture noisily sttocc, match(sex) number(20) nodots
capture noisily list _set _case _time id sex, sepby(_set)
restore

* competing risks
clear
set obs 60
gen id = _n
gen x = mod(_n*7,11)/5 - 1
gen byte grp = 1 + mod(_n,3)
gen double t = 1 + mod(_n*13,23) + mod(_n,4)/2
gen byte ev = mod(_n,5) != 0
replace ev = 2 if ev == 1 & mod(_n,7) == 0
stset t, failure(ev==1) id(id)
capture noisily stcrreg x, compete(ev==2)
capture noisily stcrreg x i.grp, compete(ev==2) noshr
capture noisily stcrreg x, compete(ev==2) nolog noheader
capture noisily ereturn list
capture noisily predict double xbhat, xb
capture noisily predict double shrhat, shr
capture noisily summarize xbhat shrhat
capture noisily stcurve, cif at(x=0) nodraw
capture noisily stcurve, cif at1(x=-1) at2(x=1) nodraw

* ties in stcox
capture noisily stcox x i.grp, exactm
capture noisily stcox x i.grp, exactp
capture noisily stcox x i.grp, exactm nolog
capture noisily stcox x, efron nolog
capture noisily stcox x, exactm vce(robust)

* streg
stset t, failure(ev==1) id(id)
capture noisily streg x, distribution(ggamma)
capture noisily streg x, distribution(ggamma) time
capture noisily streg x, distribution(ggamma) nolog nohr
capture noisily ereturn list
capture noisily predict double mu, xb
capture noisily streg x, distribution(weibull) frailty(gamma)
capture noisily streg x, distribution(weibull) frailty(invgaussian)
capture noisily gen g = ceil(id/4)
capture noisily streg x, distribution(weibull) frailty(gamma) shared(g)
capture noisily streg x, distribution(exponential) frailty(invgaussian) shared(g) time
capture noisily ereturn list
capture noisily streg x, distribution(weibull) frailty(gamma) nolog
capture noisily stcurve, survival at(x=0) nodraw
capture noisily stcurve, hazard at1(x=-1) at2(x=1) nodraw
capture noisily stcurve, cumhaz nodraw
capture noisily stcurve, survival alpha1 nodraw
capture noisily stcurve, survival unconditional nodraw

erase odta_rates.dta

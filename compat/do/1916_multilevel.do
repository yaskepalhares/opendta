* Multilevel models, second part: the me commands with nested levels of random effects,
* weights (observation level and the fweight()/pweight() options of the levels) and
* vce(robust)/vce(cluster); mestreg with nested levels, weights and predict median, mean and
* hazard; predict, reffects after mixed with R. equations. Every command whose behavior the
* manual does not settle runs under capture noisily, so that one Stata run shows the whole case.
clear
set seed 1916
set obs 30
gen s = _n
gen us = invnormal(runiform())*.8
expand 4
bysort s: gen c = _n
gen uc = invnormal(runiform())*.5
expand 8
sort s c, stable
gen x = invnormal(runiform())
gen yb = runiform() < invlogit(-.2 + .7*x + us + uc)
gen yg = 1 + .5*x + us + uc + invnormal(runiform())
gen double lam = exp(.3 + .5*x + .6*us + .4*uc)
gen double pu = runiform()
gen double pk = exp(-lam)
gen double cp = pk
gen yc = 0
forvalues k = 1/30 {
    quietly replace yc = yc + 1 if pu > cp
    quietly replace pk = pk*lam/`k'
    quietly replace cp = cp + pk
}
drop lam pu pk cp
gen fw = 1 + mod(_n, 3)
gen ws = 1 + mod(s, 2)
gen cl = ceil(s/5)

* nested levels
melogit yb x || s: || c:
capture noisily estat icc
capture noisily estat group
predict r1 r2, reffects
predict m
capture noisily summarize r1 r2 m
ereturn list
capture noisily meglm yg x || s: || c:
capture noisily mixed yg x || s: || c:
capture noisily mepoisson yc x || s: || c:
capture noisily mepoisson yc x || s: x || c:
capture noisily meprobit yb x || s: || c:, intpoints(9)
capture noisily melogit yb x || s: || c:, intmethod(ghermite)
capture noisily meqrlogit yb x || s: || c:

* vce(robust) and vce(cluster)
melogit yb x || s:, vce(robust)
ereturn list
capture noisily melogit yb x || s:, vce(cluster cl)
capture noisily melogit yb x || s:, vce(cluster c)
capture noisily mepoisson yc x || s: || c:, vce(robust)

* weights
capture noisily melogit yb x [fw=fw] || s:
capture noisily melogit yb x [pw=fw] || s:
capture noisily melogit yb x [iw=fw] || s:
capture noisily meglm yg x || s:, fweight(ws)
capture noisily meglm yg x || s:, pweight(ws)
capture noisily meglm yg x || s:, pweight(fw)
capture noisily melogit yb x [pw=fw] || s:, pweight(ws)
capture noisily melogit yb x [aw=fw] || s:

* mestreg
gen double tw = (-ln(runiform())/exp(-1 + .5*x + us + uc))^(1/1.4)
gen dw = tw < 2
replace tw = min(tw, 2)
stset tw, failure(dw)
mestreg x || s: || c:, distribution(weibull)
predict md, median
predict mn, mean
predict hz, hazard
summarize md mn hz
capture noisily mestreg x || s:, distribution(exponential) vce(robust)
capture noisily mestreg x || s:, distribution(lognormal)
capture noisily predict md2, median
capture noisily predict mn2, mean
capture noisily predict hz2, hazard
capture noisily summarize md2 mn2 hz2
capture noisily mestreg x || s:, distribution(loglogistic)
capture noisily predict md3, median
capture noisily predict hz3, hazard
capture noisily summarize md3 hz3
capture noisily mestreg x || s:, distribution(gamma)
capture noisily predict md4, median
capture noisily predict mn4, mean
capture noisily predict hz4, hazard
capture noisily summarize md4 mn4 hz4
stset tw [pw=fw], failure(dw)
capture noisily mestreg x || s:, distribution(weibull)

* predict, reffects after mixed with R. equations
clear
set seed 2907
set obs 200
gen a = ceil(runiform()*8)
gen b = ceil(runiform()*6)
gen x = invnormal(runiform())
gen double ua = .
gen double ub = .
forvalues j = 1/8 {
    quietly replace ua = `j'*.2 - .9 if a == `j'
}
forvalues j = 1/6 {
    quietly replace ub = `j'*.15 - .5 if b == `j'
}
gen y = 1 + .6*x + ua + ub + invnormal(runiform())
mixed y x || a: || _all: R.b
capture noisily predict rr, reffects
predict ra rb, reffects
summarize ra rb
mixed y x || _all: R.a || _all: R.b
predict sa sb, reffects
summarize sa sb

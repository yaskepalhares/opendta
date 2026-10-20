* Phase 16: meglm, mecloglog, menbreg, meologit, meoprobit, meqrlogit, meqrpoisson, mestreg and crossed effects in mixed
clear
set seed 1905
set obs 40
gen g = _n
gen u = invnormal(runiform())*.7
gen v = invnormal(runiform())*.4
expand 15
gen x = invnormal(runiform())
gen yg = .3 + .5*x + u + v*x + invnormal(runiform())*1.2
gen yb = runiform() < 1 - exp(-exp(-.2 + .6*x + u + v*x))
* counts and positive responses are built by inversion (the rpoisson()/rgamma() draws differ from Stata's)
gen double lam = exp(.3 + .5*x + u*.6 + v*x*.5)
gen double pu = runiform()
gen double pk = exp(-lam)
gen double cp = pk
gen yc = 0
forvalues k = 1/30 {
    quietly replace yc = yc + 1 if pu > cp
    quietly replace pk = pk*lam/`k'
    quietly replace cp = cp + pk
}
gen double ya = lam*exp(.5*invnormal(runiform()))
gen double lat = .8*x + u + v*x + logit(runiform())
gen yo = 1 + (lat > -.7) + (lat > .8)
gen yk = 0
forvalues k = 1/6 {
    quietly replace yk = yk + (runiform() < invlogit(-.3 + .7*x + u))
}
gen nt = 6
drop lam pu pk cp lat

* meglm: families and links
meglm yg x || g:
meglm yb x || g:, family(bernoulli) link(logit)
meglm yc x || g:, family(poisson) link(log)
capture noisily meglm yk x || g:, family(binomial 6) link(probit)
capture noisily meglm yk x || g:, family(binomial nt) link(logit)
capture noisily meglm ya x || g:, family(gamma) link(log)
capture noisily meglm yc x || g:, family(nbinomial mean) link(log)
capture noisily meglm yc x || g:, family(nbinomial constant) link(log)
capture noisily meglm yo x || g:, family(ordinal) link(logit)
capture noisily meglm yg x || g:, link(log)
capture noisily meglm yc x || g:, family(poisson) link(log) exposure(nt)
capture noisily meglm yc x || g:, family(poisson) link(log) intmethod(laplace)
capture noisily meglm yb x || g:, family(bernoulli) link(cloglog) intpoints(9)
capture noisily meglm yg x || g: x, covariance(unstructured) stddeviations
capture noisily meglm yb x || g:, family(bernoulli) link(logit) or
capture noisily meglm yc x || g:, family(poisson) irr
capture noisily estat icc
capture noisily estat icc
predict mu
predict xbv, xb
predict eta, eta
predict re, reffects
capture noisily predict mufe, mu conditional(fixedonly)
capture noisily predict mufe2, fixedonly
capture noisily summarize mu xbv eta re mufe mufe2
capture noisily predict mrs, mu conditional(ebmeans)
capture noisily estat recovariance
capture noisily estat group
ereturn list

* mecloglog, menbreg, meologit, meoprobit
mecloglog yb x || g:
capture noisily mecloglog yb x || g:, intpoints(9)
capture noisily mecloglog yb x || g:, intmethod(ghermite)
capture noisily mecloglog yb x || g: x, covariance(independent)
capture noisily estat icc
predict pc
summarize pc
ereturn list
menbreg yc x || g:
capture noisily menbreg yc x || g:, dispersion(constant)
capture noisily menbreg yc x || g:, irr
capture noisily menbreg yc x || g:, exposure(nt)
capture noisily estat icc
predict nb
summarize nb
ereturn list
meologit yo x || g:
capture noisily meologit yo x || g:, or
capture noisily meologit yo x || g:, intpoints(9)
capture noisily meologit yo x || g: x, covariance(unstructured)
capture noisily predict po1 po2 po3, pr
capture noisily predict ol, xb
capture noisily predict ore1 ore2, reffects
capture noisily predict opr, outcome(2)
capture noisily summarize po1 po2 po3 ol ore1 ore2
ereturn list
meoprobit yo x || g:
capture noisily meoprobit yo x || g:, intpoints(9)
capture noisily estat icc
ereturn list

* meqrlogit and meqrpoisson
meqrlogit yb x || g:
capture noisily meqrlogit yb x || g:, or
capture noisily meqrlogit yb x || g: x, covariance(unstructured)
capture noisily meqrlogit yb x || g:, intpoints(9)
capture noisily meqrlogit yb x || g:, laplace
capture noisily meqrlogit yb x || g:, refineopts(iterate(0))
capture noisily estat icc
predict qp
capture noisily predict qre, reffects
capture noisily predict qxb, xb
summarize qp
ereturn list
meqrpoisson yc x || g:
capture noisily meqrpoisson yc x || g:, irr
capture noisily meqrpoisson yc x || g:, exposure(nt)
capture noisily meqrpoisson yc x || g: x, covariance(independent)
predict qm
summarize qm
ereturn list

* mestreg
gen double tw = (-ln(runiform())/exp(-1 + .5*x + u + v*x))^(1/1.4)
gen double tl = exp(1 - .4*x + u + .7*invnormal(runiform()))
gen double tc = tl*3
gen dw = tw < 2
replace tw = min(tw, 2)
gen dl = tl < 6
replace tl = min(tl, 6)
stset tw, failure(dw)
mestreg x || g:, distribution(weibull)
capture noisily mestreg x || g:, distribution(weibull) nohr
capture noisily mestreg x || g:, distribution(weibull) time
capture noisily mestreg x || g:, distribution(weibull) tratio
capture noisily mestreg x || g:, distribution(exponential)
capture noisily mestreg x || g:, distribution(exponential) nohr
capture noisily mestreg x || g:, distribution(exponential) time
capture noisily mestreg x || g: x, distribution(weibull) covariance(unstructured)
capture noisily mestreg x || g:, distribution(weibull) intpoints(9)
capture noisily mestreg x || g:, distribution(weibull) level(90)
predict xbw, xb
capture noisily predict xbfe, fixedonly
capture noisily predict medw, median time
capture noisily predict hzw, hazard
capture noisily predict mnw, mean time
capture noisily predict ref, reffects
summarize xbw
ereturn list
* streg without random effects: the constant that Stata adds to the likelihood (VERIFY)
capture noisily streg x, distribution(weibull) nohr
capture noisily mestreg x, distribution(weibull)
stset tl, failure(dl)
capture noisily mestreg x || g:, distribution(lognormal)
capture noisily mestreg x || g:, distribution(loglogistic)
capture noisily mestreg x || g:, distribution(gamma)
capture noisily mestreg x || g:, distribution(lognormal) tratio
capture noisily mestreg x || g:, distribution(lognormal) nohr
capture noisily mestreg x || g:, distribution(lognormal) intpoints(9)
capture noisily streg x, distribution(lognormal)
capture noisily estat icc
ereturn list

* crossed random effects in mixed
clear
set seed 2907
set obs 200
gen a = ceil(runiform()*8)
gen b = ceil(runiform()*6)
gen x = invnormal(runiform())
bysort a: gen double ua = invnormal(runiform())*.8 if _n == 1
bysort a: replace ua = ua[1]
bysort b: gen double ub = invnormal(runiform())*.5 if _n == 1
bysort b: replace ub = ub[1]
gen y = 1 + .6*x + ua + ub + invnormal(runiform())
mixed y x || _all: R.a || _all: R.b
capture noisily mixed y x || _all: R.a || b:
capture noisily mixed y x || _all: R.a || _all: R.b, reml
capture noisily mixed y x || _all: R.a || _all: R.b, stddeviations
capture noisily mixed y x || a: || b:
capture noisily mixed y x || _all: R.a
capture noisily mixed y x || _all: R.a R.b
capture noisily mixed y x || a: || _all: R.b
capture noisily estat group
capture noisily estat icc
capture noisily estat recovariance
capture noisily predict yhat, fitted
capture noisily predict rr, reffects
capture noisily predict res, residuals
summarize yhat res
ereturn list

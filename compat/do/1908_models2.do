* Phase 22: models, part 2 (systems and constraints, choice models, exact methods, nonlinear models,
* frontier and treatment effects). Every command runs under capture noisily so that a difference in one
* command does not stop the rest of the case; the layouts marked VERIFY in the code are settled by the
* reference log. The data are built from rounded draws and from comparisons of runiform() draws (counts
* are sums of Bernoulli draws, because rpoisson() and rbinomial() do not reproduce Stata's draws), with
* real overdispersion and inefficiency so that no parameter sits at the boundary. Sums of rounded
* terms avoid halves (even tenths as coefficients) and exact zeros, so that the rounding does not depend
* on the last bit of a floating-point sum.
program define mksys
    clear
    set seed 1908
    set obs 400
    gen double x1 = round(invnormal(runiform()), .01)
    gen double x2 = round(invnormal(runiform()), .01)
    gen double x3 = round(.5*x1 + invnormal(runiform()), .01)
    gen double z1 = round(invnormal(runiform()), .01)
    gen double z2 = round(invnormal(runiform()), .01)
    gen double e1 = round(invnormal(runiform()), .01)
    gen double e2 = round(.5*e1 + invnormal(runiform()), .01)
    gen double y1 = round(1 + .8*x1 - .6*x2 + e1, .01)
    gen double y2 = round(-1 + .6*x1 + .4*x3 + e2, .01)
    gen double w = round(.5*z1 + .4*z2 + .3*x1 + .5*e1 + invnormal(runiform()), .01)
    gen double y3 = round(1 + .8*w + .4*x1 + e1, .01)
    gen double sdv = round(.5 + runiform(), .01)
    gen double y4 = round(1 + .8*x1 - .5*x2 + sdv*invnormal(runiform()), .01)
    gen double xt = round(x1 + .5*invnormal(runiform()), .01)
    gen double y5 = round(1 + .8*x1 + 1.2*e1, .01)
    replace y5 = y5 + 6 if _n <= 12
end

mksys
* ---- systems and constraints
capture noisily sureg (y1 x1 x2) (y2 x1 x3)
capture noisily sureg (y1 x1 x2) (y2 x1 x3), small dfk
capture noisily sureg (y1 x1 x2) (y2 x1 x3), isure
capture noisily reg3 (y1 x1 x2) (y2 x1 x3)
capture noisily reg3 (y1 x1 x2) (y2 x1 x3), 2sls
capture noisily reg3 (y3 w x1) (y1 x1 x2), inst(z1 z2 x1 x2)
capture noisily constraint define 1 x1 = x2
capture noisily constraint list
capture noisily cnsreg y1 x1 x2 x3, constraints(1)
capture noisily eivreg y4 x1 x2, r(x1 .9)
capture noisily vwls y4 x1 x2, sd(sdv)
capture noisily rreg y5 x1 x2
capture noisily rreg y5 x1 x2, genwt(rw)
* ---- choice models
clear
set seed 190801
set obs 300
gen int id = _n
gen double inc = round(invnormal(runiform()), .01)
expand 3
bysort id: gen byte alt = _n
gen double price = round(1 + runiform()*2, .01)
gen double qual = round(invnormal(runiform()), .01)
gen double u = -.8*price + .5*qual + .3*inc*(alt==2) - .4*inc*(alt==3) - ln(-ln(runiform()))
egen double umax = max(u), by(id)
gen byte choice = u == umax
egen byte rk = rank(u), by(id)
capture noisily clogit choice price qual, group(id)
capture noisily clogit choice price qual, group(id) or
capture noisily clogit choice price qual, group(id) vce(robust)
capture noisily rologit rk price qual, group(id)
capture noisily asclogit choice price qual, case(id) alternatives(alt) casevars(inc)
capture noisily asclogit choice price qual, case(id) alternatives(alt) casevars(inc) basealternative(2)
capture noisily asmprobit choice price qual, case(id) alternatives(alt) casevars(inc) intpoints(100)
clear
set seed 190802
set obs 500
gen double x1 = round(invnormal(runiform()), .01)
gen double x2 = round(invnormal(runiform()), .01)
gen double ua = .5*x1 - ln(-ln(runiform()))
gen double ub = .6*x2 - .4*x1 - ln(-ln(runiform()))
gen double uc = .3*x1 + .5*x2 - ln(-ln(runiform()))
gen byte yy = 1 + (ub > ua & ub > uc) + 2*(uc > ua & uc > ub)
capture noisily slogit yy x1 x2, dim(1)
capture noisily mprobit yy x1 x2
capture noisily mprobit yy x1 x2, baseoutcome(2)
* ---- exact methods
clear
set seed 190803
set obs 70
gen byte ge = runiform() < .5
gen byte cov = ceil(3*runiform())
gen byte ye = runiform() < invlogit(-.6 + .9*ge + .3*cov)
capture noisily exlogistic ye ge cov
capture noisily exlogistic ye ge cov, test(score)
gen double mc = exp(-.3 + .5*ge + .2*cov)
gen int ce = 0
forvalues j = 1/20 {
    replace ce = ce + (runiform() < mc/40)
}
capture noisily expoisson ce ge cov
* ---- nonlinear models
mksys
gen double yn = round(2 + 3*exp(-.8*x1) + .5*e1, .01)
capture noisily nl (yn = {b0} + {b1}*exp(-{b2}*x1)), initial(b0 1 b1 2 b2 .5)
capture noisily nlsur (y1 = {a0} + {a1}*x1 + {a2}*x2) (y2 = {c0} + {c1}*x1 + {c2}*x3)
capture noisily gmm (y1 - {b1}*x1 - {b2}*x2 - {b0}), instruments(x1 x2 z1 z2)
capture noisily gmm (y1 - {b1}*x1 - {b2}*x2 - {b0}), instruments(x1 x2 z1 z2) onestep
* ---- stochastic frontier
clear
set seed 190804
set obs 500
gen double x1 = round(invnormal(runiform()), .01)
gen double x2 = round(invnormal(runiform()), .01)
gen double v = round(.3*invnormal(runiform()), .001)
gen double ui = round(abs(.6*invnormal(runiform())), .001)
gen double yf = round(1 + .6*x1 + .3*x2 + v - ui, .001)
capture noisily frontier yf x1 x2
capture noisily frontier yf x1 x2, distribution(exponential)
* ---- treatment effects
clear
set seed 190805
set obs 900
gen double x = round(invnormal(runiform()), .01)
gen double z = round(invnormal(runiform()), .01)
gen double u = round(invnormal(runiform()), .01)
gen byte t = (.2 + .7*z + .5*x + u + .0001*runiform()) > 0
gen double yt = round(1 + x + 1.5*t + .5*t*x + .6*u + .7*invnormal(runiform()), .01)
capture noisily teffects ra (yt x) (t)
capture noisily teffects ra (yt x) (t), atet
capture noisily teffects ra (yt x) (t), pomeans
capture noisily teffects ipw (yt) (t x z)
capture noisily teffects ipw (yt) (t x z, probit)
capture noisily teffects ipwra (yt x) (t x z)
capture noisily teffects aipw (yt x) (t x z)
capture noisily teffects nnmatch (yt x) (t)
capture noisily teffects nnmatch (yt x) (t), atet nneighbors(2)
capture noisily teffects psmatch (yt) (t x z)
capture noisily eteffects (yt x) (t x z)
capture noisily etregress yt x, treat(t = z x)
capture noisily etregress yt x, treat(t = z x) twostep
gen double mu = exp(.3 + .4*x + .5*t + .5*(.6*u + .8*round(invnormal(runiform()), .01)))
gen int cy = 0
forvalues j = 1/30 {
    replace cy = cy + (runiform() < mu/30)
}
capture noisily etpoisson cy x, treat(t = z x)
gen double tt = round(exp(1 + .3*x - .5*t + .5*ln(-ln(runiform()))), .001) + .001
gen byte fl = tt < 6
replace tt = min(tt, 6)
stset tt, failure(fl)
capture noisily stteffects ra (x, weibull) (t)
capture noisily stteffects ipwra (x, exponential) (t x z)
capture noisily stteffects ipw (t x z), lognormal

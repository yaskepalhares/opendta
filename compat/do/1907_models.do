* Phase 21: models, part 1 (binary and fractional response, count data, censored, truncated and
* selection models, models with endogenous regressors). Every command runs under capture noisily so
* that a difference in one command does not stop the rest of the case; the layouts marked VERIFY in
* the code are settled by the reference log. The data are built from rounded draws, so that the last
* digit of the draws does not matter.
program define mkmod
    clear
    set seed 1907
    set obs 600
    gen int id = _n
    gen byte g = ceil(3*runiform())
    gen double x1 = round(invnormal(runiform()), .01)
    gen double x2 = round(invnormal(runiform()), .01)
    gen double x3 = round(.5*x1 + invnormal(runiform()), .01)
    gen double u = round(invnormal(runiform()), .01)
    gen double xs = round(runiform() - .5, .01)
    gen double ystar = round(.3 + .8*x1 - .5*x2 + u, .01)
    gen byte yb = ystar > 0
    gen byte yb2 = (.2 - .6*x1 + .5*x3 + round(invnormal(runiform()), .01)) > 0
    gen double mu = exp(.4 + .3*x1 - .2*x2)
    gen int cnt = rpoisson(mu)
    gen double mu2 = exp(.6 + .3*x1)
    gen int cnt2 = rpoisson(mu2)
    replace cnt2 = 0 if runiform() < .3 + .2*(x2 > 0)
    gen int ctr = cnt if cnt > 0
    gen double yfr = invlogit(.2 + .5*xs*2 - .3*x2 + .3*u)
    replace yfr = round(yfr, .001)
    gen byte yr = runiform() < invlogit(-1 + .4*xs)
    gen byte yrr = runiform() < .2 + .3*xs
    gen double yc = round(1 + .7*x1 + .5*x2 + 1.2*u, .01)
    gen double ylo = floor(yc)
    gen double yhi = ceil(yc)
    replace ylo = . if yc < -1
    replace yhi = . if yc > 3
    replace yhi = yhi + 1 if ylo == yhi & !missing(yhi)
    gen double yt = yc if yc > -1
    gen byte sel = (.3 + .6*x1 + .7*x3 + round(invnormal(runiform()), .01)) > 0
    gen double ysel = round(1 + .8*x2 + .6*u + .6*round(invnormal(runiform()), .01), .01) if sel
    gen byte ysb = (.2 + .5*x2 + .5*u + .6*round(invnormal(runiform()), .01)) > 0 if sel
    gen double z1 = round(invnormal(runiform()), .01)
    gen double z2 = round(invnormal(runiform()), .01)
    gen double v = round(.5*u + .8*invnormal(runiform()), .01)
    gen double w2 = round(.6*z1 + .4*z2 + .3*x1 + v, .01)
    gen double yiv = round(1 + .7*w2 + .4*x1 + u, .01)
    gen byte yivb = yiv > 1.5
    gen double yivt = max(yiv, 0)
    gen int yivc = rpoisson(exp(.2 + .3*w2 + .2*x1 + .3*u))
end

* ---- binary and fractional response
mkmod
capture noisily cloglog yb x1 x2
capture noisily cloglog yb x1 x2, eform
capture noisily cloglog yb x1 x2, vce(robust)
capture noisily predict double pcl, pr
capture noisily scobit yb x1 x2
capture noisily hetprobit yb x1 x2, het(x3)
capture noisily fracreg logit yfr xs x2
capture noisily fracreg probit yfr xs x2
capture noisily fracreg cloglog yfr xs x2
capture noisily binreg yr xs, or
capture noisily binreg yr xs, rr
capture noisily binreg yrr xs, rd
capture noisily binreg yrr xs, hr
mkmod
gen int pop = 20 + g
gen int pos = rbinomial(pop, invlogit(-.3 + .5*xs))
capture noisily glogit pos pop xs
* ---- count data
capture noisily zip cnt2 x1, inflate(x2)
capture noisily zip cnt2 x1, inflate(x2) irr
capture noisily zip cnt2 x1, inflate(x2) probit
capture noisily zip cnt2 x1, inflate(x2) vuong
capture noisily zinb cnt2 x1, inflate(x2)
capture noisily tpoisson ctr x1 x2, ll(0)
capture noisily tpoisson ctr x1 x2, ll(0) irr
capture noisily tnbreg ctr x1 x2, ll(0)
capture noisily gnbreg cnt x1 x2, lnalpha(x3)
* ---- interval, truncated and selection models
capture noisily intreg ylo yhi x1 x2
capture noisily intreg ylo yhi x1 x2, het(x3)
capture noisily truncreg yt x1 x2, ll(-1)
capture noisily heckman ysel x2, select(sel = x1 x3)
capture noisily heckman ysel x2, select(sel = x1 x3) twostep
capture noisily heckman ysel x2, select(sel = x1 x3) vce(robust)
capture noisily heckprobit ysb x2, select(sel = x1 x3)
capture noisily biprobit yb yb2 x1 x3
capture noisily biprobit (yb = x1 x2) (yb2 = x1 x3)
capture noisily predict double p11, p11
* ---- endogenous regressors
capture noisily ivregress liml yiv x1 (w2 = z1 z2)
capture noisily ivregress liml yiv x1 (w2 = z1 z2), vce(robust)
capture noisily ivregress gmm yiv x1 (w2 = z1 z2)
capture noisily ivregress gmm yiv x1 (w2 = z1 z2), wmatrix(unadjusted)
capture noisily ivprobit yivb x1 (w2 = z1 z2)
capture noisily ivtobit yivt x1 (w2 = z1 z2), ll(0)
capture noisily ivpoisson gmm yivc x1 (w2 = z1 z2)
capture noisily ivpoisson cfunction yivc x1 (w2 = z1 z2)

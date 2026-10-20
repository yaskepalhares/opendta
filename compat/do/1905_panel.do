* Phase 19: panel data (random-effects models by quadrature, linear panel estimators, frontier,
* unit-root tests, dynamic panel GMM and panel description). Every command runs under capture
* noisily so that a difference in one command does not stop the rest of the case; the layouts
* marked VERIFY in the code are settled by the reference log.
* The data are built from rounded draws, so that the last digit of the draws does not matter.
program define mkpanel
    clear
    set seed 1905
    set obs 480
    gen int id = ceil(_n/8)
    gen int t = _n - (id-1)*8 + 1999
    xtset id t
    gen double u = round(invnormal(runiform()), .01)
    bysort id (t): replace u = u[1]
    gen double e = round(invnormal(runiform()), .01)
    gen double x1 = round(invnormal(runiform()), .01)
    gen double z = round(invnormal(runiform()), .01)
    gen double x2 = round(.5*z + invnormal(runiform()) + .5*u, .01)
    gen double w = round(id/20, .01)
    gen double y = round(1 + .5*x1 + .3*x2 + .2*w + u + e, .01)
    gen byte yb = runiform() < invlogit(-.2 + .8*x1 + .7*u)
    gen double lat = round(.5*x1 + .6*u + e, .01)
    gen byte yo = (lat > -.4) + (lat > .6)
    gen double yt = max(0, lat)
    gen double ylo = floor(lat)
    gen double yhi = ylo + 1
    gen double yc = floor(-ln(runiform()) * exp(.3*x1 + .5*u) * 3)
    gen double ly = round(1 + .5*x1 + .3*x2 + .3*e - .3*abs(invnormal(runiform())), .01)
end

program define mkroot
    clear
    set seed 19051
    set obs 360
    gen int id = ceil(_n/30)
    gen int t = _n - (id-1)*30
    xtset id t
    gen double e = round(invnormal(runiform()), .01)
    gen double rw = e
    bysort id (t): replace rw = round(rw[_n-1] + e, .01) if _n > 1
    gen double st = e
    bysort id (t): replace st = round(.5*st[_n-1] + e, .01) if _n > 1
end

program define mkdyn
    clear
    set seed 19052
    set obs 1200
    gen int id = ceil(_n/6)
    gen int t = _n - (id-1)*6
    xtset id t
    gen double u = round(invnormal(runiform()), .01)
    bysort id (t): replace u = u[1]
    gen double x = round(invnormal(runiform()), .01)
    gen double w = round(invnormal(runiform()), .01)
    gen double y = round(u + x + invnormal(runiform()), .01)
    bysort id (t): replace y = round(.5*y[_n-1] + x + .3*w + u + invnormal(runiform()), .01) if _n > 1
end

* ---- random effects by quadrature
mkpanel
capture noisily xtprobit yb x1, re
capture noisily xtprobit yb x1, re intmethod(ghermite) intpoints(8)
capture noisily quadchk
capture noisily xtprobit yb x1, re vce(robust)
capture noisily xtcloglog yb x1, re
capture noisily xtologit yo x1
capture noisily xtoprobit yo x1
capture noisily xttobit yt x1, ll(0)
capture noisily xtintreg ylo yhi x1
capture noisily xtnbreg yc x1, re
capture noisily xtnbreg yc x1, fe
capture noisily xtnbreg yc x1, re irr

* ---- linear panel estimators
mkpanel
capture noisily xtgls y x1 x2, panels(iid)
capture noisily xtgls y x1 x2, panels(hetero)
capture noisily xtgls y x1 x2, panels(correlated)
capture noisily xtgls y x1 x2, panels(hetero) corr(ar1)
capture noisily xtgls y x1 x2, panels(iid) corr(psar1) igls
capture noisily xtpcse y x1 x2
capture noisily xtpcse y x1 x2, corr(ar1)
capture noisily xtpcse y x1 x2, hetonly
capture noisily xtregar y x1 x2, fe
capture noisily xtregar y x1 x2, re
capture noisily xtrc y x1 x2
capture noisily xtivreg y x1 (x2 = z), fe
capture noisily xtivreg y x1 (x2 = z)
capture noisily xtivreg y x1 (x2 = z), ec2sls
capture noisily xtivreg y x1 (x2 = z), be
capture noisily xtivreg y x1 (x2 = z), fd
capture noisily xthtaylor y x1 x2 w, endog(x2)
capture noisily xthtaylor y x1 x2 w, endog(x2) amacurdy

* ---- stochastic frontier
mkpanel
capture noisily xtfrontier ly x1 x2
capture noisily xtfrontier ly x1 x2, tvd
capture noisily xtfrontier ly x1 x2, cost

* ---- unit-root tests
mkroot
capture noisily xtunitroot ht rw
capture noisily xtunitroot ht st
capture noisily xtunitroot breitung rw
capture noisily xtunitroot breitung st, trend
capture noisily xtunitroot hadri rw
capture noisily xtunitroot hadri st, trend
capture noisily xtunitroot fisher rw, lags(1)
capture noisily xtunitroot fisher st, lags(1) trend
capture noisily xtunitroot llc rw
capture noisily xtunitroot ips rw

* ---- dynamic panel-data estimation
mkdyn
capture noisily xtabond y x
capture noisily estat abond
capture noisily estat sargan
capture noisily xtabond y x, lags(2) pre(w) maxldep(3)
capture noisily xtabond y x, twostep vce(robust)
capture noisily estat abond
capture noisily xtdpdsys y x
capture noisily estat abond
capture noisily xtdpdsys y x, twostep vce(robust)
capture noisily xtdpd y L.y x, dgmmiv(y) div(x)

* ---- panel description
mkpanel
capture noisily xttab yo
capture noisily xtdata y x1, fe clear
summarize y x1
mkpanel
capture noisily xtdata y x1, be clear
summarize y x1

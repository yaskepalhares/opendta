* Phase 18: time series (data tools, filters and smoothers, spectral analysis, unit-root tests,
* prais, VAR postestimation, vec, svar, irf extensions, arch, arfima, arima extensions,
* state-space models and mgarch). Every command runs under capture noisily so that a difference
* in one command does not stop the rest of the case; the layouts marked VERIFY in the code are
* settled by the reference log.
* The data are built from rounded draws, so that the last digit of the draws does not matter.
program define mkts
    clear
    set seed 1904
    set obs 160
    gen int t = _n + 599
    format t %tm
    tsset t
    gen double e1 = round(invnormal(runiform()), .01)
    gen double e2 = round(invnormal(runiform()), .01)
    gen double e3 = round(invnormal(runiform()), .01)
    gen double x = round(invnormal(runiform()), .01)
    gen double w = round(sum(e3), .01)
    gen double c1 = round(w + .5*e1, .01)
    gen double c2 = round(.5*w + .5*e2, .01)
    gen double y = 0
    gen double z = 0
    forvalues i = 2/160 {
        quietly replace y = round(.5*y[`i'-1] + .2*z[`i'-1] + e1, .01) in `i'
        quietly replace z = round(.1*y[`i'-1] + .4*z[`i'-1] + e2, .01) in `i'
    }
    gen double u = 0
    forvalues i = 2/160 {
        quietly replace u = round(.6*u[`i'-1] + e3 + .3*e3[`i'-1], .01) in `i'
    }
    gen double r = round(2 + .5*x + u, .01)
    * GARCH(1,1) series
    gen double h = 1
    gen double g = 0
    forvalues i = 2/160 {
        quietly replace h = .2 + .2*g[`i'-1]^2 + .5*h[`i'-1] in `i'
        quietly replace g = round(sqrt(h)*e1, .01) in `i'
    }
    * GARCH(1,1) series with a shock correlated with e1
    gen double h2 = 1
    gen double g2 = 0
    forvalues i = 2/160 {
        quietly replace h2 = .1 + .2*g2[`i'-1]^2 + .6*h2[`i'-1] in `i'
        quietly replace g2 = round(sqrt(h2)*(.5*e1 + .87*e2), .01) in `i'
    }
    * seasonal series (quarterly pattern)
    gen double sa = round(10 + 2*mod(_n, 4) + .1*_n + e2, .01)
    * the series below are for the commands that need data with the feature they estimate
    * fat-tailed GARCH(1,1) series (arch with t errors, egarch)
    gen double eh = round(e1*(1 + 2*(e3 > 1.2)), .01)
    gen double hh0 = 1
    gen double gt = 0
    forvalues i = 2/160 {
        quietly replace hh0 = .2 + .2*gt[`i'-1]^2 + .5*hh0[`i'-1] in `i'
        quietly replace gt = round(sqrt(hh0)*eh, .01) in `i'
    }
    * cycle (AR(2) with complex roots)
    gen double cy = 0
    forvalues i = 3/160 {
        quietly replace cy = round(1.2*cy[`i'-1] - .6*cy[`i'-2] + e2, .01) in `i'
    }
    * local linear trend plus noise
    gen double lv = 0
    gen double sl = 0
    forvalues i = 2/160 {
        quietly replace sl = sl[`i'-1] + .1*e2 in `i'
        quietly replace lv = lv[`i'-1] + sl[`i'-1] + .3*e1 in `i'
    }
    gen double tr = round(lv + e3, .01)
    * stochastic seasonal series
    gen double ss = 0
    forvalues i = 5/160 {
        quietly replace ss = round(ss[`i'-4] + e1 + .5*e1[`i'-4], .01) in `i'
    }
    * two series whose correlation moves (dynamic conditional correlation)
    gen double q11 = 1
    gen double q12 = .4
    gen double q22 = 1
    gen double z1 = 0
    gen double z2 = 0
    gen double hd1 = 1
    gen double hd2 = 1
    gen double d1 = 0
    gen double d2 = 0
    forvalues i = 2/160 {
        local j = `i' - 1
        quietly replace q11 = .05 + .15*z1[`j']^2 + .8*q11[`j'] in `i'
        quietly replace q12 = .02 + .15*z1[`j']*z2[`j'] + .8*q12[`j'] in `i'
        quietly replace q22 = .05 + .15*z2[`j']^2 + .8*q22[`j'] in `i'
        quietly replace z1 = e1 in `i'
        quietly replace z2 = q12/sqrt(q11*q22)*e1 + sqrt(1 - q12^2/(q11*q22))*e2 in `i'
        quietly replace hd1 = .1 + .2*d1[`j']^2 + .6*hd1[`j'] in `i'
        quietly replace hd2 = .15 + .25*d2[`j']^2 + .5*hd2[`j'] in `i'
        quietly replace d1 = round(sqrt(hd1)*z1, .01) in `i'
        quietly replace d2 = round(sqrt(hd2)*z2, .01) in `i'
    }
end
set more off
set linesize 255
mkts

* --- data tools ---
preserve
drop if mod(t, 7) == 0 in 20/60
capture noisily tsreport
capture noisily tsreport y z, detail
capture noisily tsfill
capture noisily tsreport
restore
preserve
capture noisily tsappend, add(6)
capture noisily list t y in -8/l
restore
preserve
capture noisily rolling _b _se, window(80) step(40) clear: regress y x
capture noisily list
restore
preserve
capture noisily rolling mean=r(mean), window(60) recursive clear: summarize y
capture noisily list in 1/5
restore

* --- filters and smoothers ---
capture noisily tsfilter hp yhp = y, smooth(100) trend(ytr)
capture noisily summarize yhp ytr
capture noisily tsfilter bk ybk = y, minperiod(6) maxperiod(32) smaorder(4)
capture noisily summarize ybk
capture noisily tsfilter cf ycf = y, minperiod(6) maxperiod(32)
capture noisily summarize ycf
capture noisily tsfilter bw ybw = y, maxperiod(32) order(2)
capture noisily summarize ybw
capture noisily tssmooth ma ma1 = y, window(2 1 2)
capture noisily tssmooth exponential ex1 = y, parms(.4)
capture noisily tssmooth exponential ex2 = y
capture noisily tssmooth dexponential dx1 = y, parms(.3)
capture noisily tssmooth hwinters hw1 = sa, parms(.4 .1)
capture noisily tssmooth shwinters sw1 = sa, period(4) parms(.4 .1 .1)
capture noisily tssmooth nl nl1 = y, smoother(4253h)
capture noisily summarize ma1 ex1 ex2 dx1 hw1 sw1 nl1

* --- spectral analysis and white-noise tests ---
capture noisily pergram y, gen(pg) nograph
capture noisily list pg in 1/5
capture noisily cumsp y, gen(cs)
capture noisily xcorr y z, lags(4) table noplot
capture noisily wntestb y, table
capture noisily wntestq y
capture noisily wntestq y, lags(5)
capture noisily arima r x, ar(1) ma(1) nolog
capture noisily psdensity pd omega
capture noisily summarize pd omega

* --- unit roots and prais ---
capture noisily dfgls w, maxlag(4)
capture noisily dfgls w, maxlag(4) notrend
capture noisily pperron w
capture noisily pperron w, trend regress lags(3)
capture noisily pperron y, noconstant
capture noisily prais r x
capture noisily prais r x, corc
capture noisily prais r x, rhotype(dw) twostep
capture noisily prais r x, ssesearch
capture noisily prais r x, rhotype(theil) nodw

* --- VAR postestimation ---
capture noisily var y z, lags(1/2)
capture noisily varlmar
capture noisily varlmar, mlag(3)
capture noisily varnorm
capture noisily varnorm, jbera skewness kurtosis
capture noisily varstable, amat(compmat)
capture noisily matrix list compmat
capture noisily varwle
capture noisily fcast compute f_, step(6)
capture noisily list t y f_y f_y_SE in -8/l
capture noisily varbasic y z, lags(1/2) nograph

* --- vecrank, vec ---
capture noisily vecrank c1 c2, lags(2)
capture noisily vecrank c1 c2, lags(2) trend(rconstant) max
capture noisily vecrank c1 c2, lags(2) ic
capture noisily vec c1 c2, rank(1) lags(2)
capture noisily vec c1 c2, rank(1) lags(2) trend(rconstant) noidtest
capture noisily vec c1 c2, rank(1) lags(2) alpha
capture noisily fcast compute v_, step(4)
capture noisily list t c1 v_c1 in -5/l
capture noisily irf create vecirf, set(odta_vec) replace step(6)
capture noisily irf table irf, irf(vecirf) impulse(c1) response(c2)

* --- svar and irf extensions ---
matrix A1 = (1, 0 \ ., 1)
matrix B1 = (., 0 \ 0, .)
capture noisily svar y z, lags(1/2) aeq(A1) beq(B1)
capture noisily irf create sv, set(odta_sv) replace step(6)
capture noisily irf table oirf, irf(sv) impulse(y) response(z)
capture noisily var y z, lags(1/2)
capture noisily irf create vv, set(odta_sv) replace step(6) bs reps(50) nodots
capture noisily irf table oirf, irf(vv) impulse(y) response(z)
capture noisily irf create vv2, set(odta_sv) replace step(6)
capture noisily irf ctable (vv2 y z irf) (vv2 y z oirf), noci
capture noisily irf describe
capture noisily irf drop vv2
capture erase odta_vec.irf
capture erase odta_sv.irf

* --- arima extensions ---
capture noisily arima ss, arima(1,0,0) sarima(0,1,1,4) nolog
capture noisily arima sa, ar(1) mar(1, 4) nolog
capture noisily arima r x, ar(1) ma(1) condition nolog
capture noisily arima r x, ar(1) diffuse nolog

* --- arch family ---
* (the periods added by fcast and tsappend are dropped: the sample of these commands is the data)
keep in 1/160
capture noisily arch g, arch(1) garch(1) nolog
capture noisily arch gt, arch(1) garch(1) distribution(t) nolog
capture noisily arch g, arch(1/2) nolog
capture noisily arch gt, earch(1) egarch(1) nolog
capture noisily arch g, arch(1) garch(1) tarch(1) nolog
capture noisily arch g, arch(1) garch(1) archm nolog
capture noisily arch r x, arch(1) garch(1) ar(1) nolog
capture noisily predict hh, variance

* --- arfima ---
capture noisily arfima g, nolog
capture noisily arfima r x, ar(1) nolog

* --- state-space models ---
capture noisily ucm y, model(llevel) nolog
capture noisily ucm y, model(rwalk) nolog
capture noisily ucm tr, model(lltrend) nolog
capture noisily ucm ss, model(llevel) seasonal(4) nolog
capture noisily ucm y x, model(ntrend) nolog
capture noisily ucm cy, model(ntrend) cycle(1) nolog
capture noisily sspace (z1 L.z1, state noconstant) (y z1, noconstant), covstate(identity) nolog
capture noisily dfactor (y z = , noconstant) (f = , ar(1)), nolog

* --- mgarch ---
capture noisily mgarch ccc (g g2 = , noconstant), arch(1) garch(1) nolog
capture noisily mgarch dcc (d1 d2 = , noconstant), arch(1) garch(1) nolog
capture noisily mgarch dvech (g g2 = , noconstant), arch(1) garch(1) nolog

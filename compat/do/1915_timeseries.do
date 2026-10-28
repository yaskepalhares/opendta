* Time series, second part: prais with vce(robust), vce(hc2), vce(hc3) and vce(cluster); arch with
* constraints() and from(); arfima with vce(robust); svar with aconstraints() and
* bconstraints(); sspace and dfactor with constraints() and predict; ucm with the frequency()
* suboption of cycle(), predict and psdensity. Every command whose behavior the manual does not
* settle runs under capture noisily, so that one Stata run shows the whole case.
clear
set seed 1915
set obs 300
gen t = _n
tsset t
gen x = round(invnormal(runiform()),.01)
gen e = round(invnormal(runiform()),.001)
gen u = e
replace u = .6*u[_n-1] + e in 2/l
gen y = 1 + 2*x + u
gen g = ceil(_n/25)

* prais
capture noisily prais y x, vce(robust)
capture noisily ereturn list
capture noisily prais y x, vce(cluster g)
capture noisily prais y x, vce(hc2) corc
capture noisily prais y x, vce(hc3) twostep

* arch: constraints() and from()
gen h = 1
gen ya = 0
replace h = .2 + .3*(e[_n-1]*sqrt(h[_n-1]))^2 + .4*h[_n-1] in 2/l
replace ya = 1 + .5*x + e*sqrt(h)
capture noisily arch ya x, arch(1) garch(1)
matrix b0 = e(b)
constraint define 1 [ARCH]L.arch = .3
capture noisily arch ya x, arch(1) garch(1) constraints(1)
capture noisily matrix list e(Cns)
capture noisily arch ya x, arch(1) garch(1) from(b0)
capture noisily arch ya x, arch(1) garch(1) from(b0, copy)

* arfima
capture noisily arfima u, ar(1) vce(robust)
capture noisily ereturn list

* svar with constraints
gen e2 = round(invnormal(runiform()),.01)
gen y2 = e2
replace y2 = .3*u[_n-1] + .4*y2[_n-1] + .5*e + e2 in 2/l
matrix A = (1,0\.,1)
matrix B = (.,0\0,.)
capture noisily svar u y2, lags(1) aeq(A) beq(B)
constraint define 11 [a_1_1]_cons = 1
constraint define 12 [a_1_2]_cons = 0
constraint define 13 [a_2_2]_cons = 1
constraint define 14 [b_1_2]_cons = 0
constraint define 15 [b_2_1]_cons = 0
capture noisily svar u y2, lags(1) aconstraints(11/13) bconstraints(14 15)
capture noisily ereturn list
capture noisily svar u y2, lags(1) aeq(A) bconstraints(14 15)

* sspace and dfactor: constraints() and predict
constraint define 21 [u]s = 1
capture noisily sspace (s L.s, state noconstant) (u s, noerror noconstant), constraints(21)
capture noisily predict uh
capture noisily predict ur, residuals
capture noisily predict urs, rstandard
capture noisily summarize uh ur urs
gen u2 = round(u + invnormal(runiform()), .01)
capture noisily dfactor (u u2 = , noconstant) (f = , ar(1))
capture noisily predict f1
capture noisily predict f2, equation(u2)
capture noisily predict rs2, rstandard equation(u2)
capture noisily summarize f1 f2 rs2

* ucm: frequency(), predict, psdensity
gen yc = round(10 + 2*cos(_pi*_n/8) + invnormal(runiform()), .01)
capture noisily ucm yc, model(rwalk) cycle(1, frequency(.4))
capture noisily predict yh
capture noisily predict yr, residuals
capture noisily predict yrs, rstandard
capture noisily summarize yh yr yrs
capture noisily psdensity sd w
capture noisily summarize sd w
capture noisily psdensity sp w2, pspectrum
capture noisily summarize sp w2
capture noisily psdensity sd3 w3, range(.1 1)
capture noisily ucm yc, model(rwalk) cycle(1, frequency(4))

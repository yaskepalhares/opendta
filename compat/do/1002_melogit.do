* Phase 10: melogit, meprobit and mepoisson (adaptive Gauss-Hermite quadrature)
clear
set seed 1003
set obs 40
gen g = _n
gen u = invnormal(runiform())*.9
gen v = invnormal(runiform())*.3
expand 15
gen x = invnormal(runiform())
gen yb = runiform() < invlogit(-.4 + .8*x + u + v*x)
* Poisson counts by inversion of the distribution function (rpoisson() draws differ from Stata's)
gen double lam = exp(.2 + .5*x + u/2)
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
melogit yb x || g:
melogit yb x || g:, or
estat icc
meprobit yb x || g:
mepoisson yc x || g:
melogit yb x || g: x, covariance(unstructured)
predict p
capture noisily predict re, reffects
predict rex rec, reffects
summarize p rex rec
ereturn list

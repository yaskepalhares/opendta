* Phase 10: melogit, meprobit and mepoisson (adaptive Gauss-Hermite quadrature)
clear
set seed 1003
set obs 40
gen g = _n
gen u = rnormal()*.9
gen v = rnormal()*.3
expand 15
gen x = rnormal()
gen yb = runiform() < invlogit(-.4 + .8*x + u + v*x)
gen yc = rpoisson(exp(.2 + .5*x + u/2))
melogit yb x || g:
melogit yb x || g:, or
estat icc
meprobit yb x || g:
mepoisson yc x || g:
melogit yb x || g: x, covariance(unstructured)
predict p
predict re, reffects
summarize p re
ereturn list

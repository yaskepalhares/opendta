* Phase 10: mixed (ML, REML, random slopes, covariances, nested levels),
* estat icc/group, predict
clear
set seed 1001
set obs 30
gen id = _n
gen u0 = rnormal()*2
gen u1 = rnormal()*.5
expand 8
bysort id: gen week = _n
gen x = runiform()
gen y = 10 + u0 + (1 + u1)*week + 2*x + rnormal()*1.5
mixed y week x || id:
mixed y week x || id:, reml
estat icc
mixed y week x || id: week
mixed y week x || id: week, covariance(unstructured)
mixed y week x || id: week, covariance(unstructured) stddeviations
ereturn list
predict yf, fitted
predict re, reffects
summarize yf re

clear
set seed 1002
set obs 20
gen school = _n
gen us = rnormal()*1.5
expand 4
bysort school: gen class = _n
gen uc = rnormal()*.8
expand 6
gen x = rnormal()
gen y = 5 + us + uc + 0.7*x + rnormal()
mixed y x || school: || class:
estat icc
estat group
mixed y x || school:
estat icc
predict re, reffects
predict yf, fitted
predict r, residuals
summarize re yf r

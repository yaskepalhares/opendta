* Phase 10: arima with MA terms (Stata's starting values for MA terms are undocumented)
clear
set seed 1101
set obs 150
gen t = _n
tsset t
gen e = invnormal(runiform())
gen u = e
replace u = 0.6*u[_n-1] + e + 0.3*e[_n-1] in 2/l
gen x = invnormal(runiform())
gen y = 2 + 0.5*x + u
gen w = sum(invnormal(runiform()) + .1)
gen z = 0.4*L.y + invnormal(runiform()) in 2/l
arima y x, ar(1) ma(1)
arima w, arima(1,1,1)

* Phase 10: time series: arima, var, vargranger, varsoc, dfuller, corrgram, newey
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
arima y, ar(1 2) vce(oim)
var y z, lags(1/2)
vargranger
varsoc y z
dfuller w
dfuller w, trend lags(2) regress
dfuller y, noconstant
dfuller w, drift
corrgram y, lags(8)
newey y x, lag(3)
ereturn list
predict yhat
predict res, residuals
summarize yhat res

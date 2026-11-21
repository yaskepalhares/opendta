* lincom with c.x and eform: the name c.x and what r() holds after eform (Stata 14).
clear all
set more off
set seed 20
set obs 400
gen x = round(invnormal(runiform()), .01)
gen t = mod(_n, 10)
gen y = runiform() < .3 + .02 * t
capture noisily poisson y c.t x
capture noisily lincom c.t, eform
capture noisily return list
capture noisily lincom t, eform
capture noisily lincom c.t + x

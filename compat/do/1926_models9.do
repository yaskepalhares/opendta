* Models, part 9 (block 20g): ivprobit and ivtobit with twostep (Newey's minimum chi-squared
* estimator), also with more than one endogenous regressor. No reference until the Stata run.
* Every command whose behavior the manual does not settle runs under capture noisily, so that
* one Stata run shows the whole case.
clear
set obs 3000
set seed 1926
gen x = round(invnormal(runiform()), .01)
gen z1 = round(invnormal(runiform()), .01)
gen z2 = round(invnormal(runiform()), .01)
gen v = round(invnormal(runiform()), .01)
gen v2 = round(invnormal(runiform()), .01)
gen u = round(.5*v + sqrt(.75)*invnormal(runiform()), .01)
gen y2 = round(.5*x + z1 + .5*z2 + v, .01)
gen y3 = round(-.3*x + .4*z1 + z2 + v2, .01)
gen ys = .5 + .4*y2 - .6*x + u
gen yb = ys > 0
gen yt = max(ys, 0)
gen ys3 = .5 + .4*y2 - .3*y3 - .6*x + u
gen yb3 = ys3 > 0
capture noisily ivprobit yb x (y2 = z1 z2), twostep
capture noisily ereturn list
capture noisily ivprobit yb x (y2 = z1 z2), twostep first
capture noisily ivtobit yt x (y2 = z1 z2), ll(0) twostep
capture noisily ereturn list
capture noisily ivprobit yb3 x (y2 y3 = z1 z2), twostep
capture noisily ivtobit yt x (y2 y3 = z1 z2), ll(0) twostep
capture noisily ivprobit yb3 x (y2 y3 = z1 z2)

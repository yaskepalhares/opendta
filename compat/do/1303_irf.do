* irf after var: create, table, describe, rename and drop (deterministic data)
clear
set seed 1303
set obs 200
gen t = _n
tsset t
gen e1 = invnormal(runiform())
gen e2 = invnormal(runiform())
gen y1 = 0
gen y2 = 0
replace y1 = .5*L.y1 + .2*L.y2 + e1 in 2/l
replace y2 = .3*L.y1 + .4*L.y2 + .5*e1 + e2 in 2/l
summarize y1 y2
var y1 y2, lags(1/2)
irf create base, set(odta_irf, replace) step(6)
irf table oirf, impulse(y1) response(y2)
irf table irf, impulse(y2) response(y1)
irf table fevd, impulse(y1) response(y2)
irf table cirf coirf, impulse(y1) response(y1) noci stderror
irf create rev, order(y2 y1) step(6)
irf table oirf, impulse(y1) response(y2)
irf describe
irf rename rev other
irf drop other
irf describe
erase odta_irf.irf

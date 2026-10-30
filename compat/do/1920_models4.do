* Models, part 4 (block 20b): gmm with igmm and derivatives(), ivpoisson cfunction with more
* than one endogenous regressor, ivregress gmm with wmatrix(hac kernel opt). No reference until
* the Stata run. Every command whose behavior the manual does not settle runs under capture
* noisily, so that one Stata run shows the whole case.
clear
set seed 1920
set obs 400
gen t = _n
tsset t
gen x = round(invnormal(runiform()), .01)
gen z = round(invnormal(runiform()), .01)
gen w = round(invnormal(runiform()), .01)
gen u = round(invnormal(runiform()), .01)
gen e = round(invnormal(runiform()), .01)
gen e2 = round(invnormal(runiform()), .01)
gen e3 = round(invnormal(runiform()), .01)
gen yy = 1 + x + e
gen xe = round(x + .5*e + z + .3*w, .01)
gen w2 = round(x + z + .5*e2 + .3*w, .01)
gen w3 = round(x - z + .5*e3 + .3*u, .01)
gen double lam = exp(.2 + .3*x + .2*w2 - .1*w3)
gen double pu = runiform()
gen double pk = exp(-lam)
gen double cp = pk
gen yc = 0
forvalues k = 1/40 {
    quietly replace yc = yc + 1 if pu > cp
    quietly replace pk = pk*lam/`k'
    quietly replace cp = cp + pk
}
drop lam pu pk cp

* gmm: iterated weight matrix and analytic derivatives
gmm (yy - {b1}*xe - {b0}), instruments(z w)
capture noisily ereturn list
capture noisily gmm (yy - {b1}*xe - {b0}), instruments(z w) igmm
capture noisily ereturn list
capture noisily gmm (yy - {b1}*xe - {b0}), instruments(z w) igmm igmmeps(1e-4) igmmweps(1e-4) igmmiterate(10)
capture noisily gmm (yy - {b1}*xe - {b0}), instruments(z w) derivatives(/b1 = -xe) derivatives(/b0 = -1)
capture noisily gmm (yy - {b1}*xe - {b0}), instruments(z w) derivatives(/b1 = -xe)
capture noisily gmm (yy - {b1}*xe - {b0}), instruments(z w) igmm derivatives(/b1 = -xe) derivatives(/b0 = -1)

* ivpoisson cfunction with two endogenous regressors
capture noisily ivpoisson cfunction yc x (w2 w3 = z w u)
capture noisily ereturn list
capture noisily ivpoisson cfunction yc x (w2 w3 = z w u), vce(cluster t)

* ivregress gmm with the Newey-West (1994) lag choice
capture noisily ivregress gmm yy (xe = z w), wmatrix(hac bartlett opt)
capture noisily ereturn list
capture noisily ivregress gmm yy (xe = z w), wmatrix(hac parzen opt)
capture noisily ivregress gmm yy (xe = z w), wmatrix(hac quadraticspectral opt)

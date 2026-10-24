* Models, part 3 (block 20a): tpoisson and tnbreg with ul() and variable limits, tnbreg and
* nbreg with dispersion(constant), glogit with or, offset() in intreg, truncreg, heckman and
* biprobit, heckman twostep with mills(), sureg and reg3 with constraints(), ivregress gmm with
* wmatrix(hac), center and igmm, frontier with uhet(), vhet() and cm(), glm and binreg with mu().
* Counts are built by inversion (the rpoisson() draws differ from Stata's). Every command whose
* behavior the manual does not settle runs under capture noisily, so that one Stata run shows the
* whole case.
clear
set seed 1918
set obs 400
gen t = _n
tsset t
gen x = round(invnormal(runiform()), .01)
gen z = round(invnormal(runiform()), .01)
gen o = round(.3*runiform(), .01)
gen w = round(runiform(), .01)
gen double lam = exp(.8 + .4*x)
gen double pu = runiform()
gen double pk = exp(-lam)
gen double cp = pk
gen y = 0
forvalues k = 1/40 {
    quietly replace y = y + 1 if pu > cp
    quietly replace pk = pk*lam/`k'
    quietly replace cp = cp + pk
}
drop lam pu pk cp
gen y1 = 1 + .5*x + .5*z + invnormal(runiform())
gen y2 = 2 + .5*x - z + invnormal(runiform())

* truncated counts
tpoisson y x if y > 0
capture noisily tpoisson y x if y > 0 & y < 6, ul(6)
capture noisily predict n1
capture noisily summarize n1
gen l = cond(x > 0, 1, 0)
capture noisily tpoisson y x if y > l, ll(l)
capture noisily tnbreg y x if y > 0, dispersion(constant)
capture noisily tnbreg y x if y > 0 & y < 8, ul(8)
capture noisily ereturn list
nbreg y x, dispersion(constant)
ereturn list
capture noisily nbreg y x, dispersion(constant) irr

* glogit, or
gen pop = y + 3
glogit y pop x in 1/80
glogit y pop x in 1/80, or

* offsets
gen lo = floor(y1)
gen hi = lo + 1
intreg lo hi x z, offset(o)
capture noisily truncreg y1 x z, ll(0) offset(o)
gen s = (x + invnormal(runiform())) > 0
gen yh = y1 if s
heckman yh x, select(s = x z) twostep mills(imr)
summarize imr
capture noisily heckman yh x, select(s = x z) mills(imr2)
capture noisily heckman yh x, select(s = x z) offset(o)
gen b1 = y1 > 1
gen b2 = y2 > 2
capture noisily biprobit b1 b2 x z, offset1(o)
capture noisily biprobit b1 b2 x z, offset2(o)

* constraints in sureg and reg3
constraint define 1 [y1]x = [y2]x
sureg (y1 x z) (y2 x z), constraints(1)
ereturn list
capture noisily reg3 (y1 x z) (y2 x z), constraints(1)
capture noisily reg3 (y1 x z) (y2 x z), constraints(1) 2sls

* ivregress gmm
gen e = invnormal(runiform())
gen yy = 1 + x + e + .2*w
gen xe = x + .5*e + z + .3*w
ivregress gmm yy (xe = z w), wmatrix(hac bartlett 2)
ereturn list
capture noisily ivregress gmm yy (xe = z w), wmatrix(hac parzen 3)
capture noisily ivregress gmm yy (xe = z w), wmatrix(hac quadraticspectral 2)
capture noisily ivregress gmm yy (xe = z w), center
capture noisily ivregress gmm yy (xe = z w), igmm
capture noisily ivregress gmm yy (xe = z w), wmatrix(unadjusted) center

* frontier with heteroskedasticity
gen u = abs(invnormal(runiform()))*exp(.5*w)
gen yf = 1 + .5*x + .3*invnormal(runiform()) - .5*u
frontier yf x, uhet(w)
ereturn list
capture noisily frontier yf x, vhet(w)
capture noisily frontier yf x, uhet(w) vhet(w)
capture noisily frontier yf x, distribution(tnormal) cm(w)
capture noisily frontier yf x, cm(w)

* glm and binreg mu()
gen yb = y1 > 1
gen m0 = .5
capture noisily glm yb x, family(binomial) mu(m0)
capture noisily binreg yb x, mu(m0)

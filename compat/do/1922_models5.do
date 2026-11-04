* Models, part 5 (block 20c): biprobit with partial and constraints(), frontier with
* constraints(). No reference until the Stata run. Every command whose behavior the manual does
* not settle runs under capture noisily, so that one Stata run shows the whole case.
clear
set seed 1922
set obs 2000
gen x1 = round(invnormal(runiform()), .01)
gen x2 = round(invnormal(runiform()), .01)
gen x3 = round(invnormal(runiform()), .01)
gen u1 = round(invnormal(runiform()), .01)
gen u2 = round(.3*u1 + sqrt(.91)*invnormal(runiform()), .01)
gen y1 = (.8 + x2 + .3*x1 + u1 > 0)
gen y2 = (.5 - x3 + .3*x1 + u2 > 0)
gen e = round(invnormal(runiform())*.3, .01) - abs(round(invnormal(runiform())*.4, .01))
gen yf = 1 + .5*x1 + .3*x2 + e

* biprobit with constraints(): the constraint list above the table, e(Cns), the Wald test
constraint 1 [y1]x1 = [y2]x1
capture noisily biprobit (y1 = x1 x2) (y2 = x1 x3), constraints(1)
capture noisily ereturn list
capture noisily matrix list e(Cns)
constraint 2 [athrho]_cons = .25
capture noisily biprobit y1 y2 x1 x2 x3, constraints(2)

* biprobit with partial observability: only y1*y2 is used
capture noisily biprobit (y1 = x1 x2) (y2 = x1 x3), partial
capture noisily ereturn list
capture noisily predict double p11, p11
capture noisily summarize p11

* frontier with constraints() on the frontier and on a variance
constraint 3 x1 + x2 = .8
capture noisily frontier yf x1 x2, constraints(3)
capture noisily ereturn list
constraint 4 [lnsig2v]_cons = -2
capture noisily frontier yf x1 x2, constraints(4)
capture noisily frontier yf x1 x2, distribution(exponential) constraints(3)
capture noisily frontier yf x1 x2, uhet(x3) constraints(3)

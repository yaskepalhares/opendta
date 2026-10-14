* Phase 3d: random-number generator (mt64) and r*() functions.
* OpenDTA's sequences only match Stata's if the seed is
* turned into a state in the same way (VERIFY in core/rng.py).
clear
set seed 123
display %20.17f runiform()
display %20.17f runiform()
display %20.17f rnormal()
set seed 0
display %20.17f runiform()
set seed 2147483647
display %20.17f runiform()
set obs 5
set seed 42
gen double u = runiform()
gen double n = rnormal(10, 2)
gen k = runiformint(1, 6)
gen b = rbinomial(10, .5)
gen p = rpoisson(3)
list
display c(rng)
display c(rng_current)
capture noisily set seed -1

* Phase 3c: sample (simple random sampling without replacement): percentage, count, if,
* by() and the by prefix, rounding, use of the random-number stream and sort order
clear
set seed 310
set obs 40
gen id = _n
gen g = mod(_n, 3)
preserve
sample 25
count
list id
display runiform()
restore
preserve
sample 7, count
list id
display runiform()
restore
preserve
sample 50 if g == 0
tabulate g
restore
preserve
sample 50 in 1/10
list id in 1/12
restore
preserve
sample 2, count by(g)
sort id
list id g
restore
preserve
bysort g: sample 10
count
list id g
restore
preserve
sort g id
sample 50
describe, short
restore
preserve
keep in 1/37
sample 10
count
sample 100
sample 0
count
restore
capture noisily sample 120
capture noisily sample -1, count
capture noisily sample 2.5, count
capture noisily sample
display runiform()

* Phase 8: codebook, misstable and bootstrap
* (the bootstrap replicates use OpenDTA's resampling: the standard errors
* differ from Stata's; the case checks the layout and the observed estimates)
clear
set obs 74
gen mpg = 12 + mod(_n*7, 30) + (mod(_n, 4) == 0)*mod(_n, 9)
gen price = 3000 + mod(_n*131, 9000) + mod(_n, 7)^3
gen foreign = _n > 52
label define orig 0 "Domestic" 1 "Foreign"
label values foreign orig
label variable mpg "Mileage (mpg)"
gen rep = mod(_n, 5) + 1 if mod(_n, 11) != 0
gen peso = 2000 + 20*mpg + mod(_n*17, 300) if mod(_n, 13) != 0
gen str12 nome = "carro " + string(mod(_n, 20))
gen w = mpg/10

codebook mpg foreign rep nome w
codebook mpg price, compact
codebook, problems
misstable summarize
misstable summarize, all
return list
misstable patterns
misstable patterns, frequency

bootstrap, reps(20) seed(123): regress price mpg
bootstrap m=r(mean) sd=r(sd), reps(20) seed(1) nodots: summarize mpg
regress price mpg, vce(bootstrap, reps(20) seed(5) nodots)
estat bootstrap, all

* Fase 5a: regress (tabela ANOVA, coeficientes, vce, pesos, opções).
clear
estimates clear
set seed 2024
set obs 60
gen g = ceil(_n/20)
gen h = mod(_n, 2)
gen x1 = round(runiform()*10, .1)
gen x2 = round(runiform()*5, .01)
gen y = 3 + 2*x1 - 1.5*x2 + g + round(runiform()*4, .01)
gen w = 1 + mod(_n, 3)
replace x2 = . in 7
regress y x1 x2
regress
regress y x1 x2, level(90)
regress y x1 x2 if h == 1
regress y x1 x2 in 1/40
regress y x1 x2, noconstant
regress y x1 x2, vce(robust)
regress y x1 x2, robust
regress y x1 x2, vce(cluster g)
regress y x1 x2 [aweight=w]
regress y x1 x2 [fweight=w]
regress y x1 x2 [pweight=w]
regress y x1 x2, noheader
regress y x1
gen x3 = 2*x1
regress y x1 x3 x2
regress y x1 x2 if g < 3
ereturn list
display e(N) " " e(df_m) " " e(df_r) " " e(r2)
display _b[x1] " " _se[x1] " " _b[_cons]
matrix list e(b)
matrix list e(V)
count if e(sample)
regress y
capture noisily display _b[x1]
capture noisily regress y x9
capture noisily regress

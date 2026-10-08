* Fase 1a: criação de dados, tipos, generate/replace, describe, list, count.
clear
set obs 6
gen id = _n
gen x = id * 1.5
gen y = 1 if id > 2
gen byte b = 2.7
gen byte big = 200 in 1/3
gen float f = 0.1
gen double d = 0.1
count if f == 0.1
count if d == 0.1
count if f == float(0.1)
gen s = "obs" + string(id)
replace x = . in 3
replace x = x[_n-1] if missing(x)
replace b = 1000 in 1
replace b = 1.5 in 2
replace b = b
replace s = "observação longa" in 1
label variable x "Valor de x"
describe
list
list id x s in 1/3, noobs
list id x, clean
list id x if id > 3, sep(2)
display x[2] " " _N " " s[1]
count
count if missing(y)

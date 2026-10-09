* Fase 5b: tobit (censura à esquerda e nos dois lados).
clear
estimates clear
ereturn clear
set seed 1010
set obs 60
gen id = ceil(_n/6)
gen a = id/3
gen x1 = round(runiform()*10, .1) + a
gen x2 = round(runiform()*5, .01)
gen y = 1 + .5*x1 - x2 + 2*a + round(runiform()*2, .01)
gen ys = max(0, y - 5)
tobit ys x1 x2, ll(0)
ereturn list
tobit ys x1 x2, ll(0) ul(4)
predict xb
summarize xb

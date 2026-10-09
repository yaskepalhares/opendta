* Fase 3b: collapse, contract, expand, fillin.
clear
input g h x w
1 1 10 1
1 2 12 2
1 1 . 1
2 1 20 3
2 2 25 1
2 2 30 2
end
label variable x "Renda"
preserve
collapse (mean) mx=x (sd) sx=x (count) n=x (sum) s=x (min) lo=x (max) hi=x (p50) md=x, by(g)
list
describe
restore
preserve
collapse x [fw=w], by(g h)
list
restore
preserve
collapse (first) f=x (last) l=x (firstnm) fn=x (lastnm) ln=x (percent) pc=x, by(g)
list
restore
preserve
contract g h
list
restore
preserve
contract g h, freq(n) percent(pct) cfreq(cf) cpercent(cpct) zero
list
describe
restore
preserve
expand w
list g x w
expand 2 if g == 1, generate(novo)
list g x novo
restore
drop if g == 2 & h == 1
fillin g h
* o sort do Stata desempata ao acaso (set sortseed): ordena por x antes de listar
sort g h x
list

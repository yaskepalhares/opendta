* Fase 3b: append, merge, joinby, cross. Cria e apaga arquivos na pasta atual.
clear
input id str5 nome renda
1 "ana" 10
2 "bia" .
3 "caio" 30
5 "eva" 50
end
label define sim 1 "um"
label values id sim
capture erase odta_m.dta
save odta_m
clear
input id renda byte idade
1 11 20
2 22 30
4 44 40
end
label define sim 1 "one"
label values id sim
capture erase odta_u.dta
save odta_u
use odta_m, clear
merge 1:1 id using odta_u
list
describe
use odta_m, clear
merge 1:1 id using odta_u, update
list
use odta_m, clear
merge 1:1 id using odta_u, update replace
list
use odta_m, clear
merge 1:1 id using odta_u, keep(match) keepusing(idade) nogenerate
list
use odta_m, clear
merge 1:1 id using odta_u, generate(origem) noreport
tab origem
use odta_m, clear
capture noisily merge 1:1 id using odta_u, assert(match)
use odta_m, clear
append using odta_u, generate(fonte)
list
describe
clear
input byte v
1
end
capture erase odta_a.dta
save odta_a
clear
input v
1.5
end
append using odta_a
describe
clear
input k str1 a
1 "x"
1 "y"
2 "z"
end
capture erase odta_j.dta
save odta_j
clear
input k b
1 10
1 20
3 30
end
joinby k using odta_j
list
clear
input k b
1 10
3 30
end
joinby k using odta_j, unmatched(both)
list
clear
input c
1
2
end
cross using odta_j
list
clear
input id
1
1
end
capture noisily merge 1:1 id using odta_u
erase odta_m.dta
erase odta_u.dta
erase odta_a.dta
erase odta_j.dta

* Fase 3b: reshape long/wide.
clear
input id sexo inc80 inc81 inc82 ue80 ue81 ue82
1 0 5000 5500 6000 0 1 0
2 1 2000 2200 3300 1 0 0
3 0 3000 2000 1000 0 0 1
end
reshape long inc ue, i(id) j(ano)
list, sep(3)
describe
reshape wide
list
describe
reshape long inc, i(id) j(ano)
list
clear
input id str3 tipo valor
1 "a" 10
1 "b" 20
2 "a" 30
end
reshape wide valor, i(id) j(tipo) string
list
reshape long valor, i(id) j(tipo) string
list
clear
input id x1 x2
1 5 6
1 7 8
end
capture noisily reshape long x, i(id) j(t)
clear
input id t x z
1 1 5 0
1 2 6 1
end
capture noisily reshape wide x, i(id) j(t)

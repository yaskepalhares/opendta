* Fase 1a: sort, gsort, by/bysort, _n, _N, subscritos, sum().
clear
input g x
2 20
1 10
2 40
1 30
1 50
.  5
end
sort g x
list
gsort -x
list
gsort -g +x
list
bysort g (x): gen n = _n
by g: gen N = _N
by g: gen soma = sum(x)
by g: gen ant = x[_n-1]
by g: gen prim = x[1]
list
by g: count
capture noisily by x: gen z = 1
sort g
by g: keep if _n == 1
list

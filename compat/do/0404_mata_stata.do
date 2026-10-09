* Fase 4: interface Stata <-> Mata (st_*).
clear
set obs 4
gen x = _n
gen str3 nome = "n" + string(_n)
mata
X = st_data(., "x")
X
st_data((1, 3), 1)
N = st_sdata(., "nome")
N
st_store(., "x", X :* 10)
st_addvar("double", "y")
st_store(., "y", X :^ 2)
st_nobs(), st_nvar()
st_varname(2)
st_varindex("y")
st_local("loc", "valor")
st_global("glob", "g1")
st_numscalar("sc", 42)
st_matrix("M", (1, 2 \ 3, 4))
st_matrix("M")
stata("display 1 + 1")
end
list
display "`loc' $glob " sc
matrix list M
mata: st_numscalar("sc")
mata: st_local("loc")

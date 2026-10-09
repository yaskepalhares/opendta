* Fase 4: st_view (grava nos dados) e arrays associativos.
clear
set obs 3
gen x = _n
mata
V = .
st_view(V, ., "x")
V
V[2, 1] = 99
V[., 1] = V[., 1] :* 2
V
A = asarray_create()
asarray(A, "um", 1)
asarray(A, "dois", (2, 2))
asarray(A, "dois")
asarray_contains(A, "tres")
asarray_elements(A)
asarray_keys(A)
end
list

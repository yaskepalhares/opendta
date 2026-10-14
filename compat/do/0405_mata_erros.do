* Phase 4: errors in Mata. Does the "mata" block (without a colon) continue after
* an error? And "mata:"?
mata
nada
(1, 2) * (3, 4)
"a" + 1
A = (1, 2)
A[3]
real scalar f()
{
    return((1, 2) * (3, 4))
}
f()
real scalar f()
{
    return(1)
}
x = 1
x
end
display "depois do bloco mata"
capture noisily mata: nada2
display _rc
mata:
y = 1
nada3
y = 2
end
display "fim"

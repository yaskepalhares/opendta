* Laços, condicionais e o eco de blocos.
forvalues i = 1/3 {
    display "i = `i'"
}
forvalues i = 10(-5)0 {
    display `i'
}
forvalues i = 1 3 to 7 {
    display `i'
}
forvalues i = 0(.25)1 {
    display `i'
}
foreach x in a "b c" d {
    display "[`x']"
}
foreach n of numlist 1/3 10 {
    display `n'
}
local nomes ana beto
foreach p of local nomes {
    display "`p'"
}
local k 0
while `k' < 3 {
    local ++k
    display "k = `k'"
}
forvalues i = 1/5 {
    if `i' == 2 continue
    if `i' == 4 continue, break
    display `i'
}
local x 5
if `x' > 3 {
    display "maior"
}
else {
    display "menor"
}
if `x' > 10 {
    display "A"
}
else if `x' > 4 {
    display "B"
}
else {
    display "C"
}
if `x' == 5 display "linha única"
quietly {
    display "não aparece"
    noisily display "aparece"
}
display "fim"

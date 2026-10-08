* Fase 2: program, syntax, args, gettoken, return/ereturn, temporários, preserve.
capture program drop _all
clear
set obs 5
gen idade = _n * 10
gen renda = _n
replace renda = . in 2
program define media, rclass
    syntax varname [if] [in] [, Detail]
    marksample touse
    quietly count if `touse'
    local n = r(N)
    tempvar soma
    gen double `soma' = sum(`varlist') if `touse'
    local total = `soma'[_N]
    return scalar N = `n'
    return scalar mean = `total' / `n'
    return local var "`varlist'"
    if "`detail'" != "" {
        display "n = `n', total = `total'"
    }
end
media idade
return list
media renda, d
display r(mean)
capture noisily media
capture noisily media idade renda
capture noisily media idade, xyz
program define dois
    args a b
    display "a=`a' b=`b'"
end
dois um "dois tres"
local lista `"um "dois tres" (quatro) cinco"'
gettoken p lista : lista
gettoken q lista : lista
display "`p'|`q'|`lista'"
program define e1, eclass
    ereturn clear
    ereturn scalar N = 42
    ereturn local cmd "e1"
end
e1
ereturn list
preserve
drop in 1/3
count
restore
count
capture noisily restore
program define pr
    preserve
    keep in 1
    error 459
end
capture noisily pr
count
describe, short
program drop media dois e1 pr

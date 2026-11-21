* Probe (survey): c(changed), c(sortlist), long value-label names, number and date formats.
version 14
set more off
clear
set seed 5
set obs 50
gen id = _n
gen g = mod(_n,3)
gen x = round(invnormal(runiform()),.01)
gen str20 s = "abc"
label define rotulo_valor_muito_longo_30ch 0 "zero com texto bem comprido para testar quebra" 1 "um" 2 "dois"
label values g rotulo_valor_muito_longo_30ch
save "odta_probe_state.dta", replace
display "after save: changed=" c(changed)
use id x using "odta_probe_state.dta", clear
display "use varlist using: changed=" c(changed)
use "odta_probe_state.dta", clear
display "use: changed=" c(changed)
capture noisily display "c(sortlist)=[" c(sortlist) "]"
display "  rc=" _rc
sort g g x
display "sort g g x: changed=" c(changed)
describe, short
save "odta_probe_state.dta", replace
compress
display "compress (nothing to do): changed=" c(changed)
keep if id > 0
display "keep (all kept): changed=" c(changed)
save "odta_probe_state.dta", replace
preserve
replace x = 1
restore
display "preserve/restore: changed=" c(changed)
merge 1:1 id using "odta_probe_state.dta"
display "merge: changed=" c(changed)
drop _merge
save "odta_probe_state.dta", replace
append using "odta_probe_state.dta"
display "append: changed=" c(changed)
use "odta_probe_state.dta", clear
describe
tab g
tab g, missing
tabstat x, by(g) stat(mean sd n)
list in 1/3
label list
set linesize 255
tab g s
gen double big = 123456789012 * _n
format big %15.0fc
list big in 1/3
capture noisily gen m = ym(2020,1) + _n
capture noisily format m %tm
capture noisily list m in 1/2
gen d = mdy(1,1,2020) + _n
format d %td
list d in 1/2
format d %tdDay,_DD_Mon_CCYY
list d in 1/2
format d %tdnn/dd/CCYY
list d in 1/2
erase "odta_probe_state.dta"

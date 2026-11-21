* Probe (survey): commands that must fail with a Stata error code (not abort) under capture.
* In OpenDTA several of these raised Python exceptions. Each block is independent.
version 14
set more off
clear
display "=== bare commands that need arguments"
foreach c in ksmirnov ladder prtest runtest sdtest {
    capture noisily `c'
    display "  rc=" _rc " for `c'"
}
display "=== file operations on a folder that does not exist"
set obs 5
gen x = _n
capture noisily copy "do/2006_probe_errors_in_files.do" "pasta_inexistente_xyz/b.txt", replace
display "  rc=" _rc " copy"
capture noisily file open fh using "pasta_inexistente_xyz/a.txt", write replace
display "  rc=" _rc " file open write"
capture noisily outfile x using "pasta_inexistente_xyz/a.txt", replace
display "  rc=" _rc " outfile"
capture noisily log using "pasta_inexistente_xyz/x.log", replace name(probe)
display "  rc=" _rc " log using"
capture noisily save "pasta_inexistente_xyz/a.dta", replace
display "  rc=" _rc " save"
capture noisily export delimited "pasta_inexistente_xyz/a.csv", replace
display "  rc=" _rc " export delimited"
capture noisily graph twoway scatter x x
capture noisily graph export "pasta_inexistente_xyz/g.png", replace
display "  rc=" _rc " graph export"
capture noisily postfile pf x using "pasta_inexistente_xyz/p.dta", replace
display "  rc=" _rc " postfile"
capture noisily use "pasta_inexistente_xyz/a.dta", clear
display "  rc=" _rc " use"
capture noisily mkdir "pasta_inexistente_xyz/a/b"
display "  rc=" _rc " mkdir"

* Differences seen in a survey do-file (Stata 14): memory, S_DATE/S_TIME, c(flavor), log using in a
* missing folder, use varlist using, long value label names in describe.
clear all
set more off
capture noisily log using "/nao/existe/x.log", text replace
display _rc
display c(flavor)
display length("$S_DATE") > 0
display length("$S_TIME") == 8

set obs 5
gen byte a = _n
gen byte b = _n
gen byte c = _n
sort a b c
save "_t2001.dta", replace
use a c using "_t2001.dta", clear
describe
display c(changed)
use a c using "_t2001.dta" if a > 2, clear
display c(changed)
use "_t2001.dta", clear
label define l10chars__ 1 "x"
label define l11chars___ 1 "x"
label define l15chars_______ 1 "x"
label variable a "has label"
gen byte d = 1
gen byte e = 1
gen byte f = 1
label values a l10chars__
label values d l11chars___
label values e l15chars_______
label variable d "second"
label variable e "third"
describe
erase "_t2001.dta"

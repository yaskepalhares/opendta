* Data management, second part: table with a supercolumn variable and by(), list with string()
* and display, import sasxport with describe, member() and the value labels of formats.xpf, and
* export sasxport with rename and vallabfile(). Every command whose behavior the manual does not
* settle runs under capture noisily, so that one Stata run shows the whole case.
clear
set obs 48
gen a = 1 + mod(_n, 2)
gen b = 1 + mod(_n*5, 3)
gen c = 1 + mod(_n*7, 4)
gen d = 1 + mod(_n*11, 3)
gen y = mod(_n*13, 17) + .25
label define al 1 "Male" 2 "Female"
label values a al
label variable b "Group"

* table with three variables and with by()
capture noisily table a b c
capture noisily table a b c, c(mean y) format(%9.2f)
capture noisily table a b c, col scol row
capture noisily table a b c, c(freq mean y sd y) format(%9.3f)
capture noisily table a b, by(c)
capture noisily table a b, by(c d) row col
capture noisily table a, by(c) c(freq mean y)
capture noisily table a b c, by(d) c(n y)
capture noisily table a b c d
capture noisily table a, by(b c d)

* list with string() and display
clear
input str24 make price mpg str10 origin
"AMC Concord" 4099 22 "Domestic"
"Buick LeSabre Custom" 5788 18 "Domestic"
"Honda Accord" 5799 25 "Foreign"
end
capture noisily list, string(5)
capture noisily list make origin, string(3) clean
capture noisily list, display
capture noisily list in 1, display noobs
capture noisily list in 2, display divider
set linesize 40
capture noisily list
capture noisily list, table
set linesize 255

* export sasxport with value labels, rename and vallabfile(); import sasxport with describe,
* member() and formats.xpf
clear
input id sex agegroup longvariablename
1 1 1 5.5
2 2 3 .
3 1 2 7
4 2 1 8.25
end
label define sexlab 1 "Female" 2 "Male"
label define agelab 1 "Young" 2 "Middle" 3 "Old"
label values sex sexlab
label values agegroup agelab
label variable longvariablename "A long name"
capture erase odta_s.xpt
capture erase formats.xpf
capture erase odta_s.sas
capture noisily export sasxport odta_s
capture noisily export sasxport odta_s, rename
capture noisily export sasxport odta_s, rename vallabfile(both) replace
capture noisily export sasxport odta_s, rename vallabfile(sascode) replace
capture noisily export sasxport odta_s, rename vallabfile(none) replace
capture noisily export sasxport odta_s, rename vallabfile(xpf) replace
capture noisily type odta_s.sas
capture noisily import sasxport odta_s, describe
capture noisily return list
capture noisily import sasxport odta_s, describe member(odta_s)
capture noisily import sasxport odta_s, clear
capture noisily describe
capture noisily list
capture noisily label list
capture noisily import sasxport odta_s, clear novallabels
capture noisily describe
capture noisily import sasxport odta_s, clear member(odta_s)
capture noisily describe
capture noisily import sasxport odta_s, clear member(nothere)
capture noisily import sasxport formats.xpf, clear
capture noisily describe
capture noisily list
capture erase odta_s.xpt
capture erase formats.xpf
capture erase odta_s.sas

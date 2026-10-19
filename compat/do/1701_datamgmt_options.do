* Data management options: file with binary read/write, seek and set byteorder; infile with
* byvariable(); putexcel to .xls; describe numbers, varlist, replace and using; list header(#) and
* nodotz; import and export sasxport (import sas and import spss are newer than Stata 14 and have
* pytest cases instead). Every command whose behavior the manual does not settle runs under
* capture noisily, so that one Stata run shows the whole case.
clear
capture erase odta_b.bin
capture noisily file open h using odta_b.bin, write binary
capture noisily file write h %14s "mymatout 1.0.1" %1b (byteorder())
capture noisily file write h %2b (3) %4b (-5) %8z (3.14159) %4z (.) %1b (.a)
capture noisily file write h %1b (2137) %1bs (2137) %1bu (2137) %2b (124.75) %3s "abcdef" %8s "ab"
capture noisily file seek h query
capture noisily return list
capture noisily file close h
capture noisily file open h using odta_b.bin, read binary
capture noisily file read h %14s sig
capture noisily display "`sig'"
capture noisily file read h %1b bo
local border = bo
capture noisily file set h byteorder `border'
capture noisily file read h %2b v
capture noisily display v
capture noisily file read h %4b v
capture noisily display v
capture noisily file read h %8z v
capture noisily display %9.5f v
capture noisily file read h %4z v
capture noisily display v
capture noisily file read h %1b v
capture noisily display v == .a
capture noisily file read h %1b v
capture noisily display v == .
capture noisily file read h %1bs v
capture noisily display v
capture noisily file read h %1bu v
capture noisily display v
capture noisily file read h %2b v
capture noisily display v
capture noisily file read h %3s t
capture noisily display "`t'"
capture noisily file read h %8s t
capture noisily display "`t'"
capture noisily return list
capture noisily file read h %2b v
capture noisily return list
capture noisily file seek h tof
capture noisily file read h %2s t
capture noisily display "`t'"
capture noisily file seek h 3
capture noisily file seek h eof
capture noisily file query
capture noisily file close h
capture noisily file read h %2b v
capture erase odta_b.bin
capture erase odta_rw.bin
capture noisily file open h using odta_rw.bin, read write binary
capture noisily file write h %4b (7) %2b (9)
capture noisily file seek h tof
capture noisily file read h %4b v
capture noisily display v
capture noisily file close h
capture erase odta_rw.bin

* infile with byvariable(): the file holds all the values of a variable, then of the next one
capture erase odta_time.raw
capture erase odta_time2.raw
capture noisily file open h using odta_time.raw, write text
capture noisily file write h "1980 1981 1982 1983" _n "14 17 25 30" _n "120 135 150" _n "180" _n
capture noisily file close h
capture noisily file open h using odta_time2.raw, write text
capture noisily file write h "1980 1981 1982 1983;" _n "14 17 25 ;" _n "120 135 150 180;" _n
capture noisily file close h
capture noisily infile year amount cost using odta_time, byvariable(4) clear
capture noisily list
capture noisily infile year amount cost using odta_time2, byvariable(100) clear
capture noisily list
capture noisily infile year amount cost using odta_time, byvariable(x) clear
capture noisily infile year amount cost using odta_time, byvariable(2) clear
capture noisily list
capture erase odta_time.raw
capture erase odta_time2.raw

* describe: numbers, varlist, replace and using
clear
input str8 make price mpg
"AMC" 4099 22
"Buick" 5788 18
"Cad" 7827 .
end
label variable price "Price (USD)"
label define yn 1 "one"
label values mpg yn
sort price
capture noisily describe, numbers
capture noisily describe price mpg, numbers
capture noisily describe, varlist
capture noisily return list
capture erase odta_d.dta
save odta_d, replace
clear
capture noisily describe using odta_d
capture noisily describe using odta_d, short
capture noisily describe using odta_d, simple
capture noisily describe make using odta_d, varlist
capture noisily return list
use odta_d, clear
capture noisily describe, clear
capture noisily describe, replace
capture noisily list
capture noisily char list _dta[]
capture noisily return list
capture noisily describe, replace
use odta_d, clear
capture noisily describe make price, replace
capture noisily list

* list: header(#) and nodotz
clear
input x y
1 2
3 .a
5 .
7 8
end
capture noisily list, header(2)
capture noisily list, clean header(2)
capture noisily list, nodotz
capture noisily list, header
capture noisily list, noheader header
capture noisily list, mean
capture noisily list, sum(x) N labvar(y)
capture noisily list, clean mean noobs
capture noisily list, nocompress
capture noisily list, labvar(zz) sum

* putexcel to .xls
capture erase odta_p.xls
capture noisily putexcel A1=("id") B1=("score") A2=(1) B2=(2.5) using odta_p.xls, replace
capture noisily putexcel set odta_p.xls, modify
capture noisily putexcel C1=("n") C2=(3)
capture noisily import excel using odta_p.xls, clear
capture noisily list
capture erase odta_p.xls

* import and export sasxport
clear
input id str6 name weight
1 "ana" 60.5
2 "bruno" .
3 "carla" 72
end
label variable weight "Weight (kg)"
capture erase odta_x.xpt
capture noisily export sasxport odta_x
capture noisily export sasxport odta_x
capture noisily export sasxport odta_x, replace
clear
capture noisily import sasxport odta_x
capture noisily import sasxport odta_x, clear
capture noisily describe
capture noisily list
capture noisily return list
capture erase odta_x.xpt
capture erase odta_d.dta

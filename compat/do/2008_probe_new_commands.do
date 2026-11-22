* Probe: the functions and commands added after compat 2005, and what the manual does not settle
* (layouts, limits, r() and e() results). Each section runs under capture noisily.
version 14
set more off
clear

display "=== date and time functions"
display %21.0f Cofc(1893456000000) " " %21.0f cofC(1893456027000) " " %21.0f Clock("01jan2020 00:00:00","DMYhms")
display %21.0f clock("5 March 2020 3:04:05.5 pm","DMYhms") " " %td dofw(3121) " " %td dofh(121)
display week(mdy(12,31,2020)) " " week(mdy(12,24,2020)) " " wofd(mdy(12,31,2020))
display date("1/2/99","MDY",2020) " " date("1/2/99","MD19Y") " " date("20200105","YMD") " " date("1/2/99","MDY")
display monthly("2020-Mar","YM") " " quarterly("2020q3","YQ") " " halfyearly("2020h2","YH") " " weekly("2020w2","YW")
display hh(clock("10:30:15","hms")) " " mm(clock("10:30:15","hms")) " " ss(clock("10:30:15","hms"))
display %21.0f hms(1,2,3) " " %21.0f Cmdyhms(1,1,2017,0,0,0) " " dofC(Cofd(21915))

display "=== other functions"
capture noisily display "f1: [" strtoname("") "]"
display "  rc=" _rc
capture noisily display "f2: [" strtoname("1x",0) "]"
display "  rc=" _rc
capture noisily display "f3: [" soundex("Ashcraft") "]"
display "  rc=" _rc
capture noisily display "f4: [" soundex_nara("Ashcraft") "]"
display "  rc=" _rc
capture noisily display "f5: [" plural(2,"glasses","-es") "]"
display "  rc=" _rc
capture noisily display "f6: [" plural(0,"box","+es") "]"
display "  rc=" _rc
capture noisily display "f7: [" ustrword("a b c",-1) "]"
display "  rc=" _rc
capture noisily display "f8: [" ustrright("abc",2) "]"
display "  rc=" _rc
capture noisily display "f9: [" uchar(231) "]"
display "  rc=" _rc
capture noisily display "f10: [" ustrregexra("abab","b","x") "]"
display "  rc=" _rc
capture noisily display "f11: [" autocode(11,4,0,10) "]"
display "  rc=" _rc
capture noisily display "f12: [" autocode(-1,4,0,10) "]"
display "  rc=" _rc
capture noisily display "f13: [" irecode(9,1,3,5) "]"
display "  rc=" _rc
capture noisily display "f14: [" minfloat() "]"
display "  rc=" _rc
capture noisily display "f15: [" maxfloat() "]"
display "  rc=" _rc
capture noisily display "f16: [" maxlong() "]"
display "  rc=" _rc
capture noisily display "f17: [" fileread("nada.xyz") "]"
display "  rc=" _rc
capture noisily display "f18: [" chi2den(2,1) "]"
display "  rc=" _rc
capture noisily display "f19: [" nchi2tail(1,1,1) "]"
display "  rc=" _rc
capture noisily display "f20: [" invnchi2tail(1,1,.5) "]"
display "  rc=" _rc
capture noisily display "f21: [" npnchi2(3,4,.2) "]"
display "  rc=" _rc
capture noisily display "f22: [" ntden(5,1,0) "]"
display "  rc=" _rc
capture noisily display "f23: [" nt(5,1,0) "]"
display "  rc=" _rc
capture noisily display "f24: [" nttail(5,1,0) "]"
display "  rc=" _rc
capture noisily display "f25: [" invnttail(5,1,.5) "]"
display "  rc=" _rc
capture noisily display "f26: [" Fden(2,3,1) "]"
display "  rc=" _rc
capture noisily display "f27: [" nFden(2,3,1,1) "]"
display "  rc=" _rc
capture noisily display "f28: [" nF(2,3,1,1) "]"
display "  rc=" _rc
capture noisily display "f29: [" nFtail(2,3,1,1) "]"
display "  rc=" _rc
capture noisily display "f30: [" invnFtail(2,3,1,.5) "]"
display "  rc=" _rc
capture noisily display "f31: [" gammaden(2,1,0,1) "]"
display "  rc=" _rc
capture noisily display "f32: [" gammaptail(2,1) "]"
display "  rc=" _rc
capture noisily display "f33: [" invgammap(2,.3) "]"
display "  rc=" _rc
capture noisily display "f34: [" invgammaptail(2,.3) "]"
display "  rc=" _rc
capture noisily display "f35: [" ibetatail(2,3,.4) "]"
display "  rc=" _rc
capture noisily display "f36: [" invibetatail(2,3,.4) "]"
display "  rc=" _rc
capture noisily display "f37: [" lnnormalden(0) "]"
display "  rc=" _rc
capture noisily display "f38: [" binormal(0,0,.5) "]"
display "  rc=" _rc
capture noisily display "f39: [" hypergeometric(10,5,3,2) "]"
display "  rc=" _rc
capture noisily display "f40: [" nbinomialp(5,2,.5) "]"
display "  rc=" _rc
capture noisily display "f41: [" nbinomialtail(5,2,.5) "]"
display "  rc=" _rc
capture noisily display "f42: [" invnbinomial(5,2,.3) "]"
display "  rc=" _rc
capture noisily display "f43: [" invnbinomialtail(5,2,.3) "]"
display "  rc=" _rc
capture noisily display "f44: [" invbinomial(10,4,.3) "]"
display "  rc=" _rc
capture noisily display "f45: [" invbinomialtail(10,4,.3) "]"
display "  rc=" _rc
capture noisily display "f46: [" invpoisson(3,.4) "]"
display "  rc=" _rc
capture noisily display "f47: [" invpoissontail(3,.4) "]"
display "  rc=" _rc
capture noisily display "f48: [" cloglog(.5) "]"
display "  rc=" _rc
capture noisily display "f49: [" invcloglog(0) "]"
display "  rc=" _rc
capture noisily display "f50: [" logistic(1) "]"
display "  rc=" _rc
capture noisily display "f51: [" ustrrpos("abcb","b") "]"
display "  rc=" _rc

display "=== creturn strings near the cut"
local s18 = "abcdefghijklmnopqr"
local s19 = "abcdefghijklmnopqrs"
local s20 = "abcdefghijklmnopqrst"
foreach v in s18 s19 s20 {
    capture noisily set httpproxyhost ``v''
}
capture noisily creturn list

display "=== query"
foreach c in output interface graphics efficiency network update trace mata unicode other {
    capture noisily query `c'
}
capture noisily set_defaults output
display "rc=" _rc
capture noisily set_defaults foo
display "rc=" _rc

display "=== which command widens the format of a labeled variable"
clear
set obs 6
gen g = mod(_n,3)
gen id = _n
label define lg 0 "zero com texto bem comprido para testar quebra" 1 "um" 2 "dois"
describe g
label values g lg
describe g
compress
describe g
save "odta_probe_fmt.dta", replace
use "odta_probe_fmt.dta", clear
describe g
merge 1:1 id using "odta_probe_fmt.dta", nogenerate
describe g
append using "odta_probe_fmt.dta"
describe g

display "=== tabulate label limits"
clear
set obs 4
gen g = _n
label define lt 1 "abcdefghij abcdefghij abcdefghij abcdefghi" 2 "abcdefghij abcdefghij abcdefghij abcdefghij" ///
    3 "abcdefghijabcdefghijabcdefghijabcdefghijab" 4 "abcdefghij abcdefghij ab"
label values g lt
gen s = "x"
gen id = _n
tab g
tab g s
set linesize 80
tab g s
set linesize 255
tabstat g, by(g) stat(mean n)
tabstat g id, by(g) stat(mean)
capture noisily tabstat g, by(g) stat(mean) labelwidth(30)

display "=== input"
clear
input x
.5
-.25
.a
end
list

display "=== small data commands"
clear
set obs 3
gen c = 1
gen B = 2
gen a1 = 3
gen a10 = 4
gen a2 = 5
aorder
describe, simple
renpfix a z
return list
capture noisily renpfix q
display "rc=" _rc
capture noisily move c
display "rc=" _rc
label define yn 0 "No" 1 "Yes"
numlabel yn, add
label list yn
numlabel yn, add
label list yn
numlabel yn, add force
label list yn
numlabel yn, remove
label list yn
capture noisily numlabel yn, add remove
display "rc=" _rc
capture noisily numlabel, add mask("[#] ")
label list

display "=== mkmat, svmat, numlist"
clear
input a b str5 id
1 2 "x y"
3 . "z"
5 6 "w"
end
mkmat a b, matrix(A) rownames(id)
matrix list A
mkmat a b, matrix(B) obs
matrix list B
capture noisily mkmat a b, matrix(C) nomissing rowprefix(p)
matrix list C
matrix D = (1, 2 \ 3, 4)
matrix colnames D = u v
matrix coleq D = e1 e2
svmat D, names(eqcol)
svmat D, names(matcol)
svmat D, names(col)
svmat D, names(k)
describe, simple
foreach n in "1/3" "3 1 2" "1(2)7" "1 1 2" ".5 1.5" {
    numlist "`n'"
    display "`n': [" r(numlist) "]"
}
capture noisily numlist "1 1 2", ascending
display "rc=" _rc
capture noisily numlist "1/2000"
display "rc=" _rc
numlist "3 1 2", sort
display r(numlist)

display "=== rcof"
capture noisily rcof "noisily display 1/" == 198
display "rc=" _rc
capture noisily rcof "display 1/" == 0
display "rc=" _rc

display "=== estat vif"
clear
set obs 30
gen x1 = mod(_n * 7, 11)
gen x2 = x1 + mod(_n, 4)
gen x3 = mod(_n * 3, 5)
gen y = x1 + 2 * x2 - x3 + mod(_n * 13, 7)
quietly regress y x1 x2 x3
capture noisily estat vif
return list
capture noisily vif
quietly regress y x1 x2 x3, noconstant
capture noisily estat vif
display "rc=" _rc
capture noisily estat vif, uncentered
gen yb = y > 15
quietly logit yb x1
capture noisily estat vif
display "rc=" _rc

display "=== fvset"
clear
set obs 30
gen g = mod(_n, 3) + 1
gen y = g * 2 + mod(_n * 7, 5)
fvset base 2 g
fvset design asbalanced g
capture noisily fvset report
regress y i.g
fvset base frequent g
fvset report g
fvset clear g
capture noisily fvset report
capture noisily fvset base bad g
display "rc=" _rc

display "=== sdtesti, robvar, sfrancia"
sdtesti 20 5 2 1.5
return list
sdtesti 20 5 2 15 4 3
return list
clear
set obs 40
gen g = mod(_n, 3)
gen x = mod(_n * 7, 13) * (1 + g)
robvar x, by(g)
return list
sfrancia x
return list
sfrancia x, noties
capture noisily sfrancia x if _n < 5
display "rc=" _rc

display "=== hexdump"
capture erase "odta_probe_hex.txt"
file open fh using "odta_probe_hex.txt", write replace
file write fh "Ab 1" _n "x" _tab "y" _n
file close fh
capture noisily hexdump "odta_probe_hex.txt"
capture noisily hexdump "odta_probe_hex.txt", analyze
return list
capture noisily hexdump "odta_probe_hex.txt", results
capture noisily hexdump "nada.xyz"
display "rc=" _rc

display "=== quantile regressions (headers; the bootstrap standard errors differ)"
clear
set obs 60
gen x = mod(_n*7,13)
gen z = mod(_n*5,7)
gen y = 2 + x + .5*z + mod(_n*11,9)
set seed 1
bsqreg y x z, reps(10)
ereturn list
iqreg y x z, reps(10) nodots
ereturn list
sqreg y x z, q(.25 .5 .75) reps(10) nodots
ereturn list
sqreg y x z, q(.1 .333) reps(10) nodots

display "=== pk"
clear
input id t c
1 0 0
1 .5 3.1
1 1 4.8
1 2 4.2
1 4 3.0
1 8 1.6
1 12 0.8
2 0 0
2 .5 2.5
2 1 4.1
2 2 4.5
2 4 3.4
2 8 1.9
2 12 1.1
3 0 0
3 .5 2.9
3 1 5.2
3 2 4.0
3 4 2.6
3 8 1.2
3 12 .5
end
pkexamine t c if id == 1
return list
pkexamine t c if id == 1, trapezoid fit(4)
return list
pkexamine t c if id == 1 & t > 0
pksumm id t c
pksumm id t c, stat(auc)
preserve
pkcollapse t c, id(id)
describe
list
restore
preserve
pkcollapse t c, id(id) stat(auc half cmax) keep(id)
list
restore
clear
input id seq p1 p2
1 1 10 12
2 1 11 13
3 1 9 11
4 2 12 10
5 2 13 12
6 2 11 10
end
pkshape id seq p1 p2, order(ab ba)
describe
list
label list
pkcross outcome
pkcross outcome, carryover(none)
pkequiv outcome treat period seq id
return list
pkequiv outcome treat period seq id, limit(.25) level(95)

display "=== forecast"
clear
set obs 30
gen t = _n
tsset t
gen x = mod(_n*7, 11)
gen y = 5 if _n == 1
replace y = 2 + .5*L.y + .3*x + mod(_n*3,5)/10 if _n > 1
replace y = . if _n > 25
quietly regress y L.y x
estimates store m1
capture noisily forecast query
capture noisily forecast describe
display "rc=" _rc
forecast create mymodel
capture noisily forecast
display "rc=" _rc
forecast estimates m1
forecast identity z = y + x, generate
forecast exogenous x
forecast describe
forecast list
forecast query
forecast solve
list t y f_y f_z in 24/30
forecast adjust y = y + 1
forecast solve, prefix(a_) log(off)
forecast solve, prefix(s_) static begin(20)
forecast solve, suffix(_p) periods(2)
list t a_y s_y y_p in 19/28
forecast drop
describe, simple
capture noisily forecast create other
display "rc=" _rc
forecast clear
forecast query

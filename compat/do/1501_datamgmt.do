* Every command that Stata 14 may not have, or whose behavior the manual does not settle, runs under
* capture noisily, so that one Stata run shows the whole case.
* Data management: clonevar, ds, lookfor, insobs, range, separate, ipolate,
* expandcl, stack, xpose, corr2data, inspect, labelbook, compare, cf, snapshot,
* filefilter, outfile, statsby, putmata, icd9 and bcal (layouts, messages, r())
clear
input id grp x y
1 1 1 10
2 1 2 .
3 1 4 30
4 2 1 5
5 2 3 .
6 2 5 25
end
label variable x "Dose given"
label define yn 1 "one" 2 "two"
label values grp yn
capture noisily clonevar xc = x if id < 4
capture noisily describe x xc
capture noisily ds
capture noisily return list
capture noisily ds, has(type numeric) not(varlabel)
ds, has(type numeric)
ds, not(varlabel)
capture noisily lookfor dose
capture noisily return list
capture noisily inspect x
capture noisily inspect y
capture noisily return list
capture noisily labelbook yn
labelbook, problems
capture noisily compare x xc
capture noisily by grp: ipolate y x, gen(yi)
capture noisily bysort grp: ipolate y x, gen(yi)
capture noisily list
capture noisily separate x, by(grp)
capture noisily describe x1 x2
capture noisily separate y, by(grp) generate(yy) shortlabel
preserve
capture noisily expandcl 2, cluster(grp) generate(newcl)
capture noisily sort grp id newcl
capture noisily list id grp newcl
restore
capture snapshot erase _all
capture noisily snapshot save, label("first")
capture noisily snapshot list
capture noisily drop if grp == 2
capture noisily snapshot restore 1
capture noisily count
capture erase odta_a.dta
capture noisily save odta_a, replace
replace x = 99 in 2
capture noisily cf _all using odta_a
capture noisily cf x y using odta_a, verbose
capture noisily statsby m=r(mean) n=r(N), by(grp) clear nodots: summarize x
capture noisily list
capture noisily putmata X=(grp m), replace
capture noisily mata: X
clear
set obs 3
capture noisily insobs 2, before(2)
capture noisily count
capture noisily range t 0 1 5
capture noisily list
clear
set obs 3
gen a = _n
gen b = 10*_n
capture noisily stack a b, into(z) clear
capture noisily list
clear
set obs 2
gen a = _n
gen b = 10*_n
xpose, clear varname
capture noisily list
clear
matrix C = (1, .3 \ .3, 1)
capture noisily corr2data u v, n(100) corr(C) means(5 10) sds(2 3)
capture noisily tabstat u v, statistics(mean sd)
capture noisily correlate u v
clear
input id x
1 1.5
2 2.25
end
capture erase odta_out.raw
capture erase odta_f.raw
capture noisily outfile using odta_out.raw, replace
capture noisily type odta_out.raw
capture noisily outfile using odta_out.raw, comma replace
capture noisily type odta_out.raw
capture noisily filefilter odta_out.raw odta_f.raw, from(",") to(";") replace
capture noisily type odta_f.raw
clear
input str8 dx
"250.00"
"25000"
"401.9"
"V27.0"
"E880.0"
"bad!"
"999.99"
end
capture noisily icd9 check dx
capture noisily return list
capture noisily icd9 check dx, list
capture noisily icd9 clean dx, dots generate(c)
gen c = dx
capture noisily icd9 clean c, dots
capture noisily drop if dx == "bad!"
capture noisily icd9 clean c, dots
capture noisily icd9 generate cat = c, category
capture noisily list
erase odta_a.dta
erase odta_out.raw
erase odta_f.raw

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
clonevar xc = x if id < 4
describe x xc
ds
return list
ds, has(type numeric) not(varlabel)
lookfor dose
return list
inspect x
inspect y
return list
labelbook yn
labelbook, problems
compare x xc
by grp: ipolate y x, gen(yi)
list
separate x, by(grp)
describe x1 x2
separate y, by(grp) generate(yy) shortlabel
preserve
expandcl 2, cluster(grp) generate(newcl)
list id grp newcl
restore
snapshot save, label("first")
snapshot list
drop if grp == 2
snapshot restore 1
count
save odta_a, replace
replace x = 99 in 2
capture noisily cf _all using odta_a
capture noisily cf x y using odta_a, verbose
statsby m=r(mean) n=r(N), by(grp) clear nodots: summarize x
list
putmata X=(grp m)
mata: X
clear
set obs 3
insobs 2, before(2)
count
range t 0 1 5
list
clear
set obs 3
gen a = _n
gen b = 10*_n
stack a b, into(z) clear
list
clear
set obs 2
gen a = _n
gen b = 10*_n
xpose, clear varname
list
clear
matrix C = (1, .3 \ .3, 1)
corr2data u v, n(100) corr(C) means(5 10) sds(2 3)
summarize u v
correlate u v
clear
input id x
1 1.5
2 2.25
end
outfile using odta_out.raw, replace
type odta_out.raw
outfile using odta_out.raw, comma replace
type odta_out.raw
filefilter odta_out.raw odta_f.raw, from(",") to(";") replace
type odta_f.raw
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
icd9 check dx
return list
icd9 check dx, list
icd9 clean dx, dots generate(c)
icd9 generate cat = dx, category
list
erase odta_a.dta
erase odta_out.raw
erase odta_f.raw

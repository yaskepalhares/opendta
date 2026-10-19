* Programming: unab, tsunab, fvunab, unabcmd, fvexpand, fvrevar, tsrevar, _rmcoll, mark, markin,
* markout, postfile, matlist, tabdisp, program sortpreserve, plugin, pause, sleep,
* signestimationsample, makecns and matcproc, viewsource. The dialogs of window (stopbox, fopen,
* fsave) and pause on are left out: they wait for someone at the screen. Every command whose
* behavior the manual does not settle runs under capture noisily, so that one Stata run shows the
* whole case.
clear
set obs 74
gen rep78 = 1+mod(_n,5)
gen foreign = mod(_n,2)
gen mpg = _n
gen turn = 30+_n
gen trunk = 10 + mod(_n*7,13)
gen tt = turn + trunk
gen t = _n
tsset t

* unab, tsunab, fvunab, unabcmd
capture noisily unab x : mp tur, name(myopt())
capture noisily display "`x'"
capture noisily unab x : junk
capture noisily unab x : mpg turn trunk, max(1) name(myopt())
capture noisily unab x : mpg, min(3)
capture noisily tsunab mylist : l(1/3).mpg
capture noisily display "`mylist'"
capture noisily tsunab mylist : l(1/2).(mpg turn)
capture noisily display "`mylist'"
capture noisily fvunab fl : i.rep78 mpg c.mpg#c.turn rep78#foreign
capture noisily display "`fl'"
capture noisily unabcmd gen
capture noisily return list
capture noisily unabcmd ka

capture noisily unab x : mp tu, min(5)
display _rc
capture noisily unab x : mpg turn, max(1)
display _rc
capture noisily unab x : mpg, min(3) name(myopt())
capture noisily unab x : mpg, min(2) max(2) name(myopt())
capture noisily unab x : mpg turn trunk, max(2) name(myopt())
capture noisily unab x : mpg turn trunk, min(2) max(2)
capture noisily unab x : mpg turn trunk, max(1) name(myopt())
display _rc

* fvexpand, fvrevar, tsrevar
capture noisily fvexpand i.rep78 mpg
capture noisily return list
capture noisily fvexpand mpg i.foreign
capture noisily fvexpand c.turn#i(2,3).rep78
capture noisily fvexpand c.mpg#i.foreign i.rep78#i.foreign
capture noisily fvrevar i.rep78 mpg turn
capture noisily display "`r(varlist)'"
capture noisily describe
capture noisily fvrevar i(2,3)bn.rep78 mpg
capture noisily display "`r(varlist)'"
capture noisily fvrevar i(2,3)bn.rep78 mpg, list
capture noisily display "`r(varlist)'"
capture noisily fvrevar c.turn#i(2,3).rep78 L.mpg
capture noisily display "`r(varlist)'"
capture noisily fvrevar c.turn#i(2,3).rep78 L.mpg, tsonly
capture noisily display "`r(varlist)'"
capture noisily tsrevar l.mpg d.mpg turn
capture noisily display "`r(varlist)'"
capture noisily tsrevar l.mpg d.mpg turn, list
capture noisily display "`r(varlist)'"
capture noisily list mpg `=word("`r(varlist)'",1)' in 1/4

capture noisily fvexpand mpg turn
capture noisily return list
capture noisily fvexpand L.mpg
capture noisily return list
capture noisily fvexpand i.foreign L.mpg
capture noisily return list
capture noisily fvrevar turn i.foreign mpg i.rep78, list
capture noisily display "`r(varlist)'"
capture noisily fvrevar i.foreign mpg trunk i.rep78 turn, list
capture noisily display "`r(varlist)'"
capture noisily fvrevar turn mpg, list
capture noisily display "`r(varlist)'"
capture noisily tsrevar turn mpg l.mpg, list
capture noisily display "`r(varlist)'"

* _rmcoll
capture noisily _rmcoll turn trunk tt
capture noisily display r(varlist)
capture noisily _rmcoll i.rep78
capture noisily display r(varlist)
capture noisily _rmcoll i.rep78, expand
capture noisily display r(varlist)
capture noisily _rmcoll rep78#foreign
capture noisily display r(varlist)
capture noisily _rmcoll turn trunk tt if turn > 50
capture noisily _rmcoll mpg junk

* mark, markin, markout
clear
set obs 10
gen x = _n
replace x = . in 3
replace x = .a in 4
gen w = cond(_n==2, 0, 1)
capture noisily mark m1 if x > 2
capture noisily list x m1, noobs
capture noisily mark m2 [fweight=w]
capture noisily list w m2 in 1/3, noobs
capture noisily mark m3 [fweight=w], zeroweight
capture noisily list w m3 in 1/3, noobs
capture noisily markin if x>5, name(rng)
capture noisily display "`rng'"
capture noisily markin in 2/6
capture noisily display "`in'"
capture noisily markin if x > 100
capture noisily display "`in'"
gen byte mk = 1
capture noisily markout mk x
capture noisily count if mk
gen byte mk2 = 1
capture noisily markout mk2 x, sysmissok
capture noisily count if mk2
gen mkf = 1
capture noisily markout mkf x
display _rc
capture noisily mark m1
capture noisily mark m4 [pweight=w]
capture noisily mark m5 [aweight=w] if x > 2, zeroweight
gen t = _n
tsset t
capture program drop odta_myp
program define odta_myp
    syntax [varlist(fv ts)] [if] [in]
    marksample touse
    count if `touse'
end
capture noisily odta_myp L.x
capture noisily odta_myp x
capture program drop odta_myp

* postfile
set seed 12345
capture program drop odta_lnsim
program odta_lnsim
    version 14
    tempname sim
    postfile `sim' mean var using odta_results, replace
    quietly {
        forvalues i = 1/20 {
            drop _all
            set obs 100
            generate z = exp(rnormal())
            summarize z
            post `sim' (r(mean)) (r(Var))
        }
    }
    postclose `sim'
end
capture noisily odta_lnsim
capture noisily use odta_results, clear
capture noisily describe
tempname m
capture noisily postfile `m' a b str5 c double(d e) using odta_r2, replace
capture noisily post `m' (1) (2) ("abcdefgh") (3.5) (.)
capture noisily post `m' (1) (2)
display _rc
capture noisily post `m' (1) (2) ("a") (3.5) (.)
display _rc
capture noisily postutil dir
capture noisily postclose `m'
display _rc
capture noisily use odta_r2, clear
capture noisily list
capture noisily describe
tempname m2
capture noisily postfile `m2' a b using odta_r3, replace
capture noisily post `m2' (1) (2) (3)
display _rc
capture noisily postutil dir
tempname m3
capture noisily postfile `m3' a str5 b using odta_r4, replace
capture noisily post `m3' (1) (2)
display _rc
capture noisily postutil dir
tempname m4
capture noisily postfile `m4' a using odta_r5, replace
capture noisily post `m4' (nosuchvar)
display _rc
capture noisily postutil dir
tempname m5
capture noisily postfile `m5' a using odta_r6, replace
capture noisily postfile `m5' b using odta_r7, replace
display _rc
capture noisily postfile `m5' z using odta_r6
display _rc
capture noisily postutil dir
capture noisily postclose `m5'
capture noisily postfile `m5' z using odta_r6
display _rc
capture noisily postfile `m5' a a using odta_r8, replace
display _rc
capture noisily postfile `m5' _all using odta_r8, replace
display _rc
capture noisily postfile `m5' using odta_r8, replace
display _rc
capture noisily postutil clear
capture noisily postutil dir

* matlist
matrix A = (1,2 \ 3,4 \ 5,6)
capture noisily matlist A
capture noisily matlist A, border(rows) rowtitle(rows) left(4)
capture noisily matlist 2*A, border(all) lines(none) format(%6.1f) names(rows) twidth(8) left(4) title(Guess what, a title)
matrix Htest = (12.3, 2, .0004464 \ 2.17, 1, .3533287 \ 8.81, 3, .0402262 \ 20.05, 6, .0010676)
matrix rownames Htest = trunk length weight overall
matrix colnames Htest = chi2 df p
capture noisily matlist Htest
capture noisily matlist Htest, rowtitle(Variables) title(Test results) cspec(o4& %12s | %8.0g & %5.0f & %8.4f o2&) rspec(&-&&--)
matrix S = (1,.5 \ .5,1)
matrix rownames S = a_1 b_2
matrix colnames S = a_1 b_2
capture noisily matlist S
capture noisily matlist S, nohalf underscore lines(cells)
capture noisily matlist Htest, lines(rctotal)
capture noisily matlist Htest, nonames
capture noisily matlist Htest, noblank
capture noisily matlist Htest, aligncolnames(center)
capture noisily matlist Htest, aligncolnames(right) format(%9.2f)
capture noisily matlist Htest, border(top bottom) lines(oneline)
capture noisily matlist Htest, rowtitle(R) coltitle(C)
matrix N = J(2,2,.)
capture noisily matlist N
capture noisily matlist Zzz
capture noisily matlist Htest, bogus

* tabdisp
clear
set obs 40
gen rep78 = 1+mod(_n,5)
gen foreign = mod(_n,3)==0
label define fl 0 "Domestic" 1 "Foreign"
label values foreign fl
label define rl 1 "Poor" 2 "Fair" 3 "Average" 4 "Good" 5 "Excellent"
label values rep78 rl
label var rep78 "Repair Record 1978"
label var foreign "Car type"
gen mpg = 15 + _n/3
collapse (mean) mpg, by(foreign rep78)
capture noisily list
capture noisily tabdisp foreign rep78, cell(mpg) format(%9.2f) center
capture noisily table foreign rep78, c(mean mpg) format(%9.2f) center
capture noisily tabdisp foreign rep78, cell(mpg) format(%9.2f)
capture noisily table foreign rep78, c(mean mpg) format(%9.2f)
clear
input a b c
0 1 15
0 2 26
0 3 11
1 1 14
1 2 12
1 3 7
end
capture noisily tabdisp a b, cell(c)
gen str8 make = "car" + string(_n)
gen mpg = c*2
capture noisily tabdisp make, cell(mpg c)
capture noisily tabdisp a, cell(c mpg) totals
capture noisily tabdisp a b, cell(c) concise
capture noisily tabdisp a b, cell(c) missing
capture noisily tabdisp a b, cell(c) dotz
capture noisily tabdisp a b, cell(c mpg) cellwidth(8) csepwidth(3) scsepwidth(2) stubwidth(10)
capture noisily tabdisp a b, cell(c) by(make)
capture noisily tabdisp a b, cell(c) format(%5.0f %5.1f)
capture noisily tabdisp a b
capture noisily tabdisp a b, cell(c mpg c mpg c mpg c mpg)
capture noisily tabdisp a b if a == 5, cell(c)
clear
input str6 agecat str6 sex str12 party str12 reaction
"Old" "Female" "Democrat" "Favor"
"Old" "Male" "Democrat" "Favor"
"Old" "Female" "Republican" "Indifferent"
"Old" "Male" "Republican" "Disfavor"
"Young" "Female" "Democrat" "Strongly"
"Young" "Male" "Democrat" "Indiff"
"Young" "Female" "Republican" "Disfavor"
"Young" "Male" "Republican" "Disfavor"
end
label var party "Party Affiliation and Sex"
label var agecat "Age category"
capture noisily tabdisp agecat sex party, c(reaction) center
capture noisily tabdisp agecat sex party, c(reaction)
capture noisily tabdisp agecat sex party, c(reaction) by(party)
capture noisily tabdisp agecat sex, c(reaction) left

* tabdisp: widths, titles and separations (probes)
clear
input g h v
1 1 5
1 2 123456
2 1 1.5
2 2 -3
end
gen double w = v*1000
label var v "Value of v"
gen str2 bb = "ab"
capture noisily tabdisp g h, cell(v)
capture noisily tabdisp g, cell(v)
capture noisily tabdisp g, cell(w)
capture noisily tabdisp g, cell(v w)
capture noisily tabdisp g h, cell(v w)
capture noisily tabdisp g h, cell(v) format(%7.1f)
capture noisily tabdisp g h, cell(v) format(%12.2f)
capture noisily tabdisp g bb, cell(v)
capture noisily tabdisp g h, cell(v) center
capture noisily tabdisp g h, cell(v) csepwidth(5)
capture noisily tabdisp g h, cell(v) stubwidth(12)
capture noisily tabdisp g h, cell(v) cellwidth(12) center
clear
input str5 r str5 c str5 sc str9 val
"x" "p" "Green" "a"
"x" "q" "Green" "b"
"y" "p" "Red" "c"
"y" "q" "Red" "d"
end
capture noisily tabdisp r c sc, cell(val)
capture noisily tabdisp r c sc, cell(val) scsepwidth(5) csepwidth(1)
capture noisily tabdisp r c sc, cell(val) center
capture noisily tabdisp r c, cell(val) by(sc)

* program, sortpreserve
clear
set obs 6
gen x = 7 - _n
gen y = _n
sort y
capture program drop odta_mysort
program define odta_mysort, sortpreserve
    sort x
    display "inside: `_sortindex'"
    list y, noobs clean
end
capture noisily odta_mysort
capture noisily list x y, noobs clean
capture noisily display "`: sortedby'"
capture program drop odta_baddrop
program define odta_baddrop, sortpreserve
    drop `_sortindex'
end
capture noisily odta_baddrop
display _rc
capture program drop odta_errp
program define odta_errp, sortpreserve
    sort x
    error 198
end
capture noisily odta_errp
display _rc
capture noisily list x y, noobs clean
capture program drop odta_srt
program define odta_srt, sortpreserve
    gen z = _n
    drop in 1
    expand 2
end
capture noisily odta_srt
capture noisily list, noobs clean

* plugin
capture noisily program foo, plugin using(odta_nofile.plugin)
capture noisily plugin call foo
capture noisily plugin call

* pause, sleep
capture noisily pause off
capture noisily pause here we are
capture noisily pause
capture noisily sleep 50
capture noisily sleep
capture noisily sleep abc
capture noisily pause bogus

* signestimationsample, checkestimationsample
clear
set obs 30
gen x1 = _n
gen x2 = mod(_n*7,11)
gen y = 2*x1 - x2 + mod(_n,3)
capture noisily regress y x1 x2
capture noisily signestimationsample y x1 x2
capture noisily display e(datasignaturevars)
capture noisily display e(datasignature)
capture noisily checkestimationsample
replace x1 = x1 + 1 in 3
capture noisily checkestimationsample
replace x1 = x1 - 1 in 3
capture noisily checkestimationsample
capture noisily signestimationsample
capture noisily signestimationsample junk
capture noisily regress y x1 if x2 > 3
capture noisily signestimationsample y x1
capture noisily display e(datasignature)
capture noisily summarize y
capture noisily signestimationsample y

* makecns, matcproc
quietly regress y x1 x2
matrix C = (1, -1, 0, 0)
capture noisily makecns C, r
capture noisily return list
capture noisily matrix list e(Cns)
capture noisily matcproc T a Cx
capture noisily matrix list T
capture noisily matrix list a
capture noisily matrix list Cx
capture noisily makecns
display _rc
matrix C2 = (2, 0, -3, 5 \ 0, 1, 1, 0)
capture noisily makecns C2, r
capture noisily return list
capture noisily makecns Zzz
display _rc
capture noisily makecns Zzz
matrix Cbad = (1, -1)
capture noisily makecns Cbad
capture noisily matcproc T a

* viewsource
capture noisily viewsource nothing.ado
capture noisily viewsource
capture erase odta_results.dta
capture erase odta_r2.dta
capture erase odta_r3.dta
capture erase odta_r4.dta
capture erase odta_r5.dta
capture erase odta_r6.dta
capture erase odta_r7.dta
capture erase odta_r8.dta

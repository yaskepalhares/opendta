* Probe (survey): functions and commands that research do-files use. rc 133 = unknown function, 199 = unknown command.
version 14
set more off
clear
display "=== functions"
capture noisily display ym(2020,1)
display "  rc=" _rc
capture noisily display yq(2020,1)
display "  rc=" _rc
capture noisily display yw(2020,1)
display "  rc=" _rc
capture noisily display yh(2020,1)
display "  rc=" _rc
capture noisily display mofd(21915)
display "  rc=" _rc
capture noisily display dofm(720)
display "  rc=" _rc
capture noisily display qofd(21915)
display "  rc=" _rc
capture noisily display yofd(21915)
display "  rc=" _rc
capture noisily display wofd(21915)
display "  rc=" _rc
capture noisily display week(21915)
display "  rc=" _rc
capture noisily display halfyear(21915)
display "  rc=" _rc
capture noisily display mdyhms(1,1,2020,1,1,1)
display "  rc=" _rc
capture noisily display clock("01jan2020 10:00","DMYhm")
display "  rc=" _rc
capture noisily display daily("01jan2020","DMY")
display "  rc=" _rc
capture noisily display monthly("2020m1","YM")
display "  rc=" _rc
capture noisily display hours(1)
display "  rc=" _rc
capture noisily display dhms(1,1,1,1)
display "  rc=" _rc
capture noisily display cofd(21915)
display "  rc=" _rc
capture noisily display dofc(1893456000000)
display "  rc=" _rc
capture noisily display Cofc(1893456000000)
display "  rc=" _rc
capture noisily display msofhours(1)
display "  rc=" _rc
capture noisily display recode(5,1,3,5)
display "  rc=" _rc
capture noisily display irecode(5,1,3,5)
display "  rc=" _rc
capture noisily display autocode(5,4,0,10)
display "  rc=" _rc
capture noisily display fileexists("nada.xyz")
display "  rc=" _rc
capture noisily display binomialp(10,5,.5)
display "  rc=" _rc
capture noisily display poissonp(2,2)
display "  rc=" _rc
capture noisily display double(5)
display "  rc=" _rc
capture noisily display byte(5)
display "  rc=" _rc
capture noisily display long(5)
display "  rc=" _rc
capture noisily display reldif(1,1.1)
display "  rc=" _rc
capture noisily display trigamma(1)
display "  rc=" _rc
capture noisily display strtoname("a b")
display "  rc=" _rc
capture noisily display plural(2,"cat")
display "  rc=" _rc
capture noisily display soundex("robert")
display "  rc=" _rc
capture noisily display strrpos("abcb","b")
display "  rc=" _rc
capture noisily display ustrpos("abc","b")
display "  rc=" _rc
capture noisily display usubstr("abc",1,1)
display "  rc=" _rc
capture noisily display ustrleft("abc",1)
display "  rc=" _rc
capture noisily display ustrtrim(" a ")
display "  rc=" _rc
capture noisily display ustrregexm("abc","b")
display "  rc=" _rc
capture noisily display maxbyte()
display "  rc=" _rc
capture noisily display maxint()
display "  rc=" _rc
capture noisily display logistic(0)
display "  rc=" _rc
capture noisily display tden(5,0)
display "  rc=" _rc
capture noisily display invibeta(1,1,.5)
display "  rc=" _rc
capture noisily display nchi2(1,1,1)
display "  rc=" _rc
capture noisily display invnchi2(1,1,.5)
display "  rc=" _rc
capture noisily display hypergeometricp(10,5,3,2)
display "  rc=" _rc
capture noisily display nbinomial(5,2,.5)
display "  rc=" _rc
capture noisily display indexnot("a","b")
display "  rc=" _rc
capture noisily display strdup("a",2)
display "  rc=" _rc
display "=== commands called without arguments (rc 199 = not a command)"
capture noisily mkmat
display "  rc=" _rc
capture noisily svmat
display "  rc=" _rc
capture noisily move
display "  rc=" _rc
capture noisily aorder
display "  rc=" _rc
capture noisily numlabel
display "  rc=" _rc
capture noisily numlist
display "  rc=" _rc
capture noisily bsample
display "  rc=" _rc
capture noisily iqreg
display "  rc=" _rc
capture noisily sqreg
display "  rc=" _rc
capture noisily bsqreg
display "  rc=" _rc
capture noisily treatreg
display "  rc=" _rc
capture noisily forecast
display "  rc=" _rc
capture noisily fvset
display "  rc=" _rc
capture noisily glcurve
display "  rc=" _rc
capture noisily hettest
display "  rc=" _rc
capture noisily ovtest
display "  rc=" _rc
capture noisily vif
display "  rc=" _rc
capture noisily lstat
display "  rc=" _rc
capture noisily pk
display "  rc=" _rc
capture noisily pkexamine
display "  rc=" _rc
capture noisily qladder
display "  rc=" _rc
capture noisily qqplot
display "  rc=" _rc
capture noisily robvar
display "  rc=" _rc
capture noisily sdtesti
display "  rc=" _rc
capture noisily sfrancia
display "  rc=" _rc
capture noisily svydes
display "  rc=" _rc
capture noisily svylc
display "  rc=" _rc
capture noisily rcof
display "  rc=" _rc
capture noisily renpfix
display "  rc=" _rc
capture noisily set_defaults
display "  rc=" _rc
capture noisily hexdump
display "  rc=" _rc
capture noisily prop
display "  rc=" _rc
capture noisily query memory
display "  rc=" _rc
capture noisily memory
display "  rc=" _rc

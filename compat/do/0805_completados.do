* Pieces completed for 1.0: anova (sequential, nested terms, error terms, repeated measures),
* xtgee and the pa estimator, xtreg be and mle, two-way tabulate summarize(), misstable
* tree/nested, ci variances, two-sample prtest, reshape advanced syntax and reshape error,
* old merge syntax, sort in and list sepby()
clear
input person drug score
1 1 30
1 2 28
1 3 16
1 4 34
2 1 14
2 2 18
2 3 10
2 4 22
3 1 24
3 2 20
3 3 18
3 4 30
4 1 38
4 2 34
4 3 20
4 4 44
5 1 26
5 2 28
5 3 14
5 4 30
end
anova score person drug, repeated(drug)
ereturn list
anova score drug person, sequential
gen g = person > 3
anova score g / person|g drug g#drug, repeated(drug)
anova score g / person|g /
tabulate g drug, summarize(score)
tabulate g drug, summarize(score) means

clear
set seed 805
set obs 200
gen id = ceil(_n/5)
gen t = mod(_n-1, 5) + 1
xtset id t
gen x = round(rnormal(), .01)
gen u = (mod(id, 7) - 3)/3
gen y = 1 + .5*x + u + round(rnormal(), .01)
gen yb = runiform() < invlogit(-.3 + .8*x + u)
gen yc = floor(exp(.2 + .3*x + .5*u) * runiform()*2)
xtgee y x
xtgee y x, corr(ar 1)
xtgee y x, corr(independent) vce(robust)
xtreg y x, pa
xtlogit yb x, pa
xtpoisson yc x, pa vce(robust)
ereturn list
xtreg y x, be
xtreg y x, mle
quietly xtreg y x, fe
predict uhat, u
predict ehat, e
summarize uhat ehat

generate g2 = mod(_n, 2)
prtest yb, by(g2)
prtesti 30 .4 45 .6
ci variances x y
ci variances x, sd
ci variances x, bonett
ci proportions yb, wilson
cii variances 40 2.25

replace x = . in 1/20
replace y = . in 11/40
replace yc = . in 191/200
misstable tree x y yc
misstable tree x y yc, frequency
misstable nested x y yc

clear
input id inc80 inc81 ue80 ue81
1 5000 5500 0 1
2 2000 2200 1 0
3 3000 2000 0 0
end
reshape i id
reshape j year
reshape xij inc ue
reshape
reshape long
reshape wide
reshape clear
clear
input id year inc
1 80 5
1 80 6
2 80 7
end
capture noisily reshape wide inc, i(id) j(year)
reshape error

clear
input id b
1 10
2 20
4 40
end
save odta_u, replace
clear
input id a
3 3
1 1
2 2
end
capture noisily merge id using odta_u
sort id
merge id using odta_u
list, sepby(_merge)
erase odta_u.dta

clear
input x
5
4
3
2
1
end
sort x in 2/4
list

* Epidemiology: cs, cc and ir with pool, istandard, estandard, standard(), rd, ird,
* binomial(), bd and tarone; mhodds with by(); roctab with bamber, binomial and lorenz;
* roccomp with binormal; rocreg (nonparametric and parametric, bootstrap) (layouts, r(), e())
clear
input agecat smokes deaths pyears
1 1 32 52407
2 1 104 43248
3 1 206 28612
4 1 186 12663
5 1 102 5317
1 0 2 18790
2 0 12 10673
3 0 28 5710
4 0 28 2585
5 0 31 1462
end
ir deaths smokes pyears, by(agecat) level(90) istandard pool
return list
ir deaths smokes pyears, by(agecat) level(90) istandard ird
return list
gen conswgt = 1
ir deaths smokes pyears, by(agecat) standard(conswgt)
ir deaths smokes pyears, by(agecat) estandard
ir deaths smokes pyears, by(agecat) pool
return list

clear
input case exp age pop
1 1 0 10
1 0 0 20
0 1 0 90
0 0 0 80
1 1 1 25
1 0 1 15
0 1 1 75
0 0 1 85
1 1 2 30
1 0 2 12
0 1 2 50
0 0 2 70
end
cs case exp [fw=pop], by(age) istandard
return list
cs case exp [fw=pop], by(age) estandard pool
return list
gen wgt = 1
cs case exp [fw=pop], by(age) standard(wgt) rd
return list
cc case exp [fw=pop], by(age) bd
return list
cc case exp [fw=pop], by(age) tb tarone
return list
cc case exp [fw=pop], by(age) istandard pool
return list
cc case exp [fw=pop], by(age) estandard
return list
cc case exp [fw=pop], by(age) pool woolf
return list
mhodds case exp [fw=pop], by(age)
return list
capture noisily cs case exp [fw=pop], pool
capture noisily cs case exp [fw=pop], by(age) rd
capture noisily csi 10 20 90 80, istandard

* binomial frequency records
clear
input cases n exp age
10 100 1 0
20 100 0 0
25 100 1 1
15 100 0 1
end
cs cases exp, binomial(n)
return list
cs cases exp, by(age) binomial(n)
cc cases exp, by(age) binomial(n)

* tabodds/mhodds data with a three-level exposure
clear
input case alc age pop
1 1 1 2
0 1 1 60
1 2 1 4
0 2 1 40
1 3 1 8
0 3 1 20
1 1 2 10
0 1 2 50
1 2 2 15
0 2 2 30
1 3 2 25
0 3 2 15
end
mhodds case alc [fw=pop], by(age)
return list
mhodds case alc [fw=pop], compare(3,1) by(age)

* ROC: Hanley and McNeil (1982) ratings
clear
input disease rating pop
0 1 33
0 2 6
0 3 6
0 4 11
0 5 2
1 1 3
1 2 2
1 3 2
1 4 11
1 5 33
end
roctab disease rating [fw=pop], bamber
return list
roctab disease rating [fw=pop], binomial
roctab disease rating [fw=pop], hanley binomial
return list
roctab disease rating [fw=pop], lorenz
return list
expand pop
* male before roccomp: Stata's roccomp, by() may reorder the data (compat 1601)
gen male = mod(_n, 2)
set seed 1601
gen r2 = rating + (disease==1)*(runiform()<.3)
gen grp = runiform() < .5
roccomp disease rating r2, binormal
return list
roccomp disease rating, by(grp) binormal

* rocreg: the bootstrap replicates depend on Stata's resampling order (VERIFY)
rocreg disease rating r2, nobootstrap
ereturn list
rocreg disease rating, tiecorrected nobootstrap
rocreg disease r2, ctrlcov(male) auc roc(.2 .4) invroc(.8) pauc(.3) nobootstrap
rocreg disease r2, ctrlcov(male) ctrlmodel(linear) pvc(normal) nobootstrap
rocreg disease r2, probit nobootstrap
ereturn list
rocreg disease r2, probit roccov(male) fprpts(5) nobootstrap
rocreg disease rating, breps(50) bseed(1601) nodots
ereturn list

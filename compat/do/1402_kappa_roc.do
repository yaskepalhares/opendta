* kap and kappa (interrater agreement); roctab, roccomp, rocfit, lroc and lsens
clear
input a b n
1 1 20
1 2 5
2 1 10
2 2 15
end
kap a b [fw=n]
return list
kap a b [fw=n], tab
clear
input r1 r2 r3
1 1 1
1 2 1
2 2 2
3 3 2
1 1 1
2 2 3
3 3 3
1 2 2
end
kap r1 r2
kap r1 r2, wgt(w)
kap r1 r2, wgt(w2)
kap r1 r2 r3
return list
clear
input c1 c2 c3
3 0 0
2 1 0
0 3 0
0 1 2
1 1 1
0 0 3
end
kappa c1 c2 c3
return list
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
roctab disease rating [fw=pop]
return list
roctab disease rating [fw=pop], detail
roctab disease rating [fw=pop], table
roctab disease rating [fw=pop], hanley
rocfit disease rating [fw=pop]
ereturn list
expand pop
set seed 1402
gen r2 = rating + (disease == 1)*(runiform() < .3)
gen grp = mod(_n, 2)
roccomp disease rating r2
return list
roccomp disease rating, by(grp)
clear
set seed 1402
set obs 200
gen x = invnormal(runiform())
gen z = invnormal(runiform())
gen y = runiform() < invlogit(-0.5 + 1.2*x + 0.5*z)
logit y x z, nolog
lroc, nograph
return list
lsens, nograph
probit y x, nolog
lroc, nograph

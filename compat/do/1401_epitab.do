* epitab: cs, cc, ir, mcc and their immediate forms, tabodds, mhodds, symmetry,
* dstdize and istdize (layouts, intervals, tests and r() results)
csi 7 12 9 2
return list
csi 7 12 9 2, or
csi 7 12 9 2, tb exact
cci 4 386 4 1250
return list
cci 4 386 4 1250, exact
cci 4 386 4 1250, cornfield
cci 4 386 4 1250, woolf
cci 20 10 30 40, tb
iri 41 15 28010 19017
return list
iri 41 15 28010 19017, exact
iri 10 30 28010 19017
mcci 8 8 3 8
return list
mcci 8 8 3 8, tb
clear
input case exp sex pop
1 1 0 10
1 0 0 20
0 1 0 90
0 0 0 80
1 1 1 25
1 0 1 15
0 1 1 75
0 0 1 85
end
cs case exp [fw=pop]
cs case exp [fw=pop], by(sex)
return list
cc case exp [fw=pop]
cc case exp [fw=pop], by(sex)
return list
cc case exp [fw=pop], by(sex) woolf
clear
input exp cases time age
1 10 2000 1
0 8 3000 1
1 25 1500 2
0 12 2500 2
end
ir cases exp time
ir cases exp time, by(age)
return list
clear
input casee ctrle pop
1 1 8
1 0 8
0 1 3
0 0 8
end
mcc casee ctrle [fw=pop]
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
tabodds case alc [fw=pop]
return list
tabodds case alc [fw=pop], or
tabodds case alc [fw=pop], or base(2)
tabodds case alc [fw=pop], or adjust(age)
mhodds case alc [fw=pop]
mhodds case alc age [fw=pop]
return list
mhodds case alc age [fw=pop], compare(3,1)
gen e = alc == 3
mhodds case e age [fw=pop]
clear
input r c n
1 1 20
1 2 5
1 3 2
2 1 3
2 2 30
2 3 6
3 1 1
3 2 4
3 3 25
end
symmetry r c [fw=n]
return list
symmetry r c [fw=n], contrib
symmetry r c [fw=n], notable exact mh
clear
input agecat stdpop
1 300000
2 200000
3 80000
end
gen stddeaths = cond(agecat == 1, 50, cond(agecat == 2, 140, 400))
save odta_std, replace
clear
input str6 state agecat pop deaths
Maine 1 100000 20
Maine 2 80000 60
Maine 3 40000 200
Iowa 1 150000 25
Iowa 2 90000 70
Iowa 3 30000 160
end
dstdize deaths pop agecat, by(state)
return list
matrix list r(adj)
dstdize deaths pop agecat, by(state) base(Maine)
istdize deaths pop agecat using odta_std, popvars(stddeaths stdpop) by(state)
return list
erase odta_std.dta

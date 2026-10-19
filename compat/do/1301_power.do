* power: sample size, power and effect size of the essential methods
power onemean 100 102, sd(5)
return list
power onemean 100 102, sd(5) knownsd
power onemean 100 102, sd(5) n(40)
capture noisily power onemean 100, sd(5) n(40)
power onemean 100, sd(5) n(40) power(.8)
power onemean 100 102, sd(5) onesided
power twomeans 12 14, sd(4)
return list
power twomeans 12 14, sd(4) power(.9)
power twomeans 12 14, sd1(4) sd2(6) nratio(2)
power twomeans 12 14, sd(4) n(100)
power pairedmeans 10 12, sddiff(5)
power pairedmeans 10 12, corr(.5) sd(4)
power oneproportion .3 .4
power oneproportion .3 .4, test(wald)
capture noisily power oneproportion .3 .4, test(binomial)
power oneproportion .3 .4, test(binomial) n(150)
power oneproportion .3 .4, n(150)
power twoproportions .3 .45
return list
power twoproportions .3 .45, n(200)
power twoproportions .3 .45, nratio(2)
capture noisily power twoproportions .3, n(300)
power twoproportions .3, n(300) power(.8)
power onecorrelation 0 .3
power twocorrelations .3 .5
power onemean 100 (102 104 106), sd(5)
power twomeans 12 (13 14), sd(4) n(50 100)

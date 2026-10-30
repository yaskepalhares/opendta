* power: twoway, repeated, cmh, mcc, trend, cox, exponential and logrank
* Everything the manual does not settle is under capture noisily (VERIFY in src/opendta/commands/power3.py)

* twoway: manual example (N = 132; with cell weights N = 153)
power twoway 134 143 91 \ 106 173 145, varerror(1417) power(.8)
return list
power twoway 134 143 91 \ 106 173 145, varerror(1417) power(.8) cellweights(2 2 2 \ 1 1 1)
power twoway 134 143 91 \ 106 173 145, varerror(1417) n(132)
power twoway 134 143 91 \ 106 173 145, varerror(1417) factor(column) power(.8)
power twoway 134 143 91 \ 106 173 145, varerror(1417) factor(rowcol) power(.8)
capture noisily power twoway, nrows(2) ncols(3) n(132) power(.8) varerror(1417)
capture noisily power twoway 134 143 91 \ 106 173 145, varerror(1417) n(100(20)140)
capture noisily power twoway 134 143 91 \ 106 173 145, varerror(1417) power(.8) showcellsizes

* repeated
capture noisily power repeated, ngroups(2) nrepeated(3) corr(.5) varbetween(1) n(40)
return list
capture noisily power repeated 10 12 14 \ 11 13 12, corr(.6) varerror(77) power(.9) factor(within)
capture noisily power repeated 10 12 14 \ 11 13 12, corr(.6) varerror(77) power(.9) factor(bwithin)
matrix S = (4, 2, 1 \ 2, 4, 2 \ 1, 2, 4)
capture noisily power repeated, ngroups(2) covmatrix(S) factor(within) vareffect(.4) n(60)

* cmh
power cmh .2 .1 .3, oratio(2.5)
return list
power cmh .2 .1 .3, oratio(2.5) n(240)
capture noisily power cmh .2 .1 .3, n(240) power(.8)
capture noisily power cmh .2 .1 .3, oratio(2.5) strweights(4 1 4) grratios(.47 .57 .51)
capture noisily power cmh .2 .1 .3, oratio(2.5) continuity
capture noisily power cmh .2 .1 .3, oratio(2 2.5 3)
capture noisily power cmh .2 .1 .3, oratio(2.5) onesided

* mcc
power mcc .3, oratio(2)
return list
power mcc .3, oratio(2) n(300)
capture noisily power mcc .3, oratio(2) m(2)
capture noisily power mcc .3, oratio(2) m(2) compare
capture noisily power mcc .3, n(300) power(.8)
capture noisily power mcc .3, oratio(2) corr(.2)

* trend
power trend .1 .2 .3, power(.8)
return list
power trend .1 .2 .3, n(300)
capture noisily power trend .1 .2 .3, n(300) exposure(1 2 4)
capture noisily power trend .1 .3 .2, n(300)
capture noisily power trend .1 .2 .3, power(.8) continuity
capture noisily power trend .1 .2 .3, power(.8) grweights(1 1 2)

* cox
power cox
return list
capture noisily power cox, hratio(.5) eventprob(.5) wdprob(.2)
capture noisily power cox, n(200) power(.8)
capture noisily power cox, hratio(.5) r2(.5)
capture noisily power cox, hratio(.5) sd(.5) n(66)

* exponential
power exponential .1 .05, power(.8)
return list
power exponential .1 .05, n(300)
capture noisily power exponential .5 .25, time(2) power(.8)
capture noisily power exponential .1 .05, power(.8) fperiod(5) aperiod(2)
capture noisily power exponential .1 .05, power(.8) fperiod(5) aperiod(2) lossprob(.1)
capture noisily power exponential .1 .05, power(.8) unconditional
capture noisily power exponential, hratio(.6) loghazard power(.8)

* logrank
power logrank
return list
capture noisily power logrank, hratio(.76) schoenfeld power(.9)
capture noisily power logrank .5, hratio(.7) wdprob(.1)
capture noisily power logrank .5 .3, n(200)
capture noisily power logrank .5, n(200) power(.8)
capture noisily power logrank, simpson(.8 .7 .6) hratio(.6)

* lists of values
capture noisily power cox, hratio(.5 .6) power(.8 .9)
capture noisily power logrank, hratio(.5 .6 .7) power(.8 .9)

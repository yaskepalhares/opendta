* Phase 23: other suites (sem with mlmv, adf, group() and the estat subcommands; gsem; irt; bayesmh
* and the bayes prefix; the other power methods; mi impute mvn, monotone, truncreg, intreg and nbreg;
* svy with more estimators). Every command runs under capture noisily so that a difference in one
* command does not stop the rest of the case; the layouts marked VERIFY in the code are settled by the
* reference log. The data are built from rounded invnormal(runiform()) draws and from comparisons of
* runiform() draws (rnormal(), rpoisson() and rbinomial() do not reproduce Stata's draws); counts are
* geometric draws by inversion, so that they are overdispersed. Each data set is summarized right after
* it is built, so that the log shows at once whether both programs made the same data. The draws of
* bayesmh and mi impute come from each program's own generator and are not expected to match.

* ---- sem: two factors with six indicators, missing values and two groups
clear
set seed 1909
set obs 500
gen double f = round(invnormal(runiform()), .01)
gen double g = round(.5*f + .8*invnormal(runiform()), .01)
gen double x1 = round(f + .6*invnormal(runiform()), .01)
gen double x2 = round(.8*f + .6*invnormal(runiform()), .01)
gen double x3 = round(1.2*f + .8*invnormal(runiform()), .01)
gen double y1 = round(g + .6*invnormal(runiform()), .01)
gen double y2 = round(.6*g + .6*invnormal(runiform()), .01)
gen double y3 = round(1.4*g + .8*invnormal(runiform()), .01)
gen byte grp = 1 + (runiform() < .5)
replace x2 = round(x2 + .4*f, .01) if grp == 2
gen double x2m = x2
replace x2m = . if runiform() < .1
gen double y3m = y3
replace y3m = . if runiform() < .1
gen double z1 = round(invnormal(runiform()), .01)
gen double z2 = round(invnormal(runiform()), .01)
gen double m = round(.6*z1 + invnormal(runiform()), .01)
gen double w = round(.4*m + .3*z1 + invnormal(runiform()), .01)
gen double u1 = round(invnormal(runiform()), .01)
gen double u2 = round(.4*u1 + invnormal(runiform()), .01)
gen double v1 = round(.3*u2 + .6*z1 + u1, .01)
gen double v2 = round(.4*v1 + .6*z2 + u2, .01)
summarize
capture noisily sem (F -> x1 x2m x3) (G -> y1 y2 y3m), method(mlmv) nolog
capture noisily ereturn list
capture noisily sem (F -> x1 x2 x3) (G -> y1 y2 y3), method(adf) nolog
capture noisily sem (F -> x1 x2 x3) (G -> y1 y2 y3), nolog
capture noisily estat eqgof
capture noisily estat framework
capture noisily estat residuals
capture noisily estat residuals, standardized
capture noisily estat summarize
capture noisily estat mindices
capture noisily estat scoretests
capture noisily estat eqtest
capture noisily sem (F -> x1 x2 x3 y1), group(grp) nolog
capture noisily estat ggof
capture noisily estat ginvariant
capture noisily sem (F -> x1 x2 x3 y1), group(grp) ginvariant(none) nolog
capture noisily sem (y1 <- z1@a z2@b) (y2 <- z1@a z2@b), nolog
capture noisily sem (m <- z1) (w <- m z1), nolog
capture noisily estat teffects
capture noisily estat stdize: lincom [w]m + [w]z1
capture noisily sem (v1 <- v2 z1) (v2 <- v1 z2), nolog
capture noisily estat stable

* ---- gsem: binary, ordinal, count and gamma responses, random intercepts and measurement models
clear
set seed 190901
set obs 60
gen int id = _n
gen double re = round(.8*invnormal(runiform()), .01)
expand 10
sort id, stable
gen double x = round(invnormal(runiform()), .01)
gen double l = round(invnormal(runiform()), .01)
gen byte yb = runiform() < invlogit(-.4 + .8*x + re)
gen double yo0 = .8*x + re + ln(runiform()/(1 - runiform()))
gen byte yo = (yo0 > -1) + (yo0 > .4) + (yo0 > 1.6)
gen double mu = exp(.4 + .3*x)
gen int yc = floor(ln(runiform())/ln(mu/(1 + mu)))
gen double yg = round(-exp(.5 + .3*x)*ln(runiform()*runiform()), .01) + .01
gen byte b1 = runiform() < invlogit(.2 + 1*l)
gen byte b2 = runiform() < invlogit(-.3 + 1.4*l)
gen byte b3 = runiform() < invlogit(.5 + .8*l)
gen byte b4 = runiform() < invlogit(-.1 + 1.2*l)
gen double c1 = round(l + .6*invnormal(runiform()), .01)
gen double c2 = round(.8*l + .6*invnormal(runiform()), .01)
summarize
capture noisily gsem (yb <- x, logit), nolog
capture noisily ereturn list
capture noisily gsem (yb <- x, probit), nolog
capture noisily gsem (yc <- x, poisson), nolog
capture noisily gsem (yc <- x, nbreg), nolog
capture noisily gsem (yg <- x, gamma), nolog
capture noisily gsem (yo <- x, ologit), nolog
capture noisily gsem (yb <- x M1[id], logit), nolog
capture noisily gsem (yc <- x M1[id], poisson), nolog
capture noisily gsem (L -> b1 b2 b3 b4, logit), nolog
capture noisily gsem (L -> c1 c2 b1 b2), nolog
capture noisily gsem (L -> c1 c2) (L -> b1 b2 b3, logit), nolog
capture noisily predict double pb1, mu(b1)
capture noisily predict double lhat, latent(L)
capture noisily summarize pb1 lhat

* ---- irt: six binary items and four items with four ordered categories
clear
set seed 190902
set obs 800
gen double th = round(invnormal(runiform()), .01)
gen byte q1 = runiform() < invlogit(1*(th + 1))
gen byte q2 = runiform() < invlogit(1.4*(th + .4))
gen byte q3 = runiform() < invlogit(.8*(th - 0))
gen byte q4 = runiform() < invlogit(1.2*(th - .6))
gen byte q5 = runiform() < invlogit(1.6*(th - 1))
gen byte q6 = runiform() < invlogit(1.2*(th - .2))
forvalues j = 1/6 {
    * two draws in separate variables, so that no logical operator decides whether a draw is made
    gen double ug = runiform()
    gen double uh = runiform()
    gen byte g`j' = ug < .2 | uh < invlogit((.6 + .2*`j')*(th - 1.2 + .4*`j'))
    drop ug uh
}
forvalues j = 1/4 {
    gen double s`j' = (.8 + .2*`j')*th + ln(runiform()/(1 - runiform()))
    gen byte r`j' = (s`j' > -1.6 + .2*`j') + (s`j' > -.2 + .2*`j') + (s`j' > 1.2 + .2*`j')
}
summarize q* g* r*
capture noisily irt 1pl q1-q6, nolog
capture noisily irt 2pl q1-q6, nolog
capture noisily ereturn list
capture noisily estat report
capture noisily estat report, byparm
capture noisily estat report q2 q4, sort(b)
capture noisily irtgraph icc q1 q2
capture noisily irtgraph tif
capture noisily irt 3pl g1-g6, nolog
capture noisily irt grm r1-r4, nolog
capture noisily irt pcm r1-r4, nolog
capture noisily irt gpcm r1-r4, nolog
capture noisily irt rsm r1-r4, nolog
capture noisily irt nrm r1-r4, nolog
capture noisily irt hybrid (2pl q1-q4) (grm r1-r2), nolog

* ---- bayesmh and the bayes prefix
clear
set seed 190903
set obs 200
gen double x = round(invnormal(runiform()), .01)
gen double z = round(invnormal(runiform()), .01)
gen double y = round(1 + .5*x - .3*z + invnormal(runiform()), .01)
gen byte b = runiform() < invlogit(-.2 + .8*x)
summarize
capture noisily bayesmh y x z, likelihood(normal({var})) prior({y:}, normal(0,100)) prior({var}, igamma(.01,.01)) block({var}) rseed(17)
capture noisily ereturn list
capture noisily bayesstats summary {y:x} {var}
capture noisily bayesstats ess
capture noisily bayestest interval {y:x}, lower(0)
capture noisily estimates store full
capture noisily bayesmh y x, likelihood(normal({var})) prior({y:}, normal(0,100)) prior({var}, igamma(.01,.01)) block({var}) rseed(17)
capture noisily estimates store small
capture noisily bayesstats ic full small
capture noisily bayestest model full small
capture noisily bayes, rseed(17): logit b x
capture noisily bayes, rseed(17): regress y x z
capture noisily bayes, rseed(17) hpd: poisson b x

* ---- power: the other methods, tests and graph
capture noisily power onevariance 4 9
capture noisily return list
capture noisily power onevariance 4 9, n(30)
capture noisily power onevariance 2 3, sd
capture noisily power onevariance 4, n(30) power(.8)
capture noisily power twovariances 4 2.25
capture noisily power twovariances 4 2.25, n(250)
capture noisily power twovariances 4, n(250) power(.8)
capture noisily power twovariances 4 (1.5(.5)3), n(250)
capture noisily power pairedproportions .105 .004
capture noisily power pairedproportions .105 .004, n(100)
capture noisily power pairedproportions .53 .4293, corr(.8)
capture noisily power twoproportions .015 .001, test(lrchi2)
capture noisily power twoproportions .015 .001, test(lrchi2) n(1100)
capture noisily power twoproportions .1 .2, test(fisher)
capture noisily power twoproportions .2 .6, test(fisher) n(24)
capture noisily return list
capture noisily power oneway 1 2 3, varerror(4)
capture noisily power oneway 1 2 3, varerror(4) n(60)
capture noisily power twovariances 4 (1.5 2 3), n(250) graph

* ---- mi impute: mvn, monotone, truncreg, intreg and nbreg
clear
set seed 190904
set obs 300
gen double x = round(invnormal(runiform()), .01)
gen double y1 = round(1 + .5*x + invnormal(runiform()), .01)
gen double y2 = round(-.5 + .4*x + .4*y1 + invnormal(runiform()), .01)
gen double yt = round(1 + .5*x + invnormal(runiform()), .01)
replace yt = . if yt <= 0
gen double yi = round(1 + .5*x + invnormal(runiform()), .01)
gen double lo = yi
gen double hi = yi
replace lo = floor(yi) if runiform() < .3
replace hi = lo + 1 if lo < yi
replace yi = . if lo < hi
gen double mu = exp(.4 + .3*x)
gen int yc = floor(ln(runiform())/ln(mu/(1 + mu)))
gen double rm = runiform()
replace y1 = . if rm < .15
replace y2 = . if rm < .3
replace yc = . if runiform() < .2
replace yt = . if runiform() < .15 & yt < .
summarize
mi set mlong
mi register imputed y1 y2 yt yi yc
capture noisily mi impute mvn y1 y2 = x, add(5) rseed(11)
capture noisily mi impute monotone (regress) y1 (regress) y2 = x, replace rseed(11)
capture noisily mi impute truncreg yt = x, ll(0) replace rseed(11)
capture noisily mi impute intreg yi = x, ll(lo) ul(hi) replace rseed(11)
capture noisily mi impute nbreg yc = x, replace rseed(11)
capture noisily mi describe

* ---- svy: more estimators
clear
set seed 190905
set obs 600
gen int psu = ceil(_n/10)
gen byte str = 1 + (psu > 30)
gen double wt = round(1 + 3*runiform(), .01)
gen double x = round(invnormal(runiform()), .01)
gen double z = round(invnormal(runiform()), .01)
gen byte b = runiform() < 1 - exp(-exp(-.4 + .6*x))
gen double ys = round(1 + .5*x + invnormal(runiform()), .01)
gen double ylo = floor(ys)
gen double yhi = ylo + 1
gen double yt = ys if ys > 0
gen byte sel = runiform() < normal(.3 + .6*z + .3*x)
gen double yh = round(1 + .5*x + invnormal(runiform()), .01)
replace yh = . if !sel
gen double mu = exp(1.2 + .3*x)
gen int yc = floor(ln(runiform())/ln(mu/(1 + mu)))
replace yc = 0 if runiform() < invlogit(-1 + .6*z)
gen byte b2 = runiform() < normal(-.2 + .5*z + .3*x)
gen double fr = round(invlogit(.2 + .5*x + .5*invnormal(runiform())), .01)
gen byte bh = runiform() < normal(.2 + .5*x)
replace bh = . if !sel
gen byte oc = 1 + (ys > .5) + (ys > 1.5)
gen double zz = round(.6*x + invnormal(runiform()), .01)
gen int grp = ceil(_n/3)
summarize
svyset psu [pweight=wt], strata(str)
capture noisily svy: cloglog b x
capture noisily ereturn list
capture noisily svy: intreg ylo yhi x
capture noisily svy: truncreg yt x, ll(0)
capture noisily svy: heckman yh x, select(sel = x z)
capture noisily svy: zip yc x, inflate(z)
capture noisily svy: zinb yc x, inflate(z)
capture noisily svy: gnbreg yc x
capture noisily svy: hetprobit b x, het(z)
capture noisily svy: biprobit b b2 x z
capture noisily svy: fracreg probit fr x
capture noisily svy: heckprobit bh x, select(sel = x z)
capture noisily svy: scobit b x
capture noisily svy: slogit oc x z
capture noisily svy: mprobit oc x z
capture noisily svy: ivregress 2sls ys (x = zz)
capture noisily svy: clogit b x, group(grp)
capture noisily svy, subpop(if x > -1): cloglog b x

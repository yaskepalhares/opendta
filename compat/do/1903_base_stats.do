* Phase 17: base statistics (resampling, ml, model building, postestimation, tests, agreement,
* transformations and displays). Every command runs under capture noisily so that a difference
* in one command does not stop the rest of the case; the layouts marked VERIFY in the code are
* settled by the reference log.
* Stata's example datasets are not available in OpenDTA, so the case builds its own data.
* Everything is rounded to integers (or cents), so that the last digit of the draws does not matter.
program define mkdata
    clear
    set seed 1903
    set obs 74
    gen int weight = round(3000 + 700*invnormal(runiform()))
    gen int length = round(weight/17 + 4*invnormal(runiform()))
    gen int mpg = max(12, round(44 - weight/110 + 2.5*invnormal(runiform())))
    gen int displacement = round(weight/12 + 30*invnormal(runiform()))
    gen int turn = round(length/5 + 2*invnormal(runiform()))
    gen byte foreign = runiform() < 1/(1 + exp((weight - 2800)/400))
    gen int price = round(exp(7.8 + .00025*weight - .01*mpg + .35*invnormal(runiform())))
    gen byte rep78 = 1 + floor(5*runiform())
    replace rep78 = . if runiform() < .07
    label variable mpg "Mileage (mpg)"
end
set more off
mkdata
gen hiprice = price > 6000
gen hiwt = weight > 3000
gen double lnprice0 = ln(price)

* --- means, tests and correlations ---
capture noisily ameans price mpg weight
capture noisily ameans price, level(90)
capture noisily bitest foreign == .5
capture noisily bitesti 74 22 .5
capture noisily bitesti 20 12 .3
capture noisily ksmirnov mpg, by(foreign)
capture noisily ksmirnov mpg, by(foreign) exact
capture noisily ksmirnov price = normal((price-6165.257)/2949.496)
capture noisily median mpg, by(foreign)
capture noisily median mpg, by(foreign) exact
capture noisily median mpg, by(rep78)
capture noisily nptrend mpg, by(rep78)
capture noisily nptrend mpg, by(rep78) nodetail
capture noisily runtest price
capture noisily runtest price, threshold(5000) continuity
capture noisily tetrachoric hiprice hiwt
capture noisily tetrachoric hiprice hiwt foreign
capture noisily ktau mpg weight
capture noisily ktau mpg weight price
capture noisily ktau mpg weight price, stats(taub p) print(.05)
capture noisily pcorr mpg weight length foreign
capture noisily return list

* --- agreement and effect sizes ---
preserve
clear
input r1 r2 r3 r4
9 2 5 8
6 1 3 2
8 4 6 8
7 1 2 6
10 5 6 9
6 2 4 7
end
gen target = _n
reshape long r, i(target) j(judge)
capture noisily icc r target judge
capture noisily icc r target judge, absolute
capture noisily icc r target judge, mixed
capture noisily icc r target
capture noisily return list
restore
capture noisily loneway price rep78
capture noisily loneway mpg rep78, level(90)
capture noisily esize twosample price, by(foreign)
capture noisily esize twosample mpg, by(foreign) all
capture noisily esizei 20 5 2 25 6 2.5
capture noisily return list

* --- transformations ---
capture noisily lnskew0 lnp = price
capture noisily lnskew0 lnw = weight, level(90)
capture noisily ladder price
capture noisily ladder mpg
capture noisily ladder price, generate(ladp)
capture noisily boxcox mpg, nolog
capture noisily boxcox mpg weight, nolog
capture noisily ereturn list

* --- displays ---
capture noisily cumul price, generate(cprice)
capture noisily cumul mpg, generate(cmpg) equal
capture noisily summarize cprice cmpg
capture noisily stem mpg
capture noisily stem price
capture noisily stem price, lines(2)
capture noisily lv mpg
capture noisily lv price weight
capture noisily smooth 3RSSH, twice mpg, generate(smpg)
capture noisily smooth 4253H price, generate(sprice)
capture noisily list mpg smpg price sprice in 1/10

* --- graphs (the log shows only the messages) ---
capture noisily quantile price
capture noisily symplot price
gen sdp = 300
capture noisily serrbar price sdp mpg
capture noisily dotplot mpg
capture noisily sunflower price weight
capture noisily lpoly price weight, nograph generate(xg yg) n(20) bwidth(300)
capture noisily list xg yg in 1/5
capture noisily cusum price mpg, generate(cs)
capture noisily gladder price

* --- ml ---
program define mylogit
    args lnf xb
    quietly replace `lnf' = ln(invlogit(`xb')) if $ML_y1 == 1
    quietly replace `lnf' = ln(invlogit(-`xb')) if $ML_y1 == 0
end
capture noisily ml model lf mylogit (foreign = mpg weight)
capture noisily ml maximize
capture noisily ml model lf mylogit (foreign = mpg weight), vce(robust)
capture noisily ml maximize, nolog
capture noisily ereturn list
capture noisily mlexp (foreign*ln(invlogit({xb: mpg weight _cons})) + (1-foreign)*ln(invlogit(-{xb:})))
capture noisily logit foreign mpg weight

* --- resampling ---
set seed 1903
capture noisily jackknife r(mean) r(sd): summarize price
capture noisily jackknife: regress price weight
capture noisily jackknife, cluster(rep78): regress price weight
capture noisily permute foreign diff=(r(mu_2)-r(mu_1)), reps(100) nodots: ttest price, by(foreign)
program define sim1, rclass
    drop _all
    set obs 30
    gen x = rnormal()
    summarize x
    return scalar m = r(mean)
end
capture noisily simulate m=r(m), reps(20) nodots: sim1
capture noisily summarize m

* --- model building ---
mkdata
capture noisily mkspline mpg1 20 mpg2 30 mpg3 = mpg
capture noisily summarize mpg1 mpg2 mpg3
capture noisily mkspline wsp = weight, cubic nknots(4) displayknots
capture noisily xi: regress price i.rep78 weight
capture noisily drop _I*
capture noisily nestreg: regress price (weight) (length foreign)
capture noisily nestreg, lr: logit foreign (mpg) (weight price)
capture noisily stepwise, pr(.2): regress price mpg weight length turn displacement
capture noisily stepwise, pe(.1): regress price mpg weight length turn displacement
capture noisily fp <weight>: regress price <weight> mpg
capture noisily mfp: regress price weight mpg length

* --- after estimation ---
capture noisily regress price mpg weight
capture noisily testnl _b[mpg]/_b[weight] = -.01
capture noisily linktest
capture noisily logit foreign mpg weight
capture noisily predictnl p1 = invlogit(xb()), se(sep1)
capture noisily summarize p1 sep1
capture noisily regress price i.rep78
capture noisily contrast rep78
capture noisily contrast r.rep78
capture noisily pwcompare rep78, mcompare(bonferroni)
capture noisily pwcompare rep78, mcompare(tukey) effects
capture noisily pwmean price, over(rep78) mcompare(sidak)
capture noisily logit foreign mpg
capture noisily estimates store a
capture noisily logit foreign weight
capture noisily estimates store b
capture noisily suest a b
capture noisily tobit mpg weight, ll(15)
capture noisily margins, dydx(weight)
capture noisily ologit rep78 mpg weight
capture noisily margins, predict(outcome(3))
capture noisily oprobit rep78 mpg weight
capture noisily margins, dydx(mpg) predict(outcome(5))

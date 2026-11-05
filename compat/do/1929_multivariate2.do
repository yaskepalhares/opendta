* Multivariate statistics, part 3 (block 23a): rotate and rotatemat with the gradient projection
* criteria (vgpf, entropy, tandem1, tandem2, cf(), bentler, oblimax, quartimin, target(),
* partial(), protect()), pca with vce(normal) and means, estat factors and estat rotatecompare,
* the estat subcommands of canon, mvtest covariances with compound, block() and equals(),
* mvtest correlations with equals() and by() (Jennrich), mvtest normality with univariate and
* bivariate, manova with error terms (/ and nested terms) and manovatest with test().
* Every command whose layout or results the manual does not settle runs under capture noisily,
* so that one Stata run shows the whole case.
clear
* the data come from runiform(): OpenDTA does not reproduce Stata's rnormal() (case 0307)
set seed 1929
set obs 200
gen f1 = invnormal(runiform())
gen f2 = invnormal(runiform())
gen g = 1 + mod(_n, 3)
gen h = 1 + mod(_n, 2)
forvalues i = 1/4 {
    gen x`i' = f1 + 0.4*f2 + invnormal(runiform())*.6
}
forvalues i = 5/8 {
    gen x`i' = f2 + 0.3*f1 + invnormal(runiform())*.6
}
gen y1 = x1 + g + invnormal(runiform())
gen y2 = x5 - h + invnormal(runiform())

* ---- rotatemat (the loading matrix of [MV] rotatemat)
matrix L = (601,019,388,221\372,-025,252,132\413,-117,388,144\487,-100,254,192\691,-304,-279,035\690,-409,-200,-076\677,-409,-292,084\674,-189,-099,122\697,-454,-212,-080\476,534,-486,092\558,332,-142,-090\472,508,-139,256\602,244,028,295\423,058,015,-415\394,089,097,-362\510,095,347,-249\466,197,-004,-381\515,312,152,-147\443,089,109,-150\614,-118,126,-038\589,227,057,123\608,-107,127,-038\687,-044,138,098\651,177,-212,-017)
matrix colnames L = F1 F2 F3 F4
matrix L = L/1000
capture noisily rotatemat L, format(%6.3f)
capture noisily return list
capture noisily rotatemat L, quartimax format(%6.3f) norotation blanks(0.3)
capture noisily return list
capture noisily rotatemat L, oblimin oblique format(%6.3f) blanks(0.3)
capture noisily matlist r(T)'*r(T)
capture noisily rotatemat L, promax blanks(0.3) format(%6.3f)
capture noisily rotatemat L, vgpf norotation noloading
capture noisily rotatemat L, entropy norotation noloading
capture noisily rotatemat L, tandem1 norotation noloading
capture noisily rotatemat L, tandem2 norotation noloading
capture noisily rotatemat L, cf(0.25) oblique norotation noloading
capture noisily rotatemat L, bentler oblique norotation noloading
capture noisily rotatemat L, oblimax norotation noloading
capture noisily rotatemat L, quartimin oblique norotation noloading
capture noisily rotatemat L, entropy oblique
capture noisily rotatemat L, cf(2)
capture noisily rotatemat L, varimax normalize nodisplay
capture noisily return list

* ---- rotate after factor and pca
factor x1-x8, ipf factors(2)
capture noisily rotate, oblimin
capture noisily ereturn list
capture noisily rotate, cf(0.5) oblique
capture noisily rotate, bentler oblique
capture noisily rotate, oblimax
capture noisily rotate, entropy
matrix W = (1,0\1,0\1,0\1,0\0,1\0,1\0,1\0,1)
capture noisily rotate, target(W) normalize
matrix Wt = J(8,2,1)
capture noisily rotate, partial(W Wt) oblique
set seed 123
capture noisily rotate, oblimin oblique normalize protect(5)
capture noisily estat rotatecompare
capture noisily estat common
capture noisily estat factors
capture noisily return list
capture noisily estat factors, factors(2)
pca x1-x8, comp(3)
capture noisily rotate, tandem1
capture noisily rotate, quartimin oblique
capture noisily estat rotatecompare

* ---- pca with vce(normal) and means
pca x1-x4, comp(2) vce(normal)
capture noisily ereturn list
capture noisily pca x1-x4, comp(2) vce(normal) novce
capture noisily pca x1-x4, means
capture noisily pca x1-x4, comp(2) covariance vce(normal) level(90)

* ---- estat after canon
quietly canon (x1 x2 x3 y1) (x5 x6 x7 y2)
capture noisily estat correlations
capture noisily return list
capture noisily estat loadings
capture noisily estat rotatecompare
capture noisily estat rotate
capture noisily return list
capture noisily estat rotate, stdcoefs format(%9.3f)
capture noisily estat rotate, loadings
capture noisily estat rotatecompare

* ---- mvtest covariances and correlations
mvtest covariances x1 x2 x5, compound
capture noisily return list
mvtest covariances x1 x2 x5 x6, block(x1 x2 || x5 x6)
capture noisily mvtest covariances x1 x2 x5 x6, block(x1 x2)
quietly correlate x1 x2 x5, covariance
matrix S0 = r(C)
mvtest covariances x1 x2 x5, equals(S0)
capture noisily mvtest covariances x1 x2 x5, by(g)
capture noisily return list
capture noisily mvtest covariances x1 x2, compound spherical
matrix C = (1, 0.75, 0 \ 0.75, 1, 0 \ 0, 0, 1)
mvtest correlations x1 x2 x5, equals(C)
capture noisily return list
mvtest correlations x1 x2 x5, by(g)
capture noisily mvtest correlations x1 x2 x5, by(g h)
capture noisily return list

* ---- mvtest normality
mvtest normality x1 x2 x3, univariate bivariate stats(all)
capture noisily return list
capture noisily matrix list r(Utest)
capture noisily matrix list r(Btest)
capture noisily matrix list r(U_dh)
capture noisily by h, sort: mvtest normality x1 x2 x3

* ---- manova with error terms
clear
* the data come from runiform(): OpenDTA does not reproduce Stata's rnormal() (case 0307)
set seed 2923
set obs 48
gen a = 1 + mod(_n, 2)
gen b = 1 + mod(floor((_n - 1) / 2), 3)
gen c = 1 + mod(floor((_n - 1) / 6), 2)
gen u = invnormal(runiform()) + a + 0.3*b
gen v = invnormal(runiform()) - a + 0.2*c
manova u v = a / b|a /
capture noisily ereturn list
capture noisily manova u v = a / b|a / c|b|a /, dropemptycells
capture noisily manovatest a / b|a
capture noisily manovatest a / b|a c|b|a
capture noisily return list
capture noisily manova u v = c / b|c a /

* ---- manovatest with test()
quietly manova u v = b
manovatest, showorder
matrix C = (1,-1,0,0 \ 1,0,-1,0)
capture noisily manovatest, test(C)
capture noisily return list
matrix C1 = (2,-1,-1,0)
capture noisily manovatest, test(C1)
matrix Y1 = (1,-1)
capture noisily manovatest, test(C1) ytransform(Y1)
matrix Cbad = (1,-1,0)
capture noisily manovatest, test(Cbad)

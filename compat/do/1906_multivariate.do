* Phase 20: multivariate statistics (principal components, factor analysis and rotation, scale
* reliability, multivariate tests and models, cluster analysis, discriminant analysis, correspondence
* analysis, multidimensional scaling, procrustes analysis and biplots). Every command runs under
* capture noisily so that a difference in one command does not stop the rest of the case; the
* layouts marked VERIFY in the code are settled by the reference log.
* The data are built from rounded draws, so that the last digit of the draws does not matter.
program define mkmv
    clear
    set seed 1906
    set obs 150
    gen int id = _n
    gen byte g = ceil(3*runiform())
    gen byte g2 = 1 + (runiform() < .5)
    gen double f1 = round(invnormal(runiform()), .01)
    gen double f2 = round(invnormal(runiform()), .01)
    gen double x1 = round(.8*f1 + .6*invnormal(runiform()) + .6*g, .01)
    gen double x2 = round(.7*f1 + .7*invnormal(runiform()) - .4*g, .01)
    gen double x3 = round(.6*f1 + .8*invnormal(runiform()) + .3*g, .01)
    gen double x4 = round(.8*f2 + .6*invnormal(runiform()), .01)
    gen double x5 = round(.7*f2 + .7*invnormal(runiform()), .01)
    gen double x6 = round(.6*f2 + .8*invnormal(runiform()), .01)
    gen double y1 = round(1 + .5*x1 + .3*x4 + invnormal(runiform()), .01)
    gen double y2 = round(.5*x2 - .4*x5 + invnormal(runiform()), .01)
    gen double y3 = round(.3*x3 + .2*g + invnormal(runiform()), .01)
    gen byte a = ceil(3*runiform())
    gen byte b = ceil(4*runiform())
    replace b = a if runiform() < .4
    gen byte c = ceil(2*runiform())
    replace c = 1 if a == 1 & runiform() < .5
end

* ---- principal components
mkmv
capture noisily pca x1 x2 x3 x4 x5 x6
capture noisily pca x1 x2 x3 x4 x5 x6, components(2)
capture noisily pca x1 x2 x3 x4 x5 x6, covariance
capture noisily pca x1 x2 x3 x4 x5 x6, mineigen(1)
capture noisily estat loadings
capture noisily estat kmo
capture noisily estat smc
capture noisily estat correlation
capture noisily estat residuals
capture noisily estat summarize
capture noisily screeplot, nodraw
capture noisily loadingplot, nodraw
capture noisily scoreplot, nodraw
capture noisily rotate
capture noisily rotate, quartimax
capture noisily rotate, promax
capture noisily rotate, oblimin oblique
capture noisily predict p1 p2, score
capture noisily summarize p1 p2
capture noisily ereturn list

* ---- factor analysis
capture noisily factor x1 x2 x3 x4 x5 x6
capture noisily factor x1 x2 x3 x4 x5 x6, pcf factors(2)
capture noisily factor x1 x2 x3 x4 x5 x6, ipf factors(2)
capture noisily factor x1 x2 x3 x4 x5 x6, ml factors(2)
capture noisily estat kmo
capture noisily estat smc
capture noisily estat common
capture noisily estat structure
capture noisily rotate
capture noisily rotate, promax
capture noisily rotate, oblimin oblique
capture noisily predict f1s f2s
capture noisily summarize f1s f2s
capture noisily ereturn list

* ---- scale reliability
capture noisily alpha x1 x2 x3
capture noisily alpha x1 x2 x3 x4 x5 x6, item
capture noisily alpha x1 x2 x3, std
capture noisily alpha x1 x2 x3 x4, casewise detail
capture noisily return list

* ---- multivariate tests and models
capture noisily hotelling x1 x2 x3
capture noisily hotelling x1 x2 x3, by(g2)
capture noisily return list
capture noisily mvreg y1 y2 y3 = x1 x2
capture noisily mvreg y1 y2 y3 = x1 x2, corr
capture noisily predict r1, residuals equation(y1)
capture noisily summarize r1
capture noisily ereturn list
capture noisily manova y1 y2 y3 = g2
capture noisily manova y1 y2 = g g2
capture noisily canon (x1 x2 x3) (y1 y2 y3)
capture noisily ereturn list
capture noisily mvtest means x1 x2 x3
capture noisily mvtest means x1 x2 x3, by(g2)
capture noisily mvtest covariances x1 x2 x3, by(g2)
capture noisily mvtest covariances x1 x2 x3, diagonal
capture noisily mvtest covariances x1 x2 x3, spherical
capture noisily mvtest correlations x1 x2 x3
capture noisily mvtest normality x1 x2 x3
capture noisily return list

* ---- cluster analysis
mkmv
capture noisily cluster kmeans x1 x2 x3, k(3) name(km) start(firstk)
capture noisily tabulate km g
capture noisily cluster kmedians x1 x2 x3, k(3) name(kmed) start(firstk)
capture noisily tabulate kmed g
capture noisily cluster averagelinkage x1 x2 x3, name(av)
capture noisily cluster generate av3 = groups(3), name(av)
capture noisily tabulate av3 g
capture noisily cluster wardslinkage x1 x2 x3, measure(L2squared) name(wd)
capture noisily cluster generate wd3 = groups(3), name(wd)
capture noisily tabulate wd3 g
capture noisily cluster stop wd, rule(calinski) groups(2/6)
capture noisily cluster stop wd, rule(duda) groups(2/5)
capture noisily cluster singlelinkage x1 x2 x3, name(sl)
capture noisily cluster completelinkage x1 x2 x3, name(co)
capture noisily cluster waveragelinkage x1 x2 x3, name(wa)
capture noisily cluster centroidlinkage x1 x2 x3, name(ce)
capture noisily cluster medianlinkage x1 x2 x3, name(me)
capture noisily cluster list
capture noisily cluster list av
capture noisily cluster dir
capture noisily cluster dendrogram av, nodraw
capture noisily cluster rename co co2
capture noisily cluster drop sl
capture noisily cluster dir
capture noisily matrix dissimilarity D1 = x1 x2 x3 in 1/5, L1
capture noisily matrix list D1
capture noisily matrix dissimilarity D2 = x1 x2 x3 in 1/5, L2
capture noisily matrix list D2
capture noisily matrix dissimilarity D3 = x1 x2 x3 in 1/5, correlation
capture noisily matrix list D3

* ---- discriminant analysis
mkmv
capture noisily discrim lda x1 x2 x3, group(g)
capture noisily discrim lda x1 x2 x3, group(g) lootable
capture noisily estat classtable
capture noisily estat classtable, loo
capture noisily estat errorrate
capture noisily estat errorrate, pp
capture noisily estat grdistances
capture noisily estat covariance
capture noisily estat list, misclassified
capture noisily predict cls
capture noisily predict pr1 pr2 pr3, pr
capture noisily predict md1 md2 md3, mahalanobis
capture noisily summarize cls pr1 pr2 pr3 md1 md2 md3
capture noisily ereturn list
capture noisily candisc x1 x2 x3, group(g)
capture noisily ereturn list
capture noisily discrim lda x1 x2 x3, group(g) priors(proportional) notable
capture noisily estat classtable
capture noisily discrim qda x1 x2 x3, group(g) lootable
capture noisily estat covariance
capture noisily discrim knn x1 x2 x3, group(g) k(3) lootable
capture noisily estat errorrate
capture noisily discrim logistic x1 x2 x3, group(g)
capture noisily ereturn list

* ---- correspondence analysis
mkmv
capture noisily ca a b
capture noisily estat table
capture noisily estat coordinates
capture noisily ca a b, dimensions(1) normalize(principal)
capture noisily ereturn list
capture noisily mca a b c
capture noisily estat coordinates
capture noisily ereturn list

* ---- multidimensional scaling, procrustes and biplots
mkmv
capture noisily mds x1 x2 x3 in 1/30, id(id) noplot
capture noisily estat config
capture noisily ereturn list
capture noisily mds x1 x2 x3 in 1/30, id(id) noplot method(modern)
capture noisily matrix dissimilarity DM = x1 x2 x3 in 1/8, L2
capture noisily mdsmat DM, noplot
capture noisily procrustes (x1 x2) (x4 x5)
capture noisily procrustes (x1 x2) (x4 x5), noscale
capture noisily predict pf1 pf2
capture noisily summarize pf1 pf2
capture noisily ereturn list
capture noisily biplot x1 x2 x3 x4, nodraw
capture noisily biplot x1 x2 x3 x4, std nodraw

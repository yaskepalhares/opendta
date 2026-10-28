* Multivariate statistics, second part: clustermat, manovatest, camat, the plots after ca, mca,
* mds and procrustes (cabiplot, caprojection, mcaplot, mcaprojection, mdsconfig, mdsshepard,
* procoverlay), mdslong, the binary similarity measures of cluster, cluster kmeans keepcenters
* and procrustes transform(unrestricted). The plots are drawn with nodraw (only their messages
* and errors are compared). Every command whose behavior the manual does not settle runs under
* capture noisily, so that one Stata run shows the whole case.
clear
matrix D = (0,1,5,6,8 \ 1,0,4,5,7 \ 5,4,0,1.5,3 \ 6,5,1.5,0,2.5 \ 8,7,3,2.5,0)
matrix rownames D = a b c d e
matrix colnames D = a b c d e

* clustermat
clustermat averagelinkage D, labelvar(nm)
list
cluster list
cluster generate g2 = groups(2)
list nm g2
capture noisily clustermat singlelinkage D
clustermat completelinkage D, add name(cmp)
cluster list cmp
capture noisily cluster dendrogram cmp, labels(nm) nodraw
matrix L = (0 \ 1 \ 0 \ 5 \ 4 \ 0 \ 6 \ 5 \ 1.5 \ 0)
capture noisily clustermat singlelinkage L, shape(lower) clear
capture noisily list

* manovatest
clear
set seed 1917
set obs 60
gen a = mod(_n, 3)
gen b = mod(_n, 2)
gen y1 = a + invnormal(runiform())
gen y2 = b + .5*a + invnormal(runiform())
gen y3 = invnormal(runiform())
manova y1 y2 y3 = a b a#b
manovatest a
return list
capture noisily manovatest a b
capture noisily manovatest b / a
matrix M = (1, -1, 0 \ 0, 1, -1)
capture noisily manovatest b, ytransform(M)
capture noisily return list
capture noisily manovatest c
capture noisily manovatest, showorder

* camat and the ca plots
matrix F = (20, 10, 5 \ 8, 15, 12 \ 3, 9, 25 \ 10, 10, 10)
matrix rownames F = r1 r2 r3 r4
matrix colnames F = c1 c2 c3
camat F
ereturn list
camat F, rowname(age) colname(answer) normalize(principal)
capture noisily cabiplot, nodraw
capture noisily cabiplot, origin norow nodraw
capture noisily caprojection, nodraw
capture noisily estat coordinates

* mca plots
clear
set seed 2917
set obs 90
gen a = 1 + mod(_n, 3)
gen b = 1 + mod(floor(_n/2), 4)
gen c = 1 + (runiform() > .5)
quietly mca a b c
capture noisily mcaplot, nodraw
capture noisily mcaplot a b, origin nodraw
capture noisily mcaprojection, nodraw

* mds plots and mdslong
quietly mdsmat D
capture noisily mdsconfig, nodraw
capture noisily mdsshepard, nodraw
clear
input i j d
1 2 1
1 3 5
1 4 6
1 5 8
2 3 4
2 4 5
2 5 7
3 4 1.5
3 5 3
4 5 2.5
end
mdslong d, id(i j)
ereturn list
capture noisily mdsconfig, nodraw
capture noisily mdslong d, id(i j) method(modern)

* binary similarity measures and keepcenters
clear
set seed 3917
set obs 14
forvalues j = 1/5 {
    gen b`j' = runiform() > .5
}
gen x = invnormal(runiform())
gen y = invnormal(runiform())
gen z = .8*x + .4*y + .3*invnormal(runiform())
gen w = y - .5*x + .3*invnormal(runiform())
cluster averagelinkage b1-b5, measure(Jaccard) name(jac)
cluster list jac
cluster generate gj = groups(3), name(jac)
list b1-b5 gj jac_id jac_ord jac_hgt
foreach m in matching Russell Hamann Dice antiDice Sneath Rogers Ochiai Yule Anderberg Kulczynski Pearson Gower2 {
    capture noisily cluster singlelinkage b1-b5, measure(`m') name(s_`m')
    capture noisily cluster generate g_`m' = groups(3), name(s_`m')
}
capture noisily list g_*
capture noisily matrix dissimilarity S = b1-b5, Jaccard
capture noisily matrix list S
cluster kmeans x y, k(3) start(firstk) keepcenters name(km)
list x y km
drop in 15/17

* procrustes transform(unrestricted) and procoverlay
procrustes (z w) (x y), transform(unrestricted)
ereturn list
capture noisily procoverlay, nodraw
quietly procrustes (z w) (x y)
capture noisily procoverlay, nodraw
capture noisily procrustes (z w) (x y), transform(oblique)

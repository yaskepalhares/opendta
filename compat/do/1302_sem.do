* sem: path model, confirmatory factor analysis and a model with two latent variables;
* deterministic data (runiform() reproduces Stata's sequence)
clear
set seed 1302
set obs 300
gen F = invnormal(runiform())
gen G = .5*F + .8*invnormal(runiform())
gen x1 = F + .6*invnormal(runiform())
gen x2 = .8*F + .7*invnormal(runiform())
gen x3 = 1.2*F + .5*invnormal(runiform())
gen y1 = G + .6*invnormal(runiform())
gen y2 = .9*G + .6*invnormal(runiform())
gen y3 = 1.1*G + .6*invnormal(runiform())
gen z = invnormal(runiform())
gen w = .5*z + .3*x1 + invnormal(runiform())
rename F f0
rename G g0
sem (w <- z x1)
ereturn list
sem (F -> x1 x2 x3)
estat gof
estat gof, stats(all)
sem (F -> x1 x2 x3) (G -> y1 y2 y3) (G <- F)
estat gof, stats(all)
sem (F -> x1 x2 x3) (G -> y1 y2 y3) (G <- F), standardized
sem (F -> x1 x2 x3) (G -> y1 y2 y3) (G <- F), covariance(e.x1*e.y1)
predict fhat, latent(F)
predict what, xb(x2)
summarize fhat what

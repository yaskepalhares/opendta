* Mata: moptimize() with lf, d2 and lf2 evaluators, display and results
clear
set seed 408
set obs 200
gen x1 = invnormal(runiform())
gen x2 = runiform()
gen y = runiform() < invlogit(-.5 + x1 - 2*x2)
gen yc = 1 + 2*x1 - x2 + invnormal(runiform())*1.5
logit y x1 x2
regress yc x1 x2
mata:
void lflogit(transmorphic M, real rowvector b, real colvector lnfj) {
    real colvector xb, y
    xb = moptimize_util_xb(M, b, 1)
    y = moptimize_util_depvar(M, 1)
    lnfj = y:*ln(invlogit(xb)) + (1:-y):*ln(invlogit(-xb))
}
M = moptimize_init()
moptimize_init_evaluator(M, &lflogit())
moptimize_init_evaluatortype(M, "lf")
moptimize_init_depvar(M, 1, "y")
moptimize_init_eq_indepvars(M, 1, "x1 x2")
moptimize_init_search(M, "off")
moptimize(M)
moptimize_result_display(M)
moptimize_result_coefs(M)
void d2logit(transmorphic M, real scalar todo, real rowvector b, lnf, g, H) {
    real colvector xb, y, p
    xb = moptimize_util_xb(M, b, 1)
    y = moptimize_util_depvar(M, 1)
    p = invlogit(xb)
    lnf = moptimize_util_sum(M, y:*ln(p) + (1:-y):*ln(1:-p))
    if (todo >= 1) {
        g = moptimize_util_vecsum(M, 1, y - p, lnf)
        if (todo == 2) H = -moptimize_util_matsum(M, 1, 1, p:*(1:-p), lnf)
    }
}
M2 = moptimize_init()
moptimize_init_evaluator(M2, &d2logit())
moptimize_init_evaluatortype(M2, "d2")
moptimize_init_depvar(M2, 1, "y")
moptimize_init_eq_indepvars(M2, 1, "x1 x2")
moptimize_init_tracelevel(M2, "none")
moptimize(M2)
moptimize_result_coefs(M2)
sqrt(diagonal(moptimize_result_V(M2)))'
void lf2norm(transmorphic M, real scalar todo, real rowvector b, lnfj, S, H) {
    real colvector mu, lns, y, z, s
    real matrix h11, h12, h22
    real scalar lnf
    mu = moptimize_util_xb(M, b, 1)
    lns = moptimize_util_xb(M, b, 2)
    y = moptimize_util_depvar(M, 1)
    s = exp(lns) :* J(rows(y), 1, 1)
    z = (y - mu):/s
    lnfj = ln(normalden(z)) - ln(s)
    if (todo >= 1) {
        S = (z:/s, z:*z :- 1)
        if (todo == 2) {
            lnf = moptimize_util_sum(M, lnfj)
            h11 = -moptimize_util_matsum(M, 1, 1, 1:/(s:^2), lnf)
            h12 = -moptimize_util_matsum(M, 1, 2, 2*z:/s, lnf)
            h22 = -moptimize_util_matsum(M, 2, 2, 2*z:*z, lnf)
            H = (h11, h12 \ h12', h22)
        }
    }
}
M3 = moptimize_init()
moptimize_init_evaluator(M3, &lf2norm())
moptimize_init_evaluatortype(M3, "lf2")
moptimize_init_depvar(M3, 1, "yc")
moptimize_init_eq_indepvars(M3, 1, "x1 x2")
moptimize_init_eq_name(M3, 1, "mu")
moptimize_init_eq_name(M3, 2, "lnsigma")
moptimize_init_search(M3, "off")
moptimize(M3)
moptimize_result_display(M3)
end

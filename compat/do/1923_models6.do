* Models, part 6 (block 20d): the program-based forms of nl, nlsur and gmm, and gmm with
* parameters(). No reference until the Stata run. Every command whose behavior the manual does
* not settle runs under capture noisily, so that one Stata run shows the whole case.
clear
set obs 400
set seed 1923
gen x1 = round(invnormal(runiform()), .01)
gen x2 = round(invnormal(runiform()), .01)
gen z = round(invnormal(runiform()), .01)
gen xp = round(runiform()*4 + 1, .01)
gen yp = round(2*xp^.7 + invnormal(runiform())*.2, .01)
gen y = round(exp(.3 + .5*x1 - .4*x2) + invnormal(runiform())*.5, .01)
gen y1 = round(1 + 2*x1 + invnormal(runiform())*.5, .01)
gen y2 = round(-1 + .5*x1 + 3*x2 + invnormal(runiform())*.5, .01)
capture program drop nlpow
program nlpow
    syntax varlist(min=2 max=2) if, at(name)
    local y : word 1 of `varlist'
    local x : word 2 of `varlist'
    replace `y' = `at'[1,1]*`x'^`at'[1,2] `if'
end
capture program drop nlpow2
program nlpow2, rclass
    syntax varlist(min=2 max=2) [if]
    local y : word 1 of `varlist'
    local x : word 2 of `varlist'
    return local eq "`y' = {a}*`x'^{b=1}"
    return local title "Power curve"
end
capture program drop gmm_poi
program gmm_poi
    syntax varlist if, at(name)
    replace `varlist' = y - exp(`at'[1,1]*x1 + `at'[1,2]*x2 + `at'[1,3]) `if'
end
capture program drop nlsurlin
program nlsurlin
    syntax varlist(min=3 max=3) if, at(name)
    tokenize `varlist'
    replace `1' = `at'[1,1] + `at'[1,2]*`3' `if'
    replace `2' = `at'[1,3] + `at'[1,2]*`3' `if'
end
capture program drop nlsurlin2
program nlsurlin2, rclass
    syntax varlist(min=3 max=3) [if]
    tokenize `varlist'
    return scalar n_eq = 2
    return local eq_1 "`1' = {a1} + {b}*`3'"
    return local eq_2 "`2' = {a2} + {b}*`3'"
end

* nl: function evaluator (parameters() and nparameters()) and substitutable expression program
capture noisily nl pow @ yp xp, parameters(a b) initial(b 1)
capture noisily ereturn list
capture noisily nl pow @ yp xp, nparameters(2) initial(b2 1)
capture noisily nl pow2 : yp xp
capture noisily ereturn list
capture noisily predict double yphat
capture noisily summarize yphat

* gmm: moment-evaluator program; parameters() with substitutable expressions
capture noisily gmm gmm_poi, nequations(1) instruments(x1 x2 z) parameters(b1 b2 b0)
capture noisily ereturn list
capture noisily gmm gmm_poi, equations(m1) instruments(m1: x1 x2 z) nparameters(3) onestep
capture noisily gmm (y - exp({b1}*x1 + {b2}*x2 + {b0})), instruments(x1 x2 z) parameters(b0 b1 b2)

* nlsur: function evaluator and substitutable expression program
capture noisily nlsur lin @ y1 y2 x1, nequations(2) parameters(a1 b a2)
capture noisily ereturn list
capture noisily nlsur lin2 : y1 y2 x1, nls

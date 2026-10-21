{smcl}
{title:Title}

{p 4 4 2}
{bf:cnsreg}, {bf:constraint}, {bf:sureg}, {bf:reg3}, {bf:eivreg}, {bf:vwls}, {bf:rreg} {hline 2} Systems of equations, constraints and robust or weighted regression

{title:Syntax}

{p 8 16 2}
{cmd:constraint} [{cmd:define}] # {it:exp} | {cmd:dir} | {cmd:list} | {cmd:drop} {it:numlist}|{cmd:_all}

{p 8 16 2}
{cmd:cnsreg} depvar indepvars [if] [in] [weight] {cmd:, constraints(}{it:numlist}{cmd:)}

{p 8 16 2}
{cmd:sureg} {cmd:(}depvar1 varlist1{cmd:)} {cmd:(}depvar2 varlist2{cmd:)} ... [if] [in] [{cmd:,} {cmd:isure} {cmd:small} {cmd:dfk}]

{p 8 16 2}
{cmd:reg3} {cmd:(}depvar1 varlist1{cmd:)} {cmd:(}depvar2 varlist2{cmd:)} ... [if] [in] [{cmd:,} {cmd:2sls}|{cmd:3sls}|{cmd:ols}|{cmd:sure} {cmd:inst(}varlist{cmd:)} {cmd:small}]

{p 8 16 2}
{cmd:eivreg} depvar indepvars [if] [in] [{cmd:,} {cmd:r(}{it:varname # ...}{cmd:)}]

{p 8 16 2}
{cmd:vwls} depvar indepvars [if] [in] {cmd:, sd(}varname{cmd:)}

{p 8 16 2}
{cmd:rreg} depvar indepvars [if] [in] [{cmd:,} {cmd:tune(}#{cmd:)} {cmd:genwt(}newvar{cmd:)}]

{title:Description}

{pstd}
{cmd:constraint} stores linear constraints; {cmd:cnsreg} fits least squares under them. {cmd:sureg} is Zellner's seemingly unrelated regression and {cmd:reg3} the two- and three-stage least squares of a system. {cmd:eivreg} corrects for measurement error with known reliabilities, {cmd:vwls} is weighted least squares with known variances and {cmd:rreg} is iteratively reweighted robust regression (Huber then biweight, after dropping the observations with Cook's D above 1).

{pstd}
Not available: {cmd:constraints()} in {cmd:sureg} and {cmd:reg3}. The table layouts are marked for checking against the Stata 14 log (case 1908).

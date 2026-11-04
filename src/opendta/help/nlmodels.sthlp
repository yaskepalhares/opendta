{smcl}
{title:Title}

{p 4 4 2}
{bf:nl}, {bf:nlsur}, {bf:gmm}, {bf:frontier} {hline 2} Nonlinear least squares, systems, GMM and stochastic frontier

{title:Syntax}

{p 8 16 2}
{cmd:nl} {cmd:(}depvar {cmd:=} {it:exp}{cmd:)} [if] [in] [weight] [{cmd:,} {cmd:variables(}varlist{cmd:)} {cmd:initial(}{it:p # ...}{cmd:)}]

{p 8 16 2}
{cmd:nlsur} {cmd:(}depvar1 {cmd:=} {it:exp1}{cmd:)} {cmd:(}depvar2 {cmd:=} {it:exp2}{cmd:)} ... [if] [in] [{cmd:,} {cmd:fgnls} {cmd:ifgnls}]

{p 8 16 2}
{cmd:gmm} {cmd:(}{it:exp}{cmd:)} ... [if] [in] {cmd:, instruments(}varlist [{cmd:,} noconstant]{cmd:)} [{cmd:onestep}|{cmd:twostep} {cmd:winitial(}{cmd:identity}|{cmd:unadjusted}{cmd:)} {cmd:wmatrix(}{cmd:unadjusted}|{cmd:robust}|{cmd:cluster} clustvar{cmd:)}]

{p 8 16 2}
{cmd:nl} {it:func} {cmd:@} depvar [varlist] [if] [in] [weight]{cmd:,} {cmd:parameters(}namelist{cmd:)}|{cmd:nparameters(}#{cmd:)} [{it:options}]

{p 8 16 2}
{cmd:nl} {it:sexp} {cmd::} depvar [varlist] [if] [in] [weight] [{cmd:,} {it:options}]

{p 8 16 2}
{cmd:nlsur} {it:func} {cmd:@} depvar1 ... [varlist]{cmd:,} {cmd:nequations(}#{cmd:)} {cmd:parameters(}namelist{cmd:)}|{cmd:nparameters(}#{cmd:)} [{it:options}]

{p 8 16 2}
{cmd:gmm} {it:moment_prog} [if] [in]{cmd:,} {cmd:nequations(}#{cmd:)}|{cmd:equations(}namelist{cmd:)} {cmd:parameters(}namelist{cmd:)}|{cmd:nparameters(}#{cmd:)} {cmd:instruments(}...{cmd:)} [{it:options}]

{p 8 16 2}
{cmd:frontier} depvar indepvars [if] [in] [{cmd:,} {cmd:distribution(}{cmd:hnormal}|{cmd:exponential}|{cmd:tnormal}{cmd:)} {cmd:cost} {cmd:constraints(}numlist{cmd:)}]

{title:Description}

{pstd}
Parameters are written {cmd:{c -(}}{it:name}{cmd:{c )-}}, {cmd:{c -(}}{it:name}{cmd:=}{it:init}{cmd:{c )-}} or {cmd:{c -(}}{it:lc}{cmd::} {it:varlist}{cmd:{c )-}} inside the expressions, which are evaluated with the ordinary expression evaluator; derivatives are numerical.

{pstd}
The function evaluator programs ({cmd:nl}{it:func}, {cmd:nlsur}{it:func} and the {cmd:gmm} moment program) are called as {it:prog newvars varlist} {cmd:if} {it:touse}{cmd:, at(}{it:matrix}{cmd:)} and fill the new variables; the substitutable expression programs ({cmd:nl}{it:sexp}, {cmd:nlsur}{it:sexp}) return the equations in {cmd:r(eq)} (or {cmd:r(n_eq)} and {cmd:r(eq_1)}, ...).

{pstd}
Not available: analytic derivatives from the programs, and {cmd:constraints()} with {cmd:cm()} in {cmd:frontier} ({cmd:uhet(}varlist{cmd:)}, {cmd:vhet(}varlist{cmd:)} and, with {cmd:distribution(tnormal)}, {cmd:cm(}varlist{cmd:)} make the variances and the mean of the inefficiency linear in their own variables).

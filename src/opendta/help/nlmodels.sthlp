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
{cmd:frontier} depvar indepvars [if] [in] [{cmd:,} {cmd:distribution(}{cmd:hnormal}|{cmd:exponential}|{cmd:tnormal}{cmd:)} {cmd:cost}]

{title:Description}

{pstd}
Parameters are written {cmd:{c -(}}{it:name}{cmd:{c )-}}, {cmd:{c -(}}{it:name}{cmd:=}{it:init}{cmd:{c )-}} or {cmd:{c -(}}{it:lc}{cmd::} {it:varlist}{cmd:{c )-}} inside the expressions, which are evaluated with the ordinary expression evaluator; derivatives are numerical.

{pstd}
Not available: the program-based forms of {cmd:nl}, {cmd:nlsur} and {cmd:gmm}, {cmd:derivatives()}, {cmd:parameters()} and {cmd:igmm} in {cmd:gmm}, and {cmd:constraints()} in {cmd:frontier} ({cmd:uhet(}varlist{cmd:)}, {cmd:vhet(}varlist{cmd:)} and, with {cmd:distribution(tnormal)}, {cmd:cm(}varlist{cmd:)} make the variances and the mean of the inefficiency linear in their own variables).

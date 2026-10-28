{smcl}
{title:Title}

{p 4 4 2}
{bf:ml}, {bf:mlexp}, {bf:mleval}, {bf:mlsum}, {bf:mlvecsum}, {bf:mlmatsum} {hline 2} Maximum likelihood estimation

{title:Syntax}

{p 8 16 2}
{cmd:ml model} {it:method progname} {cmd:(}[{it:eqname}{cmd::}]{it:varlist}{cmd:)} [{cmd:(}...{cmd:)}] [{cmd:/}{it:name}] [{it:if}] [{it:in}] [{cmd:,} {it:options}]

{p 8 16 2}
{cmd:ml} {c -(}{cmd:init}|{cmd:search}|{cmd:maximize}|{cmd:display}|{cmd:clear}|{cmd:score}{c )-} [{it:...}]

{p 8 16 2}
{cmd:mlexp (}{it:depvar_expression}{cmd:)} [{it:if}] [{it:in}] [{cmd:,} {it:options}]

{title:Description}

{pstd}
{cmd:ml model} accepts the methods {cmd:lf}, {cmd:lf0}, {cmd:lf1}, {cmd:lf2}, {cmd:d0}, {cmd:d1}, {cmd:d2}, {cmd:gf0}, {cmd:gf1} and {cmd:gf2}. The {cmd:lf0}-{cmd:lf2} and {cmd:gf} evaluators are called as {it:prog todo b lnfj} followed by the score variables ({cmd:lf1}, {cmd:lf2}) or the score matrix ({cmd:gf1}, {cmd:gf2}) and the Hessian ({cmd:lf2}, {cmd:gf2}). {cmd:constraints()} applies linear constraints defined with {cmd:constraint}; {cmd:svy} uses the {cmd:svyset} design weights and the linearized design variance. Derivatives that the evaluator does not supply are computed numerically (fourth-order differences). The variance estimators are {cmd:oim}, {cmd:opg}, {cmd:robust} and {cmd:cluster()}. {cmd:ml init}, {cmd:ml search}, {cmd:ml maximize} and {cmd:ml display} follow the usual workflow; the initial-value lines {cmd:initial}, {cmd:feasible} and {cmd:rescale} are accepted by {cmd:ml search}.

{pstd}
Inside evaluators, {cmd:mleval}, {cmd:mlsum}, {cmd:mlvecsum} and {cmd:mlmatsum} work as documented in the manual. {cmd:mlexp} fits a model from a log-likelihood expression with linear combinations {cmd:{c -(}xb: vars{c )-}} and scalar parameters {cmd:{c -(}name=#{c )-}}.

{title:Limitations}

{pstd}
The commands {cmd:ml check}, {cmd:query}, {cmd:report}, {cmd:graph}, {cmd:plot}, {cmd:trace}, {cmd:footnote} and {cmd:count} give error 198, and so do {cmd:constraints()} and {cmd:svy} with {cmd:mlexp}.

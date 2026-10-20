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
{cmd:ml model} accepts the methods {cmd:lf}, {cmd:d0}, {cmd:d1} and {cmd:d2}. Derivatives that the evaluator does not supply are computed numerically (fourth-order differences). The variance estimators are {cmd:oim}, {cmd:opg}, {cmd:robust} and {cmd:cluster()}. {cmd:ml init}, {cmd:ml search}, {cmd:ml maximize} and {cmd:ml display} follow the usual workflow; the initial-value lines {cmd:initial}, {cmd:feasible} and {cmd:rescale} are accepted by {cmd:ml search}.

{pstd}
Inside evaluators, {cmd:mleval}, {cmd:mlsum}, {cmd:mlvecsum} and {cmd:mlmatsum} work as documented in the manual. {cmd:mlexp} fits a model from a log-likelihood expression with linear combinations {cmd:{c -(}xb: vars{c )-}} and scalar parameters {cmd:{c -(}name=#{c )-}}.

{title:Limitations}

{pstd}
The methods {cmd:lf0}, {cmd:lf1}, {cmd:lf2} and {cmd:gf*}; the commands {cmd:ml check}, {cmd:query}, {cmd:report}, {cmd:graph}, {cmd:plot}, {cmd:trace}, {cmd:footnote} and {cmd:count}; constraints and {cmd:svy} give error 198.

{smcl}
{title:Title}

{p 4 4 2}
{bf:bootstrap} {hline 2} Bootstrap sampling and estimation

{title:Syntax}

{p 8 16 2}
{cmd:bootstrap} [{it:exp_list}] [{cmd:,} {cmd:reps(}{it:#}{cmd:)} {cmd:seed(}{it:#}{cmd:)}
{cmd:strata(}{it:varlist}{cmd:)} {cmd:cluster(}{it:varlist}{cmd:)} {cmd:size(}{it:#}{cmd:)}
{cmd:level(}{it:#}{cmd:)} {cmd:nodots} {cmd:noisily}] {cmd::} {it:command}

{p 8 16 2}
{it:estimation_command} ...{cmd:, vce(bootstrap} [{cmd:,} {it:bootstrap_options}]{cmd:)}

{p 8 16 2}
{cmd:estat bootstrap} [{cmd:,} {cmd:normal} {cmd:percentile} {cmd:bc} {cmd:all}]

{pstd}
{it:exp_list} contains {cmd:_b}, {cmd:_se}, {it:exp}, {cmd:(}{it:exp}{cmd:)},
{it:name}{cmd:=}{it:exp} or {cmd:(}{it:name}{cmd::} {it:exp}{cmd:)}. Without
{it:exp_list}, an estimation command is bootstrapped with {cmd:_b}.

{title:Description}

{pstd}
{cmd:bootstrap} runs {it:command} on the data once, then {cmd:reps()} times
(default 50) on samples drawn with replacement (whole clusters with
{cmd:cluster()}, within strata with {cmd:strata()}). The standard errors are
the standard deviations of the replicates; confidence intervals are
normal-based. {cmd:estat bootstrap} adds percentile and bias-corrected
intervals.

{title:OpenDTA notes}

{pstd}
The resampling uses OpenDTA's random-number generator. With the same seed,
the replicates differ from Stata's (the order of Stata's draws in
{cmd:bsample} is not documented), so bootstrap standard errors agree with
Stata's only up to simulation error.

{smcl}
{title:Title}

{p 4 4 2}
{bf:jackknife}, {bf:permute}, {bf:simulate} {hline 2} Jackknife, permutation and Monte Carlo simulation

{title:Syntax}

{p 8 16 2}
{cmd:jackknife} {it:exp_list} [{cmd:,} {cmd:cluster(}{it:varname}{cmd:)} {cmd:mse} {cmd:dof(}#{cmd:)} {cmd:saving(}{it:filename}{cmd:)} {cmd:level(}#{cmd:)} {cmd:nodots} {cmd:noisily}] {cmd::} {it:command}

{p 8 16 2}
{cmd:permute} {it:permvar exp_list} [{cmd:,} {cmd:reps(}#{cmd:)} {cmd:strata(}{it:varlist}{cmd:)} {cmd:saving(}{it:filename}{cmd:)} {cmd:left} {cmd:right} {cmd:nodots}] {cmd::} {it:command}

{p 8 16 2}
{cmd:simulate} {it:exp_list} {cmd:,} {cmd:reps(}#{cmd:)} [{cmd:seed(}#{cmd:)} {cmd:saving(}{it:filename}{cmd:)} {cmd:nodots}] {cmd::} {it:command}

{title:Description}

{pstd}
{cmd:jackknife} runs {it:command} leaving out one observation (or one cluster) at a time and reports jackknife standard errors, V = (n-1)/n sum (theta_i - mean)^2. When {it:command} is an estimation command and {it:exp_list} is omitted, the whole e(b) is used and t statistics use e(df_r). {cmd:permute} permutes {it:permvar} and reports Monte Carlo p-values with Clopper-Pearson confidence intervals. {cmd:simulate} runs {it:command} {cmd:reps()} times and leaves the replicates as the data in memory.

{title:Limitations}

{pstd}
Layouts that the manual does not show are marked VERIFY in the code and checked in compat/do/1903_base_stats.do. Weights are not passed through.

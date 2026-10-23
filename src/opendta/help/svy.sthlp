{smcl}
{title:Title}

{p 4 4 2}
{bf:svyset} {hline 2} Declare survey design for dataset{break}
{bf:svy} {hline 2} The survey prefix command

{title:Syntax}

{p 8 16 2}
{cmd:svyset} [{it:psu}] [{it:weight}] [{cmd:,} {cmd:strata(}{it:varname}{cmd:)}
{cmd:fpc(}{it:varname}{cmd:)} {cmd:vce(}{it:vcetype}{cmd:)} {cmd:singleunit(}{it:method}{cmd:)}
{cmd:brrweight(}{it:varlist}{cmd:)} {cmd:fay(}{it:#}{cmd:)} {cmd:bsrweight(}{it:varlist}{cmd:)} {cmd:bsn(}{it:#}{cmd:)}
{cmd:jkrweight(}{it:varlist}{cmd:)} {cmd:sdrweight(}{it:varlist}{cmd:)} {cmd:mse} {cmd:dof(}{it:#}{cmd:)}]
[{cmd:||} {it:ssu}{cmd:,} {cmd:strata(}{it:varname}{cmd:)} {cmd:fpc(}{it:varname}{cmd:)} ...]

{p 8 16 2}
{cmd:svy} [{cmd:,} {cmd:subpop(}{it:varname} [{it:if}]{cmd:)} {cmd:vce(}{it:vcetype}{cmd:)} {cmd:mse}]
{cmd::} {it:command}

{p 8 16 2}
{cmd:estat effects}

{pstd}
{it:vcetype} is {cmd:linearized} (default), {cmd:jackknife}, {cmd:brr},
{cmd:bootstrap} or {cmd:sdr}. {it:command} is one of {cmd:mean},
{cmd:proportion}, {cmd:total}, {cmd:ratio}, {cmd:tabulate}, {cmd:regress},
{cmd:logit}, {cmd:logistic}, {cmd:probit}, {cmd:cloglog}, {cmd:poisson}, {cmd:ologit},
{cmd:oprobit}, {cmd:mlogit}, {cmd:mprobit}, {cmd:slogit}, {cmd:clogit}, {cmd:nbreg}, {cmd:gnbreg},
{cmd:zip}, {cmd:zinb}, {cmd:glm}, {cmd:tobit}, {cmd:intreg}, {cmd:truncreg}, {cmd:heckman},
{cmd:heckprobit}, {cmd:hetprobit}, {cmd:scobit}, {cmd:biprobit}, {cmd:fracreg}
and {cmd:ivregress}.

{title:Description}

{pstd}
{cmd:svyset} declares the sampling weights, strata, primary sampling units
and finite population corrections; stages after {cmd:||} describe later
sampling stages. Later stages add variance only when the earlier stage is
sampled without replacement (an FPC is set).

{pstd}
With {cmd:vce(linearized)}, variances come from Taylor linearization. The
replication methods recompute the estimates for each replicate:
{cmd:jackknife} deletes one PSU at a time (or uses {cmd:jkrweight()});
{cmd:brr} uses a Hadamard matrix with two PSUs per stratum (or
{cmd:brrweight()}), with Fay's adjustment {cmd:fay()}; {cmd:bootstrap} and
{cmd:sdr} use the replicate weights in {cmd:bsrweight()} and
{cmd:sdrweight()}. With {cmd:mse}, deviations are taken from the full-sample
estimate.

{pstd}
{cmd:estat effects} reports the design effects DEFF (variance relative to
simple random sampling without replacement) and DEFT (square root of the
variance relative to sampling with replacement).

{title:OpenDTA notes}

{pstd}
BRR without replicate weights uses Hadamard matrices of order 4, 8, 16...
Stata may use other orders (12, 20...), giving a different number of
replicates for some numbers of strata.

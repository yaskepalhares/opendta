{smcl}
{title:Title}

{p 4 4 2}
{bf:sts} {hline 2} Survivor and cumulative hazard functions

{title:Syntax}

{p 8 16 2}
{cmd:sts} [{cmd:graph}] [if] [in] [, {cmd:by(}varlist{cmd:)} {cmd:failure} {cmd:na} {cmd:ci} {cmd:level(}#{cmd:)} graph_options]

{p 8 16 2}
{cmd:sts list} [if] [in] [, {cmd:by(}varlist{cmd:)} {cmd:at(}numlist{cmd:)} {cmd:failure} {cmd:level(}#{cmd:)}]

{p 8 16 2}
{cmd:sts generate} newvar {cmd:=} {c -(}{cmd:s}|{cmd:se(s)}|{cmd:lb(s)}|{cmd:ub(s)}|{cmd:na}|{cmd:se(na)}|{cmd:lb(na)}|{cmd:ub(na)}|{cmd:n}|{cmd:d}|{cmd:h}{c )-} [...] [, {cmd:by(}varlist{cmd:)} {cmd:level(}#{cmd:)}]

{p 8 16 2}
{cmd:sts test} varlist [if] [in] [, {cmd:logrank} {cmd:wilcoxon} {cmd:tware} {cmd:peto} {cmd:strata(}varlist{cmd:)}]

{title:Description}

{pstd}
{cmd:sts} estimates the survivor function by the Kaplan-Meier product-limit method (with delayed entry), with the Greenwood standard error and a confidence interval on the ln(-ln S) scale, and the Nelson-Aalen cumulative hazard. {cmd:sts graph} draws the step function; {cmd:sts list} lists it; {cmd:sts generate} creates variables; {cmd:sts test} compares groups with the log-rank (default), Wilcoxon (Breslow), Tarone-Ware or Peto-Peto test.

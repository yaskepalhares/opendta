{smcl}
{title:Title}

{p 4 4 2}
{bf:testnl}, {bf:predictnl}, {bf:linktest}, {bf:contrast}, {bf:pwcompare}, {bf:pwmean}, {bf:suest} {hline 2} Postestimation tests, contrasts and comparisons

{title:Syntax}

{p 8 16 2}
{cmd:testnl} {cmd:(}{it:exp}{cmd:=}{it:exp}{cmd:)} [{cmd:(}...{cmd:)}] [{cmd:,} {cmd:df(}#{cmd:)}]

{p 8 16 2}
{cmd:predictnl} {it:newvar} {cmd:=} {it:exp} [{cmd:,} {cmd:se(}{it:newvar}{cmd:)} {cmd:variance(}{it:newvar}{cmd:)}]

{p 8 16 2}
{cmd:linktest} [{cmd:,} {cmd:nopredict}]

{p 8 16 2}
{cmd:contrast} {it:termlist} [{cmd:,} {cmd:effects} {cmd:asbalanced} {cmd:asobserved} {cmd:mcompare(}{it:method}{cmd:)}]

{p 8 16 2}
{cmd:pwcompare} {it:termlist} [{cmd:,} {cmd:mcompare(}{c -(}{cmd:noadjust}|{cmd:bonferroni}|{cmd:sidak}|{cmd:scheffe}|{cmd:tukey}{c )-}{cmd:)} {cmd:effects} {cmd:cieffects}]

{p 8 16 2}
{cmd:pwmean} {it:varname} {cmd:,} {cmd:over(}{it:varlist}{cmd:)} [{cmd:mcompare(}{it:method}{cmd:)}]

{p 8 16 2}
{cmd:suest} {it:name} [{it:name} ...] [{cmd:,} {cmd:vce(cluster} {it:clustvar}{cmd:)}]

{title:Description}

{pstd}
{cmd:testnl} tests nonlinear hypotheses with the delta method. {cmd:predictnl} makes nonlinear predictions with delta-method standard errors. {cmd:linktest} refits the model with the prediction and its square. {cmd:contrast} accepts the operators {cmd:r.}, {cmd:a.}, {cmd:ar.}, {cmd:g.}, {cmd:h.} and {cmd:j.}. {cmd:suest} combines estimates stored by {cmd:regress}, {cmd:logit}, {cmd:probit} and {cmd:poisson}.

{title:Limitations}

{pstd}
{cmd:pwmean, unequal} gives error 198. Headers and layouts not shown in the manual are marked VERIFY.

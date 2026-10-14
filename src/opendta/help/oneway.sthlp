{smcl}
{title:Title}

{p 4 4 2}
{bf:oneway} {hline 2} One-way analysis of variance{break}
{bf:anova} {hline 2} Analysis of variance and covariance

{title:Syntax}

{p 8 16 2}
{cmd:oneway} {it:response_var} {it:factor_var} [{it:if}] [{it:in}] [{it:weight}]
[{cmd:,} {cmd:tabulate} {cmd:means} {cmd:standard} {cmd:freq} {cmd:obs} {cmd:noanova}]

{p 8 16 2}
{cmd:anova} {it:depvar} {it:termlist} [{it:if}] [{it:in}] [{it:weight}]
[{cmd:,} {cmd:continuous(}{it:varlist}{cmd:)}]

{p 8 16 2}
{cmd:tabulate} {it:var1} [{it:if}] [{it:in}] [{it:weight}]{cmd:,}
{cmd:summarize(}{it:var2}{cmd:)} [{cmd:means} {cmd:standard} {cmd:freq} {cmd:obs}]

{title:Description}

{pstd}
{cmd:oneway} reports the one-way ANOVA table and Bartlett's test for equal
variances. With {cmd:tabulate}, a table of means, standard deviations and
frequencies by group comes first.

{pstd}
{cmd:anova} fits the model with the terms of {it:termlist}. Terms are
variables (categorical by default; {cmd:c.} or {cmd:continuous()} make them
continuous), interactions {cmd:a#b} and factorial expansions {cmd:a##b}.
The table shows partial sums of squares: each term is tested with all the
others in the model. Coefficients are stored in {cmd:e(b)} and {cmd:e(V)} as
in {cmd:regress}, so {cmd:predict}, {cmd:test} and {cmd:margins} work
afterwards.

{title:OpenDTA notes}

{pstd}
{cmd:anova, sequential} and {cmd:repeated()} are not implemented yet.

{title:Stored results}

{pstd}
{cmd:oneway}: {cmd:r(N)}, {cmd:r(F)}, {cmd:r(df_m)}, {cmd:r(df_r)},
{cmd:r(mss)}, {cmd:r(rss)}, {cmd:r(chi2bart)}, {cmd:r(df_bart)}.{break}
{cmd:anova}: {cmd:e(N)}, {cmd:e(mss)}, {cmd:e(rss)}, {cmd:e(df_m)},
{cmd:e(df_r)}, {cmd:e(F)}, {cmd:e(r2)}, {cmd:e(r2_a)}, {cmd:e(rmse)} and, for
each term {it:k}, {cmd:e(term_}{it:k}{cmd:)}, {cmd:e(ss_}{it:k}{cmd:)},
{cmd:e(df_}{it:k}{cmd:)}, {cmd:e(F_}{it:k}{cmd:)}.

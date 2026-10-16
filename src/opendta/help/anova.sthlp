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
{cmd:anova} {it:depvar} {it:termlist} [{cmd:/}] [{it:termlist} {cmd:/} ...] [{it:if}] [{it:in}] [{it:weight}]
[{cmd:,} {cmd:partial} {cmd:sequential} {cmd:continuous(}{it:varlist}{cmd:)}
{cmd:repeated(}{it:varlist}{cmd:)} {cmd:bse(}{it:term}{cmd:)} {cmd:bseunit(}{it:varname}{cmd:)}
{cmd:grouping(}{it:varlist}{cmd:)}]

{p 8 16 2}
{cmd:tabulate} {it:var1} [{it:var2}] [{it:if}] [{it:in}] [{it:weight}]{cmd:,}
{cmd:summarize(}{it:var3}{cmd:)} [{cmd:means} {cmd:standard} {cmd:freq} {cmd:obs}]

{title:Description}

{pstd}
{cmd:oneway} reports the one-way ANOVA table and Bartlett's test for equal
variances. With {cmd:tabulate}, a table of means, standard deviations and
frequencies by group comes first.

{pstd}
{cmd:anova} fits the model with the terms of {it:termlist}. Terms are
variables (categorical by default; {cmd:c.} or {cmd:continuous()} make them
continuous), interactions {cmd:a#b} and factorial expansions {cmd:a##b}.
{cmd:a|b} is {cmd:a} nested within {cmd:b}. The table shows partial sums of
squares (each term tested with all the others in the model) or, with
{cmd:sequential}, sequential sums of squares (each term added in the order
written). Coefficients are stored in {cmd:e(b)} and {cmd:e(V)} as in
{cmd:regress}, so {cmd:predict}, {cmd:test} and {cmd:margins} work afterwards.

{pstd}
A slash names an error term: in {cmd:anova y t / c|t}, {cmd:t} is tested
against {cmd:c|t} instead of the residual. A slash at the end,
{cmd:anova y t / c|t /}, also tests the error term against the residual.

{pstd}
{cmd:repeated(}{it:varlist}{cmd:)} declares repeated-measures variables. After
the table, {cmd:anova} shows the between-subjects error term and, for each
repeated variable, the Huynh-Feldt, Greenhouse-Geisser and Box conservative
corrections, with the p-values of the terms that contain it. {cmd:bse()} names
the between-subjects error term when it cannot be found automatically.

{pstd}
{cmd:tabulate} {it:var1} {it:var2}{cmd:, summarize(}{it:var3}{cmd:)} shows the
mean, standard deviation and frequency of {it:var3} in each cell of the
two-way table.

{title:Example}

{phang2}{cmd:. anova score person drug, repeated(drug)}{p_end}
{phang2}{cmd:. anova output machine / operator|machine /}{p_end}

{title:Stored results}

{pstd}
{cmd:oneway}: {cmd:r(N)}, {cmd:r(F)}, {cmd:r(df_m)}, {cmd:r(df_r)},
{cmd:r(mss)}, {cmd:r(rss)}, {cmd:r(chi2bart)}, {cmd:r(df_bart)}.{break}
{cmd:anova}: {cmd:e(N)}, {cmd:e(mss)}, {cmd:e(rss)}, {cmd:e(df_m)},
{cmd:e(df_r)}, {cmd:e(F)}, {cmd:e(r2)}, {cmd:e(r2_a)}, {cmd:e(rmse)} and, for
each term {it:k}, {cmd:e(term_}{it:k}{cmd:)}, {cmd:e(ss_}{it:k}{cmd:)},
{cmd:e(df_}{it:k}{cmd:)}, {cmd:e(F_}{it:k}{cmd:)}.

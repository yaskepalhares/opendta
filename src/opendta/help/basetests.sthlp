{smcl}
{title:Title}

{p 4 4 2}
{bf:bitest}, {bf:bitesti}, {bf:ksmirnov}, {bf:median}, {bf:nptrend}, {bf:runtest}, {bf:tetrachoric}, {bf:ktau}, {bf:pcorr} {hline 2} Hypothesis tests and correlations

{title:Syntax}

{p 8 16 2}
{cmd:bitest} {it:varname} {cmd:==} {it:#p} [{it:if}] [{it:in}]    {cmd:bitesti} {it:#N} {it:#succ} {it:#p}

{p 8 16 2}
{cmd:ksmirnov} {it:varname} {cmd:,} {cmd:by(}{it:groupvar}{cmd:)} [{cmd:exact}]    {cmd:ksmirnov} {it:varname} {cmd:=} {it:exp}

{p 8 16 2}
{cmd:median} {it:varname} {cmd:,} {cmd:by(}{it:groupvar}{cmd:)} [{cmd:exact} {cmd:medianties(}{it:above|below|drop|split}{cmd:)}]

{p 8 16 2}
{cmd:nptrend} {it:varname} {cmd:,} {cmd:by(}{it:groupvar}{cmd:)} [{cmd:nodetail} {cmd:score(}{it:varname}{cmd:)}]

{p 8 16 2}
{cmd:runtest} {it:varname} [{cmd:,} {cmd:threshold(}#{cmd:)} {cmd:continuity} {cmd:drop} {cmd:split}]

{p 8 16 2}
{cmd:tetrachoric} {it:varlist} [{cmd:,} {cmd:stats(}...{cmd:)} {cmd:print} {cmd:pw}]   {cmd:ktau} {it:varlist} [{cmd:,} {cmd:stats(}...{cmd:)} {cmd:print} {cmd:pw}]

{p 8 16 2}
{cmd:pcorr} {it:varname1 varlist} [{it:weight}]

{title:Description}

{pstd}
{cmd:bitest} is the exact binomial test; the two-sided p-value adds the probabilities of all outcomes not more likely than the observed one. {cmd:ksmirnov} is the one- and two-sample Kolmogorov-Smirnov test. {cmd:median} is the k-sample median test. {cmd:nptrend} is Cuzick's trend test. {cmd:runtest} is the Wald-Wolfowitz runs test. {cmd:tetrachoric} estimates the correlation of two binary variables by maximum likelihood. {cmd:ktau} gives Kendall's tau-a, tau-b and the score, with the tie-corrected standard error. {cmd:pcorr} gives partial and semipartial correlations.

{title:Limitations}

{pstd}
{cmd:tetrachoric} does not implement {cmd:zeroadjust} and {cmd:ktau} with several variables shows tau-b only. The layouts the manual does not show are marked VERIFY in the code.

{smcl}
{title:Title}

{p 4 4 2}
{bf:exlogistic}, {bf:expoisson} {hline 2} Exact logistic and exact Poisson regression

{title:Syntax}

{p 8 16 2}
{cmd:exlogistic} depvar indepvars [if] [in] [{cmd:,} {cmd:test(}{cmd:sufficient}|{cmd:score}|{cmd:probability}{cmd:)} {cmd:level(}#{cmd:)}]

{p 8 16 2}
{cmd:expoisson} depvar indepvars [if] [in] [{cmd:,} {cmd:test(}{cmd:sufficient}|{cmd:score}|{cmd:probability}{cmd:)} {cmd:level(}#{cmd:)}]

{title:Description}

{pstd}
Estimates, confidence intervals and p-values come from the exact conditional distribution of the sufficient statistics. The conditional maximum likelihood estimate is shown, or the median unbiased estimate (marked with an asterisk) when the first is infinite.

{pstd}
Not available: {cmd:terms}, {cmd:group()}, {cmd:binomial()}, {cmd:condvars()}, {cmd:saving()}, {cmd:exposure()} and {cmd:offset()}.

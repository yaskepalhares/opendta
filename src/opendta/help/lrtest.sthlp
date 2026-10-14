{smcl}
{title:Title}

{p 4 4 2}
{bf:lrtest and hausman} {hline 2} Likelihood-ratio and Hausman specification tests

{title:Syntax}

{p 8 16 2}
{cmd:lrtest} modelspec1 [modelspec2] [, df(#) force]

{p 8 16 2}
{cmd:hausman} consistent [efficient] [, constant sigmamore sigmaless df(#)]

{title:Description}

{pstd}
LR test of nested models stored with {cmd:estimates store} ({cmd:.} is the current model).

{pstd}
Hausman test of the difference between two stored estimates.

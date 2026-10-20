{smcl}
{title:Title}

{p 4 4 2}
{bf:mkspline}, {bf:xi}, {bf:nestreg}, {bf:stepwise}, {bf:fp}, {bf:mfp} {hline 2} Building and selecting model terms

{title:Syntax}

{p 8 16 2}
{cmd:mkspline} {it:stubname}# {cmd:=} {it:varname} [{cmd:,} {cmd:cubic} {cmd:knots(}{it:numlist}{cmd:)} {cmd:nknots(}#{cmd:)} {cmd:marginal} {cmd:displayknots}]

{p 8 16 2}
{cmd:xi} [{cmd:,} {cmd:prefix(}{it:string}{cmd:)} {cmd:noomit}] {cmd::} {it:command} {it:term} ...

{p 8 16 2}
{cmd:nestreg} [{cmd:,} {cmd:lr}] {cmd::} {it:command depvar} ({it:varlist}) ({it:varlist}) ...

{p 8 16 2}
{cmd:stepwise} [{cmd:,} {cmd:pr(}#{cmd:)} {cmd:pe(}#{cmd:)} {cmd:forward} {cmd:hierarchical} {cmd:lockterm1}] {cmd::} {it:command depvar varlist}

{p 8 16 2}
{cmd:fp} [{cmd:<}{it:xvar}{cmd:>}] [{cmd:,} {cmd:powers(}...{cmd:)} {cmd:df(}#{cmd:)}] {cmd::} {it:command}   {cmd:mfp} [{cmd:,} {cmd:select(}#{cmd:)} {cmd:alpha(}#{cmd:)}] {cmd::} {it:command depvar varlist}

{title:Description}

{pstd}
{cmd:mkspline} creates linear splines (equally spaced, percentile or given knots) and restricted cubic splines. {cmd:xi} expands {cmd:i.}{it:varname} and {cmd:I.}{it:varname} terms and their interactions into dummy variables. {cmd:nestreg} fits nested models and tests each block with a Wald or LR test. {cmd:stepwise} does backward or forward selection. {cmd:fp} compares fractional polynomials and {cmd:mfp} runs the closed test procedure.

{title:Limitations}

{pstd}
Scaling of {cmd:fp} and the layouts of {cmd:fp} and {cmd:mfp} are marked VERIFY in the code.

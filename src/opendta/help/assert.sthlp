{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:assert} {hline 2}}Check that a condition is true{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:assert} {it:exp} [{it:if}] [{it:in}] [{cmd:,} {cmd:rc0} {cmd:null}]

{title:Description}

{pstd}
{cmd:assert} stops with error {cmd:r(9)} when {it:exp} is false for any
observation, and reports how many observations contradict it. It is the
simplest way to make a do-file check its own assumptions.
{cmd:rc0} reports without stopping; {cmd:null} accepts no observations.
{p_end}

{title:Example}

        {cmd:. assert idade >= 0 & idade < 120 if !missing(idade)}

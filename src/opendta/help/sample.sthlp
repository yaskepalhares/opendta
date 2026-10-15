{smcl}
{title:Title}

{p2colset 5 19 21 2}{...}
{p2col :{cmd:sample} {hline 2}}Draw a random sample{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:sample} {it:#} [{it:if}] [{it:in}] [{cmd:,} {opt count} {opt by(groupvars)}]{p_end}

{p 4 4 2}{cmd:by} is allowed.{p_end}

{title:Description}

{pstd}
{cmd:sample} draws a simple random sample, without replacement, of the
data in memory and drops the other observations. {it:#} is the percentage
of observations to keep (0 to 100); with {opt count} it is the number of
observations. Observations that do not meet {it:if} and {it:in} are kept.
{p_end}

{pstd}
The draw uses the random-number generator: set the seed first
({cmd:set seed}) to get the same sample again. The data keep their order.
{p_end}

{title:Options}

{phang}{opt count} takes {it:#} as a number of observations instead of a
percentage. A group with fewer observations is kept whole.{p_end}

{phang}{opt by(groupvars)} samples each group separately: {it:#} percent
(or {it:#} observations with {opt count}) of every group.{p_end}

{title:Examples}

{phang}{cmd:. set seed 2026}{p_end}
{phang}{cmd:. sample 10}{p_end}
{phang}{cmd:. sample 500, count}{p_end}
{phang}{cmd:. sample 20 if year == 2024}{p_end}
{phang}{cmd:. sample 50, count by(region)}{p_end}

{title:Also see}

{pstd}
{help drop}, {help keep}, {help bootstrap}
{p_end}

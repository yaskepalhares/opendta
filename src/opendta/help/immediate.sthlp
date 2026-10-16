{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:cii} {hline 2}}Immediate commands (statistics from typed numbers){p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:cii} {it:#obs} {it:#mean} {it:#sd} [{cmd:,} {opt l:evel(#)}]{p_end}
{p 8 16 2}{cmd:cii} {it:#obs} {it:#successes} [{cmd:,} {opt l:evel(#)} {opt exact} {opt wilson} {opt agresti} {opt jeffreys} {opt wald}]{p_end}
{p 8 16 2}{cmd:cii} {it:#exposure} {it:#events} {cmd:,} {opt poisson}{p_end}
{p 8 16 2}{cmd:cii variances} {it:#obs} {it:#variance} [{cmd:,} {opt l:evel(#)} {opt sd}]  ({opt sd}: the second number is a standard deviation){p_end}
{p 8 16 2}{cmd:tabi} {it:#11} {it:#12} [...] {cmd:\} {it:#21} {it:#22} [...] [{cmd:,} {it:tabulate options}]{p_end}
{p 8 16 2}{cmd:ttesti} {it:#obs} {it:#mean} {it:#sd} {it:#val}  |  {cmd:ttesti} {it:#obs1} {it:#mean1} {it:#sd1} {it:#obs2} {it:#mean2} {it:#sd2}{p_end}
{p 8 16 2}{cmd:prtesti} {it:#obs} {it:#p1} {it:#p2}  |  {cmd:prtesti} {it:#obs1} {it:#p1} {it:#obs2} {it:#p2}{p_end}

{title:Description}

{pstd}
Immediate commands take the summary numbers as arguments instead of variables. OpenDTA follows the syntax of Stata 14.0 for {cmd:cii} (the {cmd:cii means} and {cmd:cii proportions} forms came with 14.1 and are also accepted).
{p_end}

{title:Examples}

{phang2}{cmd:. cii 50 12}{p_end}
{phang2}{cmd:. cii 100 3.2 1.1}{p_end}
{phang2}{cmd:. tabi 10 20 \ 30 40, chi2 exact}{p_end}
{phang2}{cmd:. ttesti 24 62.6 15.8 75}{p_end}
{phang2}{cmd:. prtesti 30 .4 45 .6}{p_end}
{phang2}{cmd:. cii variances 40 2.25}{p_end}

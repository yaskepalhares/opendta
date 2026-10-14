{smcl}
{title:Title}

{p2colset 5 26 28 2}{...}
{p2col :{cmd:graph bar} {hline 2}}Bar, box, pie and matrix graphs{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:graph} {c -(}{cmd:bar}|{cmd:hbar}{c )-} [{cmd:(}{it:stat}{cmd:)}] {it:varlist} [{cmd:(}{it:stat}{cmd:)} {it:varlist} ...] [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:graph} {c -(}{cmd:bar}|{cmd:hbar}{c )-} {cmd:,} {cmd:over(}{it:var}{cmd:)} [{it:options}]{p_end}
{p 8 16 2}{cmd:graph} {c -(}{cmd:box}|{cmd:hbox}{c )-} {it:varlist} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:graph pie} {it:varlist} | {it:varname}{cmd:, over(}{it:var}{cmd:)} | {cmd:, over(}{it:var}{cmd:)} [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:graph matrix} {it:varlist} [{it:if}] [{it:in}] [{cmd:,} {cmd:half}]{p_end}

{pstd}
{it:stat} is {cmd:mean} (default), {cmd:median}, {cmd:p}{it:#}, {cmd:sum}, {cmd:rawsum},
{cmd:count}, {cmd:percent}, {cmd:min}, {cmd:max}, {cmd:sd}, {cmd:semean}, {cmd:iqr},
{cmd:first}, {cmd:last}. Without a {it:varlist}, bars show the percentage
of observations in each {cmd:over()} group.
{p_end}

{title:Options}

{p2colset 9 30 32 2}{...}
{p2col :{cmd:over(}{it:var}{cmd:)}}one bar or box per group; may be repeated{p_end}
{p2col :{cmd:asyvars}}treat the first {cmd:over()} groups as series (colors and legend){p_end}
{p2col :{cmd:stack}}stack the series{p_end}
{p2col :{cmd:percentages}}show each series as a percentage of the category total{p_end}
{p2col :{cmd:blabel(bar)}}write the values on the bars{p_end}
{p2col :{cmd:nooutsides}}box: hide points beyond the whiskers{p_end}
{p2col :{cmd:plabel(_all percent|sum)}}pie: labels on the slices{p_end}
{p2col :{cmd:half}}matrix: lower triangle only{p_end}
{p2colreset}{...}

{pstd}
Boxes span the 25th to 75th percentiles (computed as in {help summarize});
whiskers reach the most extreme values within 1.5 IQR of the box; values
beyond them are drawn as points.
{p_end}

{title:Examples}

{phang2}{cmd:. graph bar (mean) wage, over(occupation) blabel(bar)}{p_end}
{phang2}{cmd:. graph hbar (count), over(region)}{p_end}
{phang2}{cmd:. graph box income, over(sex)}{p_end}
{phang2}{cmd:. graph pie, over(region) plabel(_all percent)}{p_end}
{phang2}{cmd:. graph matrix price mpg weight, half}{p_end}

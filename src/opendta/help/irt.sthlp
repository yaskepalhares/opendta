{smcl}
{title:Title}

{p2colset 5 16 18 2}{...}
{p2col :{cmd:irt} {hline 2}}Item response theory models{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:irt} {cmd:1pl}|{cmd:2pl}|{cmd:3pl}|{cmd:grm}|{cmd:nrm}|{cmd:pcm}|{cmd:gpcm}|{cmd:rsm} {it:varlist} [{it:if}] [{it:in}]
[{cmd:,} {opt list:wise} {opt intm:ethod()} {opt intp:oints(#)} {opt l:evel(#)} {cmd:nolog}]{p_end}

{p 8 16 2}{cmd:irt hybrid} {cmd:(}{it:model} {it:varlist}{cmd:)} {cmd:(}{it:model} {it:varlist}{cmd:)} ... [{cmd:,} {it:options}]{p_end}

{p 8 16 2}{cmd:estat report} [{it:varlist}] [{cmd:,} {opt byparm} {opt sort(a|b)}]{p_end}

{p 8 16 2}{cmd:irtgraph} {cmd:icc}|{cmd:iif}|{cmd:tcc}|{cmd:tif} [{it:varlist}] [{cmd:,} {opt range(# #)} {it:graph_options}]{p_end}

{title:Description}

{pstd}
{cmd:irt} fits models for binary items (one-, two- and three-parameter logistic), ordered items
(graded response, partial credit, generalized partial credit, rating scale) and nominal items,
with one latent trait distributed N(0, 1). The fit is done by {help gsem}; the table reports
the discrimination and the difficulties of each item (by the delta method). {cmd:irt hybrid}
uses a different model for each group of items.
{p_end}

{pstd}
{cmd:estat report} shows the parameters again, grouped by item or by parameter
({cmd:byparm}). {cmd:irtgraph} draws the item characteristic curves ({cmd:icc}), the item and
test information functions ({cmd:iif}, {cmd:tif}) and the test characteristic curve ({cmd:tcc}).
{p_end}

{title:Examples}

{phang2}{cmd:. irt 2pl q1-q8}{p_end}
{phang2}{cmd:. estat report, byparm}{p_end}
{phang2}{cmd:. irtgraph icc q1 q2}{p_end}
{phang2}{cmd:. irt grm r1-r5}{p_end}

{title:OpenDTA notes}

{pstd}
The names of the rows of the nominal and rating-scale models follow the manual and may differ
from Stata. {cmd:group()} (multiple-group IRT) is not available.
{p_end}

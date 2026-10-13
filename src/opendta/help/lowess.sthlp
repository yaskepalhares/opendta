{smcl}
{title:Title}

{p 4 4 2}
{bf:lowess} {hline 2} Lowess smoothing

{title:Syntax}

{p 8 16 2}
{cmd:lowess} yvar xvar [if] [in] [, bwidth(#) mean noweight logit adjust generate(newvar) nograph twoway_options]

{title:Description}

{pstd}
Locally weighted regression of {it:yvar} on {it:xvar} using a fraction {cmd:bwidth()} of the data (default .8) with tricube weights. Also available as the {cmd:twoway lowess} plot type.

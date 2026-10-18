{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:symmetry} {hline 2}}Symmetry and marginal homogeneity tests{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:symmetry} {it:casevar} {it:controlvar} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt notable} {opt contrib} {opt ex:act} {opt mh} {opt cc}]{p_end}

{title:Description}

{pstd}
{cmd:symmetry} tabulates two paired categorical variables (for example, the exposure of matched cases and controls) and tests the symmetry of the table (Bowker) and the homogeneity of its margins (Stuart-Maxwell; with {opt mh}, also Bickenboller's statistic). {opt contrib} shows each pair of cells' contribution to the symmetry statistic, {opt exact} its exact significance and {opt cc} a continuity correction.
{p_end}

{title:Examples}

{phang2}{cmd:. symmetry before after}{p_end}
{phang2}{cmd:. symmetry before after [fw=n], contrib}{p_end}

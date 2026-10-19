{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:ds} {hline 2}}Listing and inspecting variables and value labels{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:ds} [{it:varlist}] [{cmd:,} {opt not} {opt a:lpha} {opt d:etail} {opt varw:idth(#)} {opt sk:ip(#)} {opt has(spec)} {opt not(spec)} {opt ins:ensitive} {opt ind:ent(#)}]{p_end}
{p 8 16 2}{cmd:lookfor} {it:string} [{it:string} ...]{p_end}
{p 8 16 2}{cmd:inspect} [{it:varlist}] [{it:if}] [{it:in}]{p_end}
{p 8 16 2}{cmd:labelbook} [{it:lblname-list}] [{cmd:,} {opt a:lpha} {opt len:gth(#)} {opt li:st(#)} {opt p:roblems} {opt d:etail}]{p_end}

{title:Description}

{pstd}
{cmd:ds} lists variable names compactly and leaves them in {cmd:r(varlist)}; {opt has()} and {opt not()} select by {cmd:type}, {cmd:varlabel}, {cmd:vallabel}, {cmd:format} or {cmd:char}. {cmd:lookfor} finds the variables whose name or label contains any of the words. {cmd:inspect} shows, per variable, the counts of negative, zero, positive, integer, noninteger and missing values, a small histogram and the number of unique values. {cmd:labelbook} describes value labels: their values, lengths, gaps, uniqueness and the variables that use them.
{p_end}

{title:Examples}

{phang2}{cmd:. ds, has(type string)}{p_end}
{phang2}{cmd:. lookfor income}{p_end}
{phang2}{cmd:. inspect age}{p_end}
{phang2}{cmd:. labelbook, problems}{p_end}

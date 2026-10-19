{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:outfile} {hline 2}}Writing data as text{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:outfile} [{it:varlist}] {cmd:using} {it:filename} [{it:if}] [{it:in}] [{cmd:,} {opt dict:ionary} {opt nolabel} {opt noquote} {opt c:omma} {opt w:ide} {opt rjs} {opt fjs} {opt mi:ssing} {opt run:together} {opt replace}]{p_end}

{title:Description}

{pstd}
{cmd:outfile} writes the data to a text file in the layout that {cmd:infile} reads, using the display formats; with {opt dictionary} it writes an {cmd:infile} dictionary at the top. Value labels are written instead of the codes unless {opt nolabel}.
{p_end}

{title:Examples}

{phang2}{cmd:. outfile using mydata.raw, replace}{p_end}
{phang2}{cmd:. outfile id age using mydata.raw, dictionary replace}{p_end}

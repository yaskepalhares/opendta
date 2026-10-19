{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:postfile} {hline 2}}Save results in a Stata dataset{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:postfile} {it:postname} {it:newvarlist} {cmd:using} {it:filename} [{cmd:,} {opt every(#)}
{opt replace}]

{p 8 16 2}
{cmd:post} {it:postname} {cmd:(}{it:exp}{cmd:)} {cmd:(}{it:exp}{cmd:)} ... {cmd:(}{it:exp}{cmd:)}

{p 8 16 2}
{cmd:postclose} {it:postname}

{p 8 16 2}
{cmd:postutil} {c -(}{cmd:dir} | {cmd:clear}{c )-}

{title:Description}

{pstd}
{cmd:postfile} declares the variables of a new dataset, which {cmd:post} fills one
observation at a time and {cmd:postclose} completes. The data in memory are not touched.
{it:newvarlist} may give types: {cmd:a b str20 c double(d e f)}; the default type is the one
of {cmd:set type}. A string longer than the type is cut. Without extension, the file is
{it:filename}{cmd:.dta}; an existing file is replaced only with {cmd:replace}.
{p_end}

{pstd}
Name {it:postname} with {cmd:tempname}: it is then forgotten when the program ends, even
after an error. {cmd:every()} is accepted and has no effect, since the rows are written when
the file is closed. {cmd:postutil dir} lists the open postfiles and {cmd:postutil clear}
closes them all.
{p_end}

{title:Example}

{phang2}{cmd:. tempname sim}{p_end}
{phang2}{cmd:. postfile `sim' mean var using results, replace}{p_end}
{phang2}{cmd:. post `sim' (r(mean)) (r(Var))}{p_end}
{phang2}{cmd:. postclose `sim'}{p_end}

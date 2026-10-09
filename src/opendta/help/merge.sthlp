{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:merge} {hline 2}}Combine datasets{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:merge} {c -(}{cmd:1:1}|{cmd:m:1}|{cmd:1:m}|{cmd:m:m}{c )-} {it:varlist} {cmd:using} {it:filename} [{cmd:,} {it:options}]

{p 8 16 2}{cmd:append using} {it:filename} [{it:filename} ...] [{cmd:,} {opt gen:erate(newvar)} {opt keep(varlist)} {opt nol:abel} {opt force}]{p_end}
{p 8 16 2}{cmd:joinby} [{it:varlist}] {cmd:using} {it:filename} [{cmd:,} {opt un:matched(none|both|master|using)} {opt update} {opt replace}]{p_end}
{p 8 16 2}{cmd:cross using} {it:filename}{p_end}

{title:Options of merge}

{p2colset 9 30 32 2}{...}
{p2col :{opt keepus:ing(varlist)}}variables to take from the using data{p_end}
{p2col :{opt gen:erate(name)}, {opt nogen:erate}}name of the result variable, or none{p_end}
{p2col :{opt upd:ate}, {opt rep:lace}}fill missing values, or overwrite them{p_end}
{p2col :{opt keep(results)}, {opt assert(results)}}{cmd:master}, {cmd:using}, {cmd:match}, {cmd:match_update}, {cmd:match_conflict} or 1-5{p_end}
{p2col :{opt nol:abel}, {opt norep:ort}, {opt force}}{p_end}
{p2colreset}{...}

{title:Description}

{pstd}
{cmd:merge} matches observations by key variables (or by {cmd:_n}) and
reports how many came from each side in {cmd:_merge}: 1 master only,
2 using only, 3 matched, 4 missing updated, 5 nonmissing conflict.
{cmd:append} adds observations, {cmd:joinby} forms all pairs within
groups and {cmd:cross} forms all pairs.
{p_end}

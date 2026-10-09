{smcl}
{title:Title}

{p2colset 5 17 19 2}{...}
{p2col :{cmd:egen} {hline 2}}Extensions to generate{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
[{cmd:by} {it:varlist}{cmd::}] {cmd:egen} [{it:type}] {it:newvar} {cmd:=} {it:fcn}{cmd:(}{it:arguments}{cmd:)} [{it:if}] [{it:in}] [{cmd:,} {it:options}]

{title:Functions}

{pstd}Computed within groups of {cmd:by} or {opt by()}:{p_end}
{p2colset 9 30 32 2}{...}
{p2col :{cmd:count(}{it:exp}{cmd:)}}nonmissing values{p_end}
{p2col :{cmd:mean}, {cmd:sd}, {cmd:total}, {cmd:min}, {cmd:max}}{opt missing} keeps {cmd:total} missing when all are missing{p_end}
{p2col :{cmd:median}, {cmd:pctile}, {cmd:iqr}}{opt p(#)} sets the percentile{p_end}
{p2col :{cmd:mode}}{opt minmode}, {opt maxmode} break ties{p_end}
{p2col :{cmd:skew}, {cmd:kurt}, {cmd:mad}, {cmd:mdev}}shape and dispersion{p_end}
{p2col :{cmd:std(}{it:exp}{cmd:)}}standardized values; {opt mean(#)} {opt sd(#)}{p_end}
{p2col :{cmd:rank(}{it:exp}{cmd:)}}{opt field}, {opt track}, {opt unique}{p_end}
{p2col :{cmd:seq()}}sequence; {opt from(#)} {opt to(#)} {opt block(#)}{p_end}
{p2colreset}{...}

{pstd}Other functions:{p_end}
{p2colset 9 30 32 2}{...}
{p2col :{cmd:group(}{it:varlist}{cmd:)}}one number per combination; {opt label} {opt missing}{p_end}
{p2col :{cmd:tag(}{it:varlist}{cmd:)}}1 in the first observation of each combination{p_end}
{p2col :{cmd:fill(}{it:numlist}{cmd:)}}continues a sequence or a repeating pattern{p_end}
{p2col :{cmd:cut(}{it:exp}{cmd:)}}{opt at(numlist)} or {opt group(#)}; {opt icodes}{p_end}
{p2col :{cmd:concat(}{it:varlist}{cmd:)}}joins values as text; {opt punct()} {opt decode}{p_end}
{p2col :{cmd:ends(}{it:strvar}{cmd:)}}{opt head}, {opt last}, {opt tail}; {opt punct()}{p_end}
{p2col :{cmd:diff(}{it:varlist}{cmd:)}}1 when the variables differ{p_end}
{p2colreset}{...}

{pstd}Across variables in each observation:{p_end}
{p2colset 9 30 32 2}{...}
{p2col :{cmd:rowmean}, {cmd:rowtotal}, {cmd:rowsd}}{p_end}
{p2col :{cmd:rowmin}, {cmd:rowmax}, {cmd:rowmedian}, {cmd:rowpctile}}{p_end}
{p2col :{cmd:rownonmiss}, {cmd:rowmiss}, {cmd:rowfirst}, {cmd:rowlast}}{p_end}
{p2col :{cmd:anycount}, {cmd:anymatch}, {cmd:anyvalue}}{opt values(numlist)}{p_end}
{p2colreset}{...}

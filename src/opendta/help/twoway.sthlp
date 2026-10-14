{smcl}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:twoway} {hline 2}}Two-way graphs: scatter, line, fits and more{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}[{cmd:graph}] {cmd:twoway} {it:plot} [{it:if}] [{it:in}] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}[{cmd:graph}] {cmd:twoway} {cmd:(}{it:plot}{cmd:)} {cmd:(}{it:plot}{cmd:)} ... [{cmd:,} {it:options}]{p_end}
{p 8 16 2}[{cmd:graph}] {cmd:twoway} {it:plot} {cmd:||} {it:plot} {cmd:||} ... [{cmd:||,} {it:options}]{p_end}

{pstd}
{cmd:scatter}, {cmd:line}, {cmd:connected} and {cmd:tsline} may be typed without {cmd:twoway}.
{p_end}

{p2colset 9 30 32 2}{...}
{p2col :{it:plot}}draws{p_end}
{p2line}
{p2col :{cmd:scatter} {it:y}... {it:x}}markers; {cmd:[aweight]} scales marker size{p_end}
{p2col :{cmd:line} {it:y}... {it:x}}lines in data order (add {cmd:sort}){p_end}
{p2col :{cmd:connected} {it:y}... {it:x}}lines with markers{p_end}
{p2col :{cmd:tsline} {it:y}...}lines against the {help tsset} time variable{p_end}
{p2col :{cmd:lfit}, {cmd:qfit} {it:y x}}linear or quadratic fit{p_end}
{p2col :{cmd:lfitci}, {cmd:qfitci} {it:y x}}fit with confidence band ({cmd:level()}){p_end}
{p2col :{cmd:area}, {cmd:bar}, {cmd:spike}, {cmd:dropline} {it:y x}}filled areas, bars, spikes{p_end}
{p2col :{cmd:rcap}, {cmd:rarea}, {cmd:rbar}, {cmd:rspike}, {cmd:rline}, {cmd:rconnected}, {cmd:rscatter}, {cmd:rcapsym} {it:y1 y2 x}}ranges{p_end}
{p2col :{cmd:function} [{it:y} {cmd:=}] {it:exp}}a function of {cmd:x}; {cmd:range(}{it:# #}{cmd:)} (default 0 1), {cmd:n(#)} (default 300){p_end}
{p2col :{cmd:histogram} {it:x}}see {help histogram}{p_end}
{p2col :{cmd:kdensity} {it:x}}see {help histogram}{p_end}
{p2col :{cmd:scatteri} {it:# #} ...}markers at the given y x pairs{p_end}
{p2line}
{p2colreset}{...}

{title:Plot options}

{pstd}
{cmd:mcolor()}, {cmd:msymbol()} ({cmd:O}, {cmd:D}, {cmd:S}, {cmd:T}, {cmd:X}, {cmd:+}, {cmd:p},
hollow {cmd:Oh}, {cmd:Dh}, {cmd:Sh}, {cmd:Th}), {cmd:msize()}, {cmd:mlabel(}{it:var}{cmd:)},
{cmd:lcolor()}, {cmd:lpattern()} ({cmd:solid}, {cmd:dash}, {cmd:dot}, {cmd:dash_dot},
{cmd:shortdash}, {cmd:longdash}), {cmd:lwidth()}, {cmd:color()}, {cmd:fcolor()},
{cmd:barwidth()}, {cmd:base()}, {cmd:horizontal}, {cmd:sort}. Colors are names
({cmd:navy}, {cmd:maroon}, {cmd:gs8}...) or "{it:R G B}".
{p_end}

{pstd}
Graph options ({cmd:title()}, axes, {cmd:legend()}, {cmd:by()}, {cmd:name()}...) are described in
{help graph}; they may appear in the options of any plot.
{p_end}

{title:Defaults}

{pstd}
With one plotted series, the axis titles are the variable labels (or names)
and there is no legend. With two or more, the y title is left blank and a
legend identifies the series.
{p_end}

{title:Examples}

{phang2}{cmd:. scatter price weight}{p_end}
{phang2}{cmd:. twoway (scatter price weight) (lfitci price weight), by(foreign)}{p_end}
{phang2}{cmd:. line le_male le_female year, sort legend(order(1 "Men" 2 "Women"))}{p_end}
{phang2}{cmd:. twoway function y = normalden(x), range(-4 4)}{p_end}

{smcl}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:graph} {hline 2}}Draw and manage graphs{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:graph} {it:family} ...{p_end}

{p2colset 9 32 34 2}{...}
{p2col :{it:family}}see{p_end}
{p2line}
{p2col :{cmd:twoway}}{help twoway}: scatter, line, fits, areas, bars, ranges, functions{p_end}
{p2col :{cmd:bar}, {cmd:hbar}, {cmd:dot}}{help graph_bar}{p_end}
{p2col :{cmd:box}, {cmd:hbox}}{help graph_bar}{p_end}
{p2col :{cmd:pie}}{help graph_bar}{p_end}
{p2col :{cmd:matrix}}{help graph_bar}{p_end}
{p2col :{cmd:combine}}below{p_end}
{p2line}
{p2colreset}{...}

{p 8 16 2}{cmd:histogram} and {cmd:kdensity}: see {help histogram}.{p_end}

{title:Managing graphs}

{p 8 16 2}{cmd:graph dir} [{it:pattern}] [{cmd:,} {cmd:memory}]{p_end}
{p 8 16 2}{cmd:graph drop} {it:name} [{it:name} ...] | {cmd:_all}{p_end}
{p 8 16 2}{cmd:graph rename} [{it:oldname}] {it:newname} [{cmd:,} {cmd:replace}]{p_end}
{p 8 16 2}{cmd:graph copy} [{it:oldname}] {it:newname} [{cmd:,} {cmd:replace}]{p_end}
{p 8 16 2}{cmd:graph display} [{it:name}] [{cmd:,} {cmd:scheme(}{it:name}{cmd:)} {cmd:xsize(#)} {cmd:ysize(#)}]{p_end}
{p 8 16 2}{cmd:graph export} {it:filename}{cmd:.}{it:suffix} [{cmd:,} {cmd:as(}{it:suffix}{cmd:)} {cmd:name(}{it:name}{cmd:)} {cmd:replace} {cmd:width(#)} {cmd:height(#)}]{p_end}
{p 8 16 2}{cmd:graph save} [{it:name}] {it:filename} [{cmd:,} {cmd:replace}]{p_end}
{p 8 16 2}{cmd:graph use} {it:filename} [{cmd:,} {cmd:name(}{it:name}[{cmd:, replace}]{cmd:)} {cmd:nodraw}]{p_end}
{p 8 16 2}{cmd:graph combine} {it:name} [{it:name} ...] [{cmd:,} {cmd:cols(#)} {cmd:rows(#)} {it:title_options}]{p_end}
{p 8 16 2}{cmd:graph query, schemes}{p_end}

{title:Description}

{pstd}
Every graph is kept in memory under a name: {cmd:Graph} unless option
{cmd:name()} gives another. Drawing a new graph with a name that is already
used requires {cmd:name(}{it:name}{cmd:, replace)}; the default name
{cmd:Graph} is always replaced. In the interface, graphs appear in the Graph
window, one tab per name. In batch mode they are only kept in memory, and can
be written with {cmd:graph export} or {cmd:graph save}.
{p_end}

{pstd}
{cmd:graph export} writes the image. The format comes from the suffix:
{cmd:.png}, {cmd:.pdf}, {cmd:.svg}, {cmd:.eps}, {cmd:.ps}, {cmd:.tif} or
{cmd:.jpg}. {cmd:width()} and {cmd:height()} set the size in pixels of
bitmap formats.
{p_end}

{pstd}
{cmd:graph save} writes a {cmd:.gph} file. OpenDTA's {cmd:.gph} files are
its own format (a JSON description of the graph, not an image); Stata's
{cmd:.gph} files cannot be read, and Stata cannot read OpenDTA's.
{p_end}

{title:Common options}

{p2colset 5 30 32 2}{...}
{p2col :{cmd:title()}, {cmd:subtitle()}}text above the graph; several quoted strings make several lines{p_end}
{p2col :{cmd:note()}, {cmd:caption()}}text below the graph{p_end}
{p2col :{cmd:xtitle()}, {cmd:ytitle()}}axis titles{p_end}
{p2col :{cmd:xlabel()}, {cmd:ylabel()}}ticks: a {it:numlist}, value-"label" pairs, {cmd:#}{it:k}, {cmd:none}; suboptions {cmd:angle()}, {cmd:grid}, {cmd:nogrid}, {cmd:valuelabel}{p_end}
{p2col :{cmd:xscale()}, {cmd:yscale()}}{cmd:log}, {cmd:range()}, {cmd:reverse}, {cmd:off}{p_end}
{p2col :{cmd:xline()}, {cmd:yline()}}reference lines{p_end}
{p2col :{cmd:legend()}}{cmd:off}, {cmd:on}, {cmd:order()}, {cmd:label()}, {cmd:rows()}, {cmd:cols()}{p_end}
{p2col :{cmd:by()}}one panel per group; suboptions {cmd:total}, {cmd:cols()}, {cmd:rows()}, {cmd:note()}{p_end}
{p2col :{cmd:name()}}name in memory{p_end}
{p2col :{cmd:saving()}}also save to a {cmd:.gph} file{p_end}
{p2col :{cmd:scheme()}}visual scheme{p_end}
{p2col :{cmd:nodraw}}keep the graph in memory without showing it{p_end}
{p2colreset}{...}

{title:Schemes}

{pstd}
OpenDTA has its own schemes: {cmd:opendta} (color), {cmd:opendta_mono}
(grayscale) and {cmd:opendta_dark}. Stata scheme names such as
{cmd:s2color}, {cmd:s2mono} or {cmd:s1mono} are accepted and mapped to the
closest OpenDTA scheme, so do-files run unchanged; the look is OpenDTA's, not
Stata's. {cmd:set scheme} {it:name} changes the default; {cmd:c(scheme)} returns it.
{cmd:set graphics off} keeps graphs in memory without opening the Graph window.
{p_end}

{title:Examples}

{phang2}{cmd:. scatter mpg weight, name(a)}{p_end}
{phang2}{cmd:. histogram mpg, name(b)}{p_end}
{phang2}{cmd:. graph combine a b, cols(2) title("Two graphs")}{p_end}
{phang2}{cmd:. graph export combined.png, width(1600)}{p_end}
{phang2}{cmd:. graph dir, memory}{p_end}
{phang2}{cmd:. graph drop _all}{p_end}

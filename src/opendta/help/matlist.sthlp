{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:matlist} {hline 2}}Display a matrix and control its format{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:matlist} {it:matrix_exp} [{cmd:,} {it:style_options} {it:general_options}]

{p 8 16 2}
{cmd:matlist} {it:matrix_exp}{cmd:,} {opt cspec(cspec)} {opt rspec(rspec)} [{it:general_options}]

{synoptset 24}{...}
{synopthdr:style_options}
{synoptline}
{synopt :{opt lines(lstyle)}}{cmd:oneline} (default), {cmd:eq}, {cmd:rowtotal}, {cmd:coltotal},
{cmd:rctotal}, {cmd:rows}, {cmd:columns}, {cmd:cells} or {cmd:none}{p_end}
{synopt :{opt border(bspec)}}{cmd:none}, {cmd:all}, {cmd:rows}, {cmd:columns}, {cmd:left},
{cmd:right}, {cmd:top}, {cmd:bottom}{p_end}
{synopt :{opt format(%fmt)}}display format; default {cmd:%9.0g}{p_end}
{synopt :{opt twidth(#)}}width of the row names; default 12{p_end}
{synopt :{opt left(#)}}indent of the table{p_end}
{synopt :{opt right(#)}}right margin, which affects where columns wrap{p_end}
{synoptline}

{synopthdr:general_options}
{synoptline}
{synopt :{opt title(string)}}title above the table{p_end}
{synopt :{opt tindent(#)}}indent of the title{p_end}
{synopt :{opt rowtitle(string)}}header of the row names{p_end}
{synopt :{opt names(rows|columns|all|none)}}which names to show; {opt nonames} is {cmd:names(none)}{p_end}
{synopt :{opt aligncolnames(ralign|lalign|center)}}alignment of the column names{p_end}
{synopt :{opt noblank}}no blank line before the table{p_end}
{synopt :{opt nohalf}}the whole of a symmetric matrix, not only its lower triangle{p_end}
{synopt :{opt nodotz}}show {cmd:.z} as blank{p_end}
{synopt :{opt underscore}}show underscores in names as blanks{p_end}
{synopt :{opt linesize(#)}}width of the page{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:matlist} shows a matrix, or a matrix expression, as a table with the names as headers.
Columns that do not fit in the page wrap, repeating the row names. With {cmd:cspec()} and
{cmd:rspec()} each column has its own format: {cmd:cspec()} lists separators ({cmd:&} no line,
{cmd:|} a line, {cmd:o}{it:#} spaces) and fields (qualifiers {cmd:L R C w}{it:#}, then
{cmd:%}{it:#}{cmd:.}{it:#}{cmd:f}, {cmd:g} or {cmd:%}{it:#}{cmd:s}); {cmd:rspec()} has a
character for each horizontal boundary, {cmd:-} or {cmd:|} for a line and {cmd:&} for none.
{p_end}

{pstd}
Matrices in OpenDTA carry no equation names yet, so {cmd:showcoleq()}, {cmd:roweqonly},
{cmd:coleqonly}, {cmd:colorcoleq()} and {cmd:keepcoleq} are not available.
{p_end}

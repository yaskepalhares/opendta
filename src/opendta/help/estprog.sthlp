{smcl}
{title:Title}

{p2colset 5 24 26 2}{...}
{p2col :{cmd:signestimationsample} {hline 2}}Tools for programmers of estimation commands{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:signestimationsample} {it:varlist}{break}
{cmd:checkestimationsample}{break}
{cmd:makecns} {it:matname} [{cmd:,} {opt nocnsnotes} {opt r}]{break}
{cmd:matcproc} {it:T} {it:a} {it:C}

{title:Description}

{pstd}
{cmd:signestimationsample} stores in {cmd:e(datasignaturevars)} and {cmd:e(datasignature)} the
signature of {it:varlist} (and of the weight) over the observations of {cmd:e(sample)};
{cmd:checkestimationsample} fails with {cmd:r(459)} if the data in the estimation sample have
changed since.
{p_end}

{pstd}
{cmd:makecns} {it:matname} stores in {cmd:e(Cns)} the matrix of constraints {it:matname}, which has
the columns of {cmd:e(b)} and one more, so that {it:R} {it:b}' = {it:r} for {it:matname} = ({it:R}, {it:r}).
{cmd:matcproc} returns {it:T} and {it:a} such that every {it:b} = {it:bc} {it:T}' + {it:a} obeys the
constraints, and {it:C}, the constraint matrix.
{p_end}

{pstd}
{cmd:makecns} {it:clist} builds {cmd:e(Cns)} from the constraints numbered in {it:clist}, defined
with {cmd:constraint}. A constraint that cannot be read, or that is redundant with or contradicts the
earlier ones, is left out with a note ({cmd:nocnsnotes} hides the notes). {cmd:displaycns} lists the
constraints used; {cmd:r} returns them in {cmd:r(cns}{it:#}{cmd:)}, with {cmd:r(k)} and
{cmd:r(clist)}. {cmd:matrix dispCns} lists the constraints of {cmd:e(Cns)}; with {cmd:, r} it
returns them instead.
{p_end}

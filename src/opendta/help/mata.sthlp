{smcl}
{title:Title}

{p2colset 5 17 19 2}{...}
{p2col :{cmd:mata} {hline 2}}The matrix programming language{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:mata}{space 6}(block, ended by {cmd:end}; errors are reported and the block goes on){p_end}
{p 8 16 2}{cmd:mata:}{space 5}(block, ended by {cmd:end}; the first error stops it){p_end}
{p 8 16 2}{cmd:mata:} {it:statement}{p_end}
{p 8 16 2}{cmd:mata} {c -(}{cmd:clear}|{cmd:describe}|{cmd:drop} {it:name}|{it:name}{cmd:()}{c )-}{p_end}

{title:Description}

{pstd}
Mata works with matrices of real numbers, complex numbers, strings and
pointers. Typing an expression shows its value; {cmd:x = }{it:exp} stores it.
Variables typed at the Mata prompt live until {cmd:mata clear}.
{p_end}

{title:Operators}

{p2colset 9 28 30 2}{...}
{p2col :{cmd:+ - * / ^}}arithmetic ({cmd:*} is the matrix product){p_end}
{p2col :{cmd::+ :- :* :/ :^}}element by element{p_end}
{p2col :{cmd:'}}transpose; {cmd:#} Kronecker product{p_end}
{p2col :{cmd:,}  {cmd:\}}join columns, join rows: {cmd:(1,2\3,4)}{p_end}
{p2col :{cmd:a..b}  {cmd:a::b}}row and column sequences{p_end}
{p2col :{cmd:== != > >= < <=}}comparison; {cmd::==} etc. element by element{p_end}
{p2col :{cmd:& | !  && ||}}logic; {cmd:?:} conditional{p_end}
{p2col :{cmd:A[i,j]  A[|i,j \ k,l|]}}subscripts; {cmd:.} means all rows or columns{p_end}
{p2col :{cmd:&x  *p  s.x  p->x}}pointers and struct members{p_end}
{p2colreset}{...}

{title:Programming}

{p 8 12 2}
{it:type} {it:org} {it:name}{cmd:(}{it:args}{cmd:)} {cmd:{c -(}} ... {cmd:{c )-}}, with
{it:type} {cmd:real}, {cmd:string}, {cmd:complex}, {cmd:pointer}, {cmd:transmorphic}
or {cmd:void} and {it:org} {cmd:scalar}, {cmd:vector}, {cmd:rowvector},
{cmd:colvector} or {cmd:matrix}. Arguments after {cmd:|} are optional
({cmd:args()} counts them). Arguments are passed by reference.
{p_end}
{p 8 12 2}
{cmd:if}/{cmd:else}, {cmd:for}, {cmd:while}, {cmd:do}...{cmd:while}, {cmd:break},
{cmd:continue}, {cmd:return()}, {cmd:struct}.
{p_end}

{title:Functions}

{pstd}
Size and type: {cmd:rows cols length eltype orgtype isreal isstring}.
Creation: {cmd:J I e range rangen runiform rnormal}.
Math: the Stata functions ({cmd:sqrt ln exp normal invnormal round mod}...)
element by element; {cmd:sum rowsum colsum max min rowmax colmax mean
variance correlation cross runningsum}.
Matrices: {cmd:invsym luinv cholinv pinv cholesky det rank trace diag
diagonal lusolve cholsolve qrsolve symeigensystem eigensystem svd norm
sort order uniqrows select selectindex rowshape colshape vec vech}.
Missing: {cmd:missing nonmissing hasmissing editmissing editvalue}.
Strings: {cmd:strlen substr strupper strlower strtrim strpos subinstr
strofreal strtoreal tokens invtokens}; output: {cmd:printf sprintf display
errprintf}.
{p_end}

{pstd}
Stata interface: {cmd:st_data st_sdata st_view st_store st_sstore st_addvar
st_addobs st_nobs st_nvar st_varindex st_varname st_local st_global
st_numscalar st_strscalar st_matrix st_matrixrowstripe st_matrixcolstripe
st_vartype st_varformat st_varlabel stata _stata}.
{cmd:st_view()} returns a copy of the data in OpenDTA: changing it does not
change the dataset; use {cmd:st_store()}.
{p_end}

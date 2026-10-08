{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:matrix} {hline 2}}Matrices{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmdab:mat:rix} [{cmd:define}] {it:A} {cmd:=} {it:matrix expression}

{p 8 16 2}
{cmdab:mat:rix} {cmd:list} {it:A} [{cmd:,} {cmd:format(}{it:%fmt}{cmd:)} {cmd:title(}{it:text}{cmd:)} {cmd:nonames}]

{p 8 16 2}
{cmdab:mat:rix} {c -(}{cmd:rownames}|{cmd:colnames}{c )-} {it:A} {cmd:=} {it:names}{space 4}{cmdab:mat:rix} {cmd:drop} {it:A}|{cmd:_all}

{title:Expressions}

{pstd}
{cmd:A + B}, {cmd:A - B}, {cmd:A * B} (matrix or scalar product), {cmd:A # B}
(Kronecker), {cmd:A'} (transpose), {cmd:A, B} (side by side),
{cmd:A \ B} (one above the other), {cmd:(1,2\3,4)}, submatrices
{cmd:A[1..2, 2...]} and the functions {cmd:J()}, {cmd:I()}, {cmd:inv()},
{cmd:invsym()}, {cmd:cholesky()}, {cmd:diag()}, {cmd:vecdiag()},
{cmd:hadamard()}, {cmd:corr()}, {cmd:vec()} and {cmd:nullmat()}.
{p_end}

{pstd}
In ordinary expressions: {cmd:A[i,j]}, {cmd:el(A,i,j)}, {cmd:rowsof()},
{cmd:colsof()}, {cmd:trace()}, {cmd:det()}, {cmd:issymmetric()}.
{p_end}

{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:generate} {hline 2}}Create or change the contents of a variable{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmdab:g:enerate} [{it:type}] {it:newvar} {cmd:=} {it:exp} [{it:if}] [{it:in}]

{p 8 16 2}
{cmd:replace} {it:oldvar} {cmd:=} {it:exp} [{it:if}] [{it:in}] [{cmd:,} {cmdab:nop:romote}]

{title:Description}

{pstd}
{cmd:generate} creates a new variable from an expression; observations
excluded by {it:if} or {it:in} get missing. The default type is
{cmd:float}; use {cmd:double} for identifiers and large integers.
{p_end}

{pstd}
{cmd:replace} changes the values of an existing variable. When the new
values do not fit the storage type, the type is promoted (for example,
{cmd:byte} to {cmd:int}) unless {cmd:nopromote} is given.
Expressions may use {cmd:_n}, {cmd:_N} and subscripts such as {cmd:x[_n-1]}.
{p_end}

{title:Examples}

        {cmd:. generate imc = peso / (altura/100)^2}
        {cmd:. replace imc = . if altura == 0}
        {cmd:. bysort familia (idade): generate ordem = _n}

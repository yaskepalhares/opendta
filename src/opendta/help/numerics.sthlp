{smcl}
{title:Title}

{p2colset 5 26 28 2}{...}
{p2col :{cmd:set numerics} {hline 2}}Stata arithmetic or more precise arithmetic{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:set numerics} {c -(}{cmd:stata}|{cmd:precise}{c )-} [{cmd:,} {cmdab:perm:anently}]

{title:Description}

{pstd}
With {cmd:set numerics stata} (the default), OpenDTA computes the way Stata 14
does, including its rounding, so that results match Stata digit for digit.
{p_end}

{pstd}
{cmd:set numerics precise} keeps the same estimators and the same output, but
computes them more exactly:
{p_end}

{p 8 12 2}
- {cmd:regress} solves the least-squares problem by a QR decomposition of X
instead of inverting X'X; the error grows with the condition number of X, not
with its square (polynomials and variables on very different scales gain
several digits);{p_end}
{p 8 12 2}
- {cmd:D.} and {cmd:S.} of float variables stay in double precision (Stata
rounds them to float);{p_end}
{p 8 12 2}
- {cmd:logit}, {cmd:probit}, {cmd:poisson} and Mata's {cmd:optimize()} take
extra Newton steps after convergence, so the estimates sit at the exact
maximum (the iteration log does not change);{p_end}
{p 8 12 2}
- numerical derivatives in {cmd:optimize()} use Richardson extrapolation;{p_end}
{p 8 12 2}
- new variables are {cmd:double} when the default type ({cmd:set type}) is
{cmd:float}, so 0.1 is stored as 0.1 and not as 0.1000000015.{p_end}

{pstd}
Results usually agree with Stata to 6 or more digits; they differ only where
Stata's arithmetic loses precision. {cmd:permanently} remembers the choice the
next time OpenDTA opens.
{p_end}

{title:Example}

        {cmd:. set numerics precise}
        {cmd:. regress y c.x##c.x##c.x}

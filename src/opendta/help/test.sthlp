{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:test} {hline 2}}Tests and combinations of coefficients{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:test} {it:coeflist}  |  {cmd:test} {it:exp} {cmd:=} {it:exp} [{cmd:=} ...] [{cmd:,} {opt a:ccumulate}]{p_end}
{p 8 16 2}{cmd:testparm} {it:varlist} [{cmd:,} {opt eq:ual}]{p_end}
{p 8 16 2}{cmd:lincom} {it:exp} [{cmd:,} {opt or} {opt irr} {opt rrr} {opt eform} {opt l:evel(#)}]{p_end}
{p 8 16 2}{cmd:nlcom} [{it:name}{cmd::}] {it:exp} [{cmd:,} {opt l:evel(#)}]{p_end}

{title:Description}

{pstd}
{cmd:test} performs Wald tests of linear hypotheses about the last estimates (F after regress, chi2 after maximum-likelihood models). {cmd:testparm} tests that all coefficients of a varlist (factor variables allowed) are zero. {cmd:lincom} estimates a linear combination with its standard error; {cmd:nlcom} a nonlinear combination by the delta method.
{p_end}

{title:Examples}

{phang2}{cmd:. regress y x1 x2 i.g}{p_end}
{phang2}{cmd:. test x1 = x2}{p_end}
{phang2}{cmd:. testparm i.g}{p_end}
{phang2}{cmd:. lincom x1 + x2}{p_end}
{phang2}{cmd:. nlcom ratio: _b[x1]/_b[x2]}{p_end}

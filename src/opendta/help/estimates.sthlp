{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:estimates} {hline 2}}Store, restore and compare estimation results{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:estimates store} {it:name}{p_end}
{p 8 16 2}{cmd:estimates restore} {it:name}{p_end}
{p 8 16 2}{cmd:estimates replay} [{it:namelist}]{p_end}
{p 8 16 2}{cmd:estimates dir}  |  {cmd:estimates drop} {it:namelist}  |  {cmd:estimates clear}{p_end}
{p 8 16 2}{cmd:estimates table} [{it:namelist}] [{cmd:,} {opt b}[{cmd:(}{it:%fmt}{cmd:)}] {opt se} {opt t} {opt p} {opt star} {opt stats(scalarlist)} {opt keep()} {opt drop()}]{p_end}

{title:Description}

{pstd}
Keeps estimation results in memory under a name, to compare models side by side with {cmd:estimates table} or to test them with {help lrtest} and {help hausman}.
{p_end}

{title:Examples}

{phang2}{cmd:. regress y x1}{p_end}
{phang2}{cmd:. estimates store m1}{p_end}
{phang2}{cmd:. regress y x1 x2}{p_end}
{phang2}{cmd:. estimates store m2}{p_end}
{phang2}{cmd:. estimates table m1 m2, se stats(N r2)}{p_end}

{smcl}
{title:Title}

{p 4 4 2}
{bf:fvset} {hline 2} Declare factor-variable settings

{title:Syntax}

{p 8 16 2}{cmd:fvset} {cmd:base} {it:base_spec} {it:varlist}{p_end}
{p 8 16 2}{cmd:fvset} {cmd:design} {it:design_spec} {it:varlist}{p_end}
{p 8 16 2}{cmd:fvset} {cmd:clear} {it:varlist}{p_end}
{p 8 16 2}{cmd:fvset} {cmd:report} [{it:varlist}] [{cmd:,} {opt base(base_spec)} {opt design(design_spec)}]{p_end}

{p 4 4 2}
{it:base_spec} is {cmd:default}, {cmd:first}, {cmd:last}, {cmd:frequent}, {cmd:none} or a level
{it:#}; {it:design_spec} is {cmd:default}, {cmd:asbalanced} or {cmd:asobserved}.

{title:Description}

{pstd}
{cmd:fvset base} sets the base level that {cmd:i.}{it:varname} uses when the term gives none
(an {cmd:ib}{it:#}{cmd:.} operator still wins); {cmd:fvset design} records the design that
{cmd:margins} reads; {cmd:fvset clear} removes both, and {cmd:fvset report} lists them. The
settings are characteristics of the variables ({cmd:_fv_base}, {cmd:_fv_design}) and are saved
with the data.
{p_end}

{title:Examples}

{phang2}{cmd:. fvset base 3 region}{p_end}
{phang2}{cmd:. regress y i.region}{p_end}
{phang2}{cmd:. fvset report}{p_end}

{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:datautil} {hline 2}}Creating and reshaping variables and observations{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:clonevar} {it:newvar} {cmd:=} {it:varname} [{it:if}] [{it:in}]{p_end}
{p 8 16 2}{cmd:insobs} {it:#} [{cmd:,} {opt b:efore(#)} | {opt a:fter(#)}]{p_end}
{p 8 16 2}{cmd:range} {it:varname} {it:#first} {it:#last} [{it:#obs}]{p_end}
{p 8 16 2}{cmd:separate} {it:varname} [{it:if}] [{it:in}] {cmd:,} {opt by(groupvar|exp)} [{opt gen:erate(stubname)} {opt seq:uential} {opt miss:ing} {opt shortl:abel}]{p_end}
{p 8 16 2}{cmd:ipolate} {it:yvar} {it:xvar} [{it:if}] [{it:in}] {cmd:,} {opt gen:erate(newvar)} [{opt e:polate}]{p_end}
{p 8 16 2}{cmd:expandcl} [{cmd:=}]{it:exp} [{it:if}] [{it:in}] {cmd:,} {opt cl:uster(varlist)} {opt gen:erate(newvar)}{p_end}
{p 8 16 2}{cmd:stack} {it:varlist} [{it:if}] [{it:in}] {cmd:,} {opt into(newvars)} | {opt g:roup(#)} [{opt clear} {opt wide} {opt keep(varlist)}]{p_end}
{p 8 16 2}{cmd:xpose} {cmd:,} {opt clear} [{opt for:mat} | {opt format(%fmt)} {opt varn:ame} {opt pro:mote}]{p_end}

{title:Description}

{pstd}
{cmd:clonevar} creates an exact copy of a variable: storage type, display format, value and variable labels, notes and characteristics. {cmd:insobs} inserts empty observations at the end or at a position; {cmd:range} creates a variable going from {it:#first} to {it:#last} in equal steps (with {it:#obs}, setting the number of observations); {cmd:separate} creates one variable per group; {cmd:ipolate} interpolates (and with {opt epolate} extrapolates) {it:yvar} linearly on {it:xvar}, also with {cmd:by}; {cmd:expandcl} duplicates whole clusters and numbers the new ones; {cmd:stack} stacks groups of variables into new variables, with {cmd:_stack} identifying the group; {cmd:xpose} transposes the data, observations becoming variables.
{p_end}

{title:Examples}

{phang2}{cmd:. clonevar income2 = income}{p_end}
{phang2}{cmd:. range t 0 1 101}{p_end}
{phang2}{cmd:. separate weight, by(sex)}{p_end}
{phang2}{cmd:. by id: ipolate bmi age, gen(bmi_i)}{p_end}
{phang2}{cmd:. stack x1 y1 x2 y2, into(x y) clear}{p_end}

{smcl}
{title:Title}

{p2colset 5 19 21 2}{...}
{p2col :{cmd:recode} {hline 2}}Recode values and convert variables{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:recode} {it:varlist} {cmd:(}{it:rule}{cmd:)} [{cmd:(}{it:rule}{cmd:)} ...] [{it:if}] [{it:in}] [{cmd:,} {opt gen:erate(newvars)} {opt pre:fix(str)} {opt lab:el(name)} {opt copy:rest}]{p_end}
{p 8 16 2}{cmd:encode} {it:strvar} {cmd:,} {opt gen:erate(newvar)} [{opt label(name)} {opt noextend}]{p_end}
{p 8 16 2}{cmd:decode} {it:numvar} {cmd:,} {opt gen:erate(newvar)}{p_end}
{p 8 16 2}{cmd:destring} [{it:varlist}] {cmd:,} {c -(}{opt gen:erate()}|{opt replace}{c )-} [{opt ignore("chars")} {opt force} {opt float} {opt percent} {opt dpcomma}]{p_end}
{p 8 16 2}{cmd:tostring} {it:varlist} {cmd:,} {c -(}{opt gen:erate()}|{opt replace}{c )-} [{opt force} {opt format(%fmt)}]{p_end}
{p 8 16 2}{cmd:split} {it:strvar} [{cmd:,} {opt gen:erate(stub)} {opt p:arse(strings)} {opt l:imit(#)} {opt destring}]{p_end}
{p 8 16 2}{cmd:mvdecode} {it:varlist} {cmd:,} {opt mv(numlist [= .a])}  |  {cmd:mvencode} {it:varlist} {cmd:,} {opt mv(#)} [{opt override}]{p_end}

{title:Rules of recode}

{p2colset 9 30 32 2}{...}
{p2col :{cmd:(1 2 = 1)}}listed values{p_end}
{p2col :{cmd:(3/5 = 2 "mid")}}a closed interval, with a value label{p_end}
{p2col :{cmd:(min/0 = 0)}, {cmd:(10/max = 9)}}smallest and largest values{p_end}
{p2col :{cmd:(missing = 0)}, {cmd:(nonmissing = 1)}}{p_end}
{p2col :{cmd:(else = .)}}everything not matched before{p_end}
{p2colreset}{...}

{pstd}
The first rule that matches an observation wins.
{p_end}

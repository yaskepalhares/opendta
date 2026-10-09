{smcl}
{title:Title}

{p2colset 5 21 23 2}{...}
{p2col :{cmd:collapse} {hline 2}}Replace the data by group statistics{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:collapse} [{cmd:(}{it:stat}{cmd:)}] {it:varlist} [[{cmd:(}{it:stat}{cmd:)}] {it:newvar}{cmd:=}{it:var} ...] [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt by(varlist)} {opt cw}]

{p 8 16 2}{cmd:contract} {it:varlist} [{cmd:,} {opt f:req(name)} {opt cf:req(name)} {opt p:ercent(name)} {opt cp:ercent(name)} {opt z:ero} {opt nomiss}]{p_end}
{p 8 16 2}{cmd:expand} [{cmd:=}]{it:exp} [{it:if}] [{it:in}] [{cmd:,} {opt g:enerate(newvar)}]{p_end}
{p 8 16 2}{cmd:fillin} {it:varlist}{p_end}

{title:Description}

{pstd}
{cmd:collapse} keeps one observation per group of {opt by()} with the
requested statistics: {cmd:mean} (default), {cmd:median}, {cmd:p1}-{cmd:p99},
{cmd:sd}, {cmd:semean}, {cmd:sebinomial}, {cmd:sepoisson}, {cmd:sum},
{cmd:rawsum}, {cmd:count}, {cmd:percent}, {cmd:min}, {cmd:max}, {cmd:iqr},
{cmd:first}, {cmd:last}, {cmd:firstnm} and {cmd:lastnm}.
{p_end}

{pstd}
{cmd:contract} keeps one observation per combination with its frequency;
{cmd:expand} duplicates observations ({it:exp} copies, new ones at the
end); {cmd:fillin} adds the missing combinations of {it:varlist} and marks
them with {cmd:_fillin}.
{p_end}

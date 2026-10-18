{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:epitab} {hline 2}}Tables for epidemiologists{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:cs} {it:casevar} {it:exposedvar} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt by(varlist)} {opt or} {opt tb} {opt w:oolf} {opt ex:act} {opt l:evel(#)} {opt nocrude} {opt nohom}]{p_end}
{p 8 16 2}{cmd:csi} {it:#a #b #c #d} [{cmd:,} {opt or} {opt tb} {opt w:oolf} {opt ex:act} {opt l:evel(#)}]{p_end}
{p 8 16 2}{cmd:cc} {it:casevar} {it:exposedvar} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt by(varlist)} {opt ex:act} {opt corn:field} {opt w:oolf} {opt tb} {opt l:evel(#)} {opt nocrude} {opt nohom}]{p_end}
{p 8 16 2}{cmd:cci} {it:#a #b #c #d} [{cmd:,} {opt ex:act} {opt corn:field} {opt w:oolf} {opt tb} {opt l:evel(#)}]{p_end}
{p 8 16 2}{cmd:ir} {it:casevar} {it:exposedvar} {it:timevar} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt by(varlist)} {opt tb} {opt ex:act} {opt l:evel(#)}]{p_end}
{p 8 16 2}{cmd:iri} {it:#a #b #N1 #N2} [{cmd:,} {opt tb} {opt ex:act} {opt l:evel(#)}]{p_end}
{p 8 16 2}{cmd:mcc} {it:var_exposed_case} {it:var_exposed_control} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt tb} {opt l:evel(#)}]{p_end}
{p 8 16 2}{cmd:mcci} {it:#a #b #c #d} [{cmd:,} {opt tb} {opt l:evel(#)}]{p_end}
{p 8 16 2}{cmd:tabodds} {it:casevar} {it:expvar} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt or} {opt adj:ust(varlist)} {opt ba:se(#)} {opt gr:aph} {opt ci:plot} {opt l:evel(#)}]{p_end}
{p 8 16 2}{cmd:mhodds} {it:casevar} {it:expvar} [{it:adjvars}] [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt bin:ary} {opt com:pare(v1,v2)} {opt l:evel(#)}]{p_end}

{title:Description}

{pstd}
{cmd:cs} (cohort studies) reports the risks, the risk difference, the risk ratio and the attributable or prevented fractions; {cmd:cc} (case-control studies) the odds ratio, by default with exact limits; {cmd:ir} (incidence-rate data) the rate difference and the rate ratio, with exact limits and mid-p values; {cmd:mcc} (matched case-control pairs) McNemar's test and the comparisons of the proportions with the factor. The commands ending in {cmd:i} take the four cells of the table on the command line: cases exposed, cases unexposed, noncases (controls) exposed, noncases unexposed; {cmd:iri} takes the cases and the person-time of the exposed and the unexposed; {cmd:mcci} takes the pairs with both exposed, case only, control only and neither.
{p_end}

{pstd}
With {opt by()}, the tables are stratified: one estimate per stratum, the crude estimate, the Mantel-Haenszel combined estimate and a test of homogeneity; {cmd:cc} adds the Mantel-Haenszel test that the combined odds ratio is 1.
{p_end}

{pstd}
{cmd:tabodds} tabulates the odds of being a case by level of the exposure, with the test of homogeneity and the score test for trend; with {opt or}, odds ratios against a base level, adjusted with {opt adjust()} by Mantel-Haenszel. {cmd:mhodds} estimates the Mantel-Haenszel odds ratio controlling for other variables; with an exposure of more than two levels it reports the score test for trend and the approximate odds ratio per unit.
{p_end}

{title:Remarks}

{pstd}
The options {opt estandard}, {opt istandard}, {opt standard()}, {opt binomial} and {opt pool} of {cmd:cs}, {cmd:cc} and {cmd:ir}, and {opt by()} of {cmd:mhodds}, are not available yet.
{p_end}

{title:Examples}

{phang2}{cmd:. cs case exposed}{p_end}
{phang2}{cmd:. csi 7 12 9 2}{p_end}
{phang2}{cmd:. cc case exposed, by(agegrp)}{p_end}
{phang2}{cmd:. iri 41 15 28010 19017}{p_end}
{phang2}{cmd:. mcci 8 8 3 8}{p_end}
{phang2}{cmd:. tabodds case alcohol, or adjust(age)}{p_end}
{phang2}{cmd:. mhodds case alcohol age}{p_end}

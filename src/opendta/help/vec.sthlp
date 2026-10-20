{smcl}
{title:Title}

{p 4 4 2}
{bf:vecrank, vec, svar} {hline 2} Cointegration and structural VAR

{title:Syntax}

{p 8 16 2}{cmd:vecrank} {it:depvarlist} [{it:if}] [{it:in}] [{cmd:,} {opt lags(#)} {opt trend(}{it:none}|{it:rconstant}|{it:constant}|{it:rtrend}|{it:trend}{cmd:)} {opt noreduce} {opt notrace} {opt max} {opt ic} {opt levela}]{p_end}
{p 8 16 2}{cmd:vec} {it:depvarlist} [{it:if}] [{it:in}] [{cmd:,} {opt rank(#)} {opt lags(#)} {opt trend(...)} {opt noreduce} {opt level(#)} {opt nobtable} {opt noidtest} {opt alpha} {opt pi(name)} {opt noptable} {opt mai} {opt noetable}]{p_end}
{p 8 16 2}{cmd:svar} {it:depvarlist} [{it:if}] [{it:in}] [{cmd:,} {opt lags(numlist)} {opt exog(varlist)} {opt noconstant} {opt aeq(matrix)} {opt beq(matrix)} {opt cholesky} {opt small} {opt dfk} {opt level(#)} {opt nolog}]{p_end}

{title:Description}

{pstd}
{cmd:vecrank} reports the Johansen trace (and maximum-eigenvalue) statistics, with the information criteria of {cmd:ic}. {cmd:vec} fits a vector error-correction model of the chosen rank, with the Johansen normalization of the cointegrating equations. {cmd:svar} fits a short-run structural VAR with the A and B matrices of {cmd:aeq()} and {cmd:beq()} (missing elements are free). {cmd:irf create} works after {cmd:var}, {cmd:svar} and {cmd:vec}.
{p_end}

{pstd}
Not available: {cmd:bconstraints()}, {cmd:aconstraints()} and {cmd:sindicators()} in {cmd:vec}; {cmd:aconstraints()}, {cmd:bconstraints()} and {cmd:lrconstraints()} in {cmd:svar}; the Osterwald-Lenum critical values are tabulated only for {cmd:trend(constant)} and {cmd:trend(rconstant)}, and the other trends show a dot in the critical-value columns.
{p_end}

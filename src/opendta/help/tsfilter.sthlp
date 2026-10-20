{smcl}
{title:Title}

{p 4 4 2}
{bf:tsfilter, tssmooth} {hline 2} Filters and smoothers for time series

{title:Syntax}

{p 8 16 2}{cmd:tsfilter} {c -(}{opt hp}|{opt bk}|{opt cf}|{opt bw}{c )-} {it:newvars} {cmd:=} {it:varlist} [{it:if}] [{it:in}] [{cmd:,} {opt smooth(#)} {opt minperiod(#)} {opt maxperiod(#)} {opt smaorder(#)} {opt order(#)} {opt stationary} {opt drift} {opt trend(newvars)} {opt gain(newvars)}]{p_end}
{p 8 16 2}{cmd:tssmooth} {c -(}{opt ma}|{opt exponential}|{opt dexponential}|{opt hwinters}|{opt shwinters}|{opt nl}{c )-} {it:newvar} {cmd:=} {it:exp} [{it:if}] [{it:in}] [{cmd:,} {opt replace} {opt forecast(#)} {opt window(#)} {opt weights(...)} {opt parms(...)} {opt samp0(#)} {opt s0(...)} {opt period(#)} {opt additive} {opt smoother(...)}]{p_end}

{title:Description}

{pstd}
{cmd:tsfilter} extracts the cyclical component with the Hodrick-Prescott, Baxter-King, Christiano-Fitzgerald or Butterworth filter; the trend and the gain go to the new variables of {cmd:trend()} and {cmd:gain()}. {cmd:tssmooth} applies moving averages, single and double exponential smoothing, Holt-Winters seasonal and nonseasonal smoothing and the nonlinear smoothers of {cmd:smooth}. The coefficients of the recursive smoothers are chosen by minimizing the sum of squared one-step errors. Both need {cmd:tsset} data.
{p_end}

{title:Examples}

{phang2}{cmd:. tsfilter hp gdp_c = gdp, smooth(1600) trend(gdp_t)}{p_end}
{phang2}{cmd:. tssmooth ma sm = gdp, window(2 1 2)}{p_end}
{phang2}{cmd:. tssmooth exponential ex = gdp}{p_end}

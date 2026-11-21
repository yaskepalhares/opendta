{smcl}
{title:Title}

{p 4 4 2}
{bf:forecast} {hline 2} Econometric model forecasting

{title:Syntax}

{p 8 16 2}{cmd:forecast create} [{it:name}] [{cmd:,} {opt replace}]{p_end}
{p 8 16 2}{cmd:forecast estimates} {it:estname} [{cmd:,} {opt names(newvar)}]{p_end}
{p 8 16 2}{cmd:forecast identity} {it:varname} {cmd:=} {it:exp} [{cmd:,} {opt gen:erate}]{p_end}
{p 8 16 2}{cmd:forecast exogenous} {it:varlist}{p_end}
{p 8 16 2}{cmd:forecast adjust} {it:varname} {cmd:=} {it:exp}{p_end}
{p 8 16 2}{cmd:forecast solve} [{cmd:,} {opt pre:fix(str)} {opt suf:fix(str)} {opt beg:in(time)} {opt end(time)} {opt per:iods(#)} {opt stat:ic} {opt act:uals} {opt log(off)} {opt double}]{p_end}
{p 8 16 2}{cmd:forecast describe}  |  {cmd:forecast list}  |  {cmd:forecast query}  |  {cmd:forecast drop} [{cmd:,} {opt pre:fix()} {opt suf:fix()}]  |  {cmd:forecast clear}{p_end}

{title:Description}

{pstd}
A forecast model gathers estimation results stored with {cmd:estimates store}, identities and
exogenous variables. {cmd:forecast solve} computes, period by period from {cmd:begin()} to
{cmd:end()}, the values of the endogenous variables (Gauss-Seidel within each period) into new
variables named with the prefix {cmd:f_}; dynamic forecasts use the forecasts of the lags, and
{cmd:static} the actual values. {cmd:forecast adjust} changes the value of an equation after it
is computed (for example {cmd:forecast adjust y = y + 1}). The data must be {cmd:tsset}.
{p_end}

{pstd}
Not available yet: {cmd:forecast coefvector}, multiple-equation estimation results, panel data,
{cmd:predict()} and stochastic simulation ({cmd:simulate()}).
{p_end}

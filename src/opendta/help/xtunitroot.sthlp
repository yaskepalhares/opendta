{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:xtunitroot} {hline 2}}Panel data: unit-root tests{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:xtunitroot} {it:test} {it:varname} [{it:if}] [{it:in}] [{cmd:,} {it:options}]{p_end}

{title:Description}

{pstd}
{cmd:xtunitroot} performs panel unit-root tests. Available tests: {cmd:ht} (Harris-Tzavalis), {cmd:breitung}, {cmd:hadri} (stationarity null) and {cmd:fisher} (augmented Dickey-Fuller tests of every panel, combined by Maddala-Wu and Choi). {cmd:llc} and {cmd:ips} are not available yet.
{p_end}

{title:Options}

{phang}
{opt trend} include a time trend.
{p_end}

{phang}
{opt noconstant} omit the panel means.
{p_end}

{phang}
{opt dfuller} augmented Dickey-Fuller tests of every panel; required by {cmd:fisher}.
{p_end}

{phang}
{opt lags(#)} lags of the ADF regressions ({cmd:fisher}).
{p_end}

{phang}
{opt robust} heteroskedasticity-robust statistic ({cmd:hadri}).
{p_end}

{title:Examples}

{phang2}{cmd:. xtset id year}{p_end}
{phang2}{cmd:. xtunitroot ht y}{p_end}
{phang2}{cmd:. xtunitroot fisher y, dfuller lags(2)}{p_end}
{phang2}{cmd:. xtunitroot hadri y, trend}{p_end}

{title:Also see}

{pstd}
See also {help dfuller}.
{p_end}

{smcl}
{title:Title}

{p 4 4 2}
{bf:pergram, cumsp, xcorr, wntestb, wntestq, psdensity} {hline 2} Spectral analysis and white-noise tests

{title:Syntax}

{p 8 16 2}{cmd:pergram} {it:varname} [{it:if}] [{it:in}] [{cmd:,} {opt gen:erate(newvar)} {opt nograph} {it:graph options}]{p_end}
{p 8 16 2}{cmd:cumsp} {it:varname} [{it:if}] [{it:in}] [{cmd:,} {opt gen:erate(newvar)} {opt nograph} {it:graph options}]{p_end}
{p 8 16 2}{cmd:xcorr} {it:var1 var2} [{it:if}] [{it:in}] [{cmd:,} {opt lags(#)} {opt gen:erate(newvar)} {opt table} {opt noplot} {it:graph options}]{p_end}
{p 8 16 2}{cmd:wntestb} {it:varname} [{it:if}] [{it:in}] [{cmd:,} {opt table} {opt level(#)} {opt noplot}]{p_end}
{p 8 16 2}{cmd:wntestq} {it:varname} [{it:if}] [{it:in}] [{cmd:,} {opt lags(#)}]{p_end}
{p 8 16 2}{cmd:psdensity} {it:newvars} [{cmd:,} {opt pspectrum} {opt range(#)} {opt cycle(#)}]{p_end}

{title:Description}

{pstd}
{cmd:pergram} plots the log sample spectral density, {cmd:cumsp} the cumulative sample spectral distribution, {cmd:xcorr} the cross-correlogram of two series, {cmd:wntestb} the Bartlett cumulative-periodogram white-noise test and {cmd:wntestq} the portmanteau (Q) test. After {cmd:arima} or {cmd:arfima}, {cmd:psdensity} computes the parametric spectral density and the implied autocorrelation. {cmd:psdensity} after {cmd:ucm} is not available.
{p_end}

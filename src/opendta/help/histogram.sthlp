{smcl}
{title:Title}

{p2colset 5 24 26 2}{...}
{p2col :{cmd:histogram} {hline 2}}Histograms and kernel density estimates{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:histogram} {it:varname} [{it:if}] [{it:in}] [{it:fweight}] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:kdensity} {it:varname} [{it:if}] [{it:in}] [{it:fweight}|{it:aweight}] [{cmd:,} {it:options}]{p_end}

{title:histogram}

{p2colset 9 30 32 2}{...}
{p2col :{cmd:bin(#)}}number of bars{p_end}
{p2col :{cmd:width(#)}}width of the bars{p_end}
{p2col :{cmd:start(#)}}lower limit of the first bar (default: the minimum){p_end}
{p2col :{cmd:discrete}}one bar per value{p_end}
{p2col :{cmd:density} {cmd:fraction} {cmd:frequency} {cmd:percent}}scale of the y axis (default {cmd:density}){p_end}
{p2col :{cmd:normal}, {cmd:kdensity}}add a normal density or a kernel density{p_end}
{p2col :{cmd:addlabels}}write the heights above the bars{p_end}
{p2col :{cmd:horizontal}}horizontal bars{p_end}
{p2colreset}{...}

{pstd}
Without {cmd:bin()} and {cmd:width()}, the number of bars is
min(sqrt(N), 10 log10(N)), truncated to an integer. The message
{cmd:(bin=}{it:k}{cmd:, start=}{it:a}{cmd:, width=}{it:w}{cmd:)} reports the bars used. With
{cmd:by()}, all panels use the bars of the whole sample.
{p_end}

{title:kdensity}

{p2colset 9 34 36 2}{...}
{p2col :{cmd:kernel(}{it:name}{cmd:)}}{cmd:epanechnikov} (default), {cmd:epan2}, {cmd:biweight}, {cmd:cosine}, {cmd:gaussian}, {cmd:parzen}, {cmd:rectangle}, {cmd:triangle}{p_end}
{p2col :{cmd:bwidth(#)}}half-width of the kernel{p_end}
{p2col :{cmd:n(#)}}number of points (default min(N, 50)){p_end}
{p2col :{cmd:at(}{it:var}{cmd:)}}estimate at the values of {it:var}{p_end}
{p2col :{cmd:generate(}{it:newx newd}{cmd:)}}store the points and the estimates{p_end}
{p2col :{cmd:normal}}add a normal density{p_end}
{p2col :{cmd:nograph}}do not draw{p_end}
{p2colreset}{...}

{pstd}
The default bandwidth is 0.9 m / n^(1/5), with m = min(sd, IQR/1.349), the
width that would minimize the mean integrated squared error if the data were
normal. {cmd:kdensity} returns {cmd:r(bwidth)}, {cmd:r(n)} and {cmd:r(kernel)}.
{p_end}

{title:Examples}

{phang2}{cmd:. histogram mpg, frequency normal}{p_end}
{phang2}{cmd:. histogram rep78, discrete percent by(foreign)}{p_end}
{phang2}{cmd:. kdensity mpg, kernel(gaussian) generate(x d) nograph}{p_end}

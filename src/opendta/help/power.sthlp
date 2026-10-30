{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:power} {hline 2}}Power and sample-size analysis{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:power onemean} {it:m0} [{it:ma}] [{cmd:,} {opt sd(#)} {opt knownsd} {it:options}]{p_end}
{p 8 16 2}{cmd:power twomeans} {it:m1} [{it:m2}] [{cmd:,} {opt sd(#)} {opt sd1(#)} {opt sd2(#)} {opt knownsds} {it:options}]{p_end}
{p 8 16 2}{cmd:power pairedmeans} {it:ma1} [{it:ma2}] {cmd:,} {opt sddiff(#)} | {opt corr(#)} {opt sd(#)} [{opt sd1(#)} {opt sd2(#)}] [{it:options}]{p_end}
{p 8 16 2}{cmd:power oneproportion} {it:p0} [{it:pa}] [{cmd:,} {opt test(score|wald|binomial)} {it:options}]{p_end}
{p 8 16 2}{cmd:power twoproportions} {it:p1} [{it:p2}] [{cmd:,} {opt continuity} {it:options}]{p_end}
{p 8 16 2}{cmd:power onecorrelation} {it:r0} [{it:ra}] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:power twocorrelations} {it:r1} [{it:r2}] [{cmd:,} {it:options}]{p_end}
{p 8 16 2}{cmd:power twoway} [{it:means}] [{cmd:,} {opt varerror(#)} {opt nrows(#)} {opt ncols(#)} {opt factor(row|column|rowcol)} {opt cellw:eights(...)} {it:options}]{p_end}
{p 8 16 2}{cmd:power repeated} [{it:means}] [{cmd:,} {opt ngroups(#)} {opt nrepeated(#)} {opt corr(#)} {opt covmatrix(matname)} {opt factor(between|within|bwithin)} {it:options}]{p_end}
{p 8 16 2}{cmd:power cmh} {it:p1_1 p1_2 ...} [{cmd:,} {opt oratio(#)} {opt strw:eights(...)} {opt grratios(...)} {opt cont:inuity} {it:options}]{p_end}
{p 8 16 2}{cmd:power mcc} {it:p0} [{cmd:,} {opt oratio(#)} {opt corr(#)} {opt m(#)} {opt compare} {it:options}]{p_end}
{p 8 16 2}{cmd:power trend} {it:p_1 p_2 ...} [{cmd:,} {opt expos:ure(...)} {opt grweights(...)} {opt cont:inuity} {it:options}]{p_end}
{p 8 16 2}{cmd:power cox} [{it:b1}] [{cmd:,} {opt hratio(#)} {opt sd(#)} {opt r2(#)} {opt eventprob(#)} {opt wdprob(#)} {it:options}]{p_end}
{p 8 16 2}{cmd:power exponential} [{it:h1} [{it:h2}]] [{cmd:,} {opt hratio(#)} {opt hdiff:erence(#)} {opt loghaz:ard} {opt time(#)} {opt fperiod(#)} {opt aperiod(#)} {opt lossprob(#)} {opt uncond:itional} {it:options}]{p_end}
{p 8 16 2}{cmd:power logrank} [{it:s1} [{it:s2}]] [{cmd:,} {opt hratio(#)} {opt schoenfeld} {opt simpson(...)} {opt wdprob(#)} {it:options}]{p_end}

{p 4 4 2}{it:options}: {opt a:lpha(#)} {opt p:ower(#)} {opt b:eta(#)} {opt n(#)} {opt n1(#)}
{opt n2(#)} {opt nr:atio(#)} {opt diff(#)} {opt ones:ided} {opt dir:ection(upper|lower)}
{opt nfrac:tional} {opt tab:le}{p_end}

{title:Description}

{pstd}
{cmd:power} computes the sample size needed to reach a power (the default, with power 0.8
and alpha 0.05), the power of a given sample size ({cmd:n()}), or the smallest effect that
reaches a power with a given sample size ({cmd:n()} without the alternative value). Each
argument and option accepts a list of values in parentheses, such as {cmd:(102 104 106)};
the results then come as a table with one row per combination.
{p_end}

{pstd}
Means use the t test ({cmd:knownsd}: z test; unequal standard deviations in
{cmd:twomeans}: Satterthwaite's t test). Proportions use the score test in
{cmd:oneproportion} (also the Wald test and the exact binomial test) and Pearson's
chi-squared test in {cmd:twoproportions}. Correlations use Fisher's z transformation.
Sample sizes are rounded up unless {cmd:nfractional} is given.
{p_end}

{title:Examples}

{phang2}{cmd:. power onemean 100 102, sd(5)}{p_end}
{phang2}{cmd:. power twomeans 12 14, sd(4) power(.9)}{p_end}
{phang2}{cmd:. power twoproportions .3 .45, n(200)}{p_end}
{phang2}{cmd:. power onecorrelation 0 (.2 .3 .4)}{p_end}
{phang2}{cmd:. power twoway 134 143 91 \ 106 173 145, varerror(1417) power(.8)}{p_end}
{phang2}{cmd:. power cmh .2 .1 .3, oratio(2.5)}{p_end}
{phang2}{cmd:. power logrank .5, hratio(.7) wdprob(.1)}{p_end}

{title:Stored results}

{pstd}
{cmd:r(N)}, {cmd:r(N1)}, {cmd:r(N2)}, {cmd:r(N_per_group)}, {cmd:r(power)}, {cmd:r(beta)},
{cmd:r(alpha)}, {cmd:r(delta)}, the values of the parameters, {cmd:r(method)} and
{cmd:r(test)}.
{p_end}

{title:OpenDTA notes}

{pstd}
{cmd:onevariance}, {cmd:twovariances}, {cmd:pairedproportions}, {cmd:oneway}, {cmd:twoway},
{cmd:repeated}, {cmd:cmh}, {cmd:mcc}, {cmd:trend}, {cmd:cox}, {cmd:exponential} and
{cmd:logrank} are available, as are the {cmd:lrchi2} and {cmd:fisher} tests of
{cmd:twoproportions}. The {cmd:graph} option works for the first group of methods. In
{cmd:twoway}, {cmd:repeated}, {cmd:cmh}, {cmd:mcc}, {cmd:trend}, {cmd:cox}, {cmd:exponential}
and {cmd:logrank}, the options {cmd:graph()}, {cmd:saving()}, {cmd:show}, {cmd:st1()} and
{cmd:n1()}/{cmd:n2()} are not available yet.
{p_end}

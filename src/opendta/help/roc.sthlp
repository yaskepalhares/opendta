{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:roc} {hline 2}}Receiver operating characteristic (ROC) analysis{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:roctab} {it:refvar} {it:classvar} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt d:etail} {opt tab:le} {opt han:ley} {opt gr:aph} {opt noref:line} {opt l:evel(#)}]{p_end}
{p 8 16 2}{cmd:roccomp} {it:refvar} {it:classvar} [{it:classvars}] [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt by(varname)} {opt t:est(matname)} {opt gr:aph} {opt l:evel(#)}]{p_end}
{p 8 16 2}{cmd:rocfit} {it:refvar} {it:classvar} [{it:if}] [{it:in}] [{it:weight}] [{cmd:,} {opt cont:inuous(#)} {opt l:evel(#)} {opt nolog}]{p_end}
{p 8 16 2}{cmd:lroc} [{it:depvar}] [{cmd:,} {opt nog:raph}]{p_end}
{p 8 16 2}{cmd:lsens} [{it:depvar}] [{cmd:,} {opt genp:rob(var)} {opt gens:ens(var)} {opt gens:pec(var)} {opt rep:lace} {opt nog:raph}]{p_end}

{title:Description}

{pstd}
{cmd:roctab} computes the nonparametric ROC curve of a classifier against a 0/1 reference variable: the area under the curve with DeLong's standard error (or Hanley and McNeil's with {opt hanley}), the sensitivity, specificity and likelihood ratios at each cutpoint ({opt detail}) and the curve ({opt graph}). {cmd:roccomp} compares the areas of several classifiers measured on the same subjects (DeLong's covariance), or of one classifier in independent groups ({opt by()}), with a chi-squared test of equality or of the contrasts in {opt test()}.
{p_end}

{pstd}
{cmd:rocfit} fits the binormal model to ordinal ratings by maximum likelihood and reports the cutpoints, the separation and scale of the diseased distribution and the indices of the fitted curve (area, delta(m), d(e), d(a)). {cmd:lroc} and {cmd:lsens} work after {cmd:logit}, {cmd:logistic} and {cmd:probit}: the area under the ROC curve of the predicted probabilities, and the sensitivity and specificity against the probability cutoff.
{p_end}

{title:Remarks}

{pstd}
The options {opt binomial}, {opt bamber} and {opt lorenz} of {cmd:roctab}, {opt binormal} of {cmd:roccomp}, {opt all} and {opt beta()} of {cmd:lroc}, and the command {cmd:rocreg}, are not available yet.
{p_end}

{title:Examples}

{phang2}{cmd:. roctab disease rating, detail graph}{p_end}
{phang2}{cmd:. roccomp disease test1 test2, graph}{p_end}
{phang2}{cmd:. rocfit disease rating}{p_end}
{phang2}{cmd:. logit low age lwt smoke}{p_end}
{phang2}{cmd:. lroc}{p_end}

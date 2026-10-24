{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:rocreg} {hline 2}}Receiver operating characteristic (ROC) regression{p_end}
{p2colreset}{...}

{title:Syntax}

{pstd}Nonparametric ROC estimation{p_end}

{p 8 16 2}{cmd:rocreg} {it:refvar} {it:classvar} [{it:classvars}] [{it:if}] [{it:in}] [{cmd:,} {opt ctrlc:ov(varlist)} {opt ctrlm:odel(strata|linear)} {opt pvc(empirical|normal)} {opt tie:corrected} {opt auc} {opt roc(numlist)} {opt invr:oc(numlist)} {opt pauc(numlist)} {it:bootstrap_options}]{p_end}

{pstd}Parametric ROC estimation (probit GLM){p_end}

{p 8 16 2}{cmd:rocreg} {it:refvar} {it:classvar} [{it:if}] [{it:in}]{cmd:,} {opt pr:obit} [{opt rocc:ov(varlist)} {opt fprp:ts(#)} {opt ctrlf:prall} {opt ctrlc:ov(varlist)} {opt ctrlm:odel(strata|linear)} {opt pvc(empirical|normal)} {opt tie:corrected} {it:bootstrap_options}]{p_end}

{pstd}Parametric ROC estimation (maximum likelihood){p_end}

{p 8 16 2}{cmd:rocreg} {it:refvar} {it:classvar} [{it:if}] [{it:in}]{cmd:,} {opt pr:obit} {opt ml} [{opt rocc:ov(varlist)} {opt ctrlc:ov(varlist)} {opt cl:uster(varname)} {opt l:evel(#)}]{p_end}

{pstd}{it:bootstrap_options}: {opt nobootstrap} {opt bseed(#)} {opt breps(#)} {opt nobstrata} {opt cl:uster(varname)} {opt nodots} {opt l:evel(#)}{p_end}

{title:Description}

{pstd}
{cmd:rocreg} estimates the ROC curve of each classifier against the 0/1 reference variable, with the controls' distribution adjusted for the covariates in {opt ctrlcov()}: within their strata ({opt ctrlmodel(strata)}, the default) or through a linear regression of the controls' values ({opt ctrlmodel(linear)}). Each case is placed by its percentile value among the controls, taken from their empirical distribution or from a normal one ({opt pvc(normal)}); {opt tiecorrected} counts ties one half.
{p_end}

{pstd}
The nonparametric estimation reports the area under the curve ({opt auc}, the default), the true-positive rate at the false-positive rates in {opt roc()}, the false-positive rate at the true-positive rates in {opt invroc()} and the partial area up to the false-positive rates in {opt pauc()}. With {opt probit}, the ROC curve is modeled as Phi^-1 ROC(t) = {cmd:i_cons} + {cmd:s_cons} Phi^-1(t) + the effects of the covariates in {opt roccov()}, fit by probit regression over {opt fprpts()} false-positive rates (10 by default) or over all those of the controls ({opt ctrlfprall}).
{p_end}

{pstd}
The standard errors and the normal (N), percentile (P) and bias-corrected (BC) intervals come from {opt breps()} bootstrap replications (1000 by default), sampling cases and controls separately unless {opt nobstrata}; {opt bseed()} sets the random-number seed first.
{p_end}

{pstd}
With {opt probit ml}, the binormal model is fit by maximum likelihood: the controls' values are normal with a mean linear in the {opt ctrlcov()} variables, and the cases' values are normal with the mean and standard deviation that give the ROC curve above. The standard errors come from the information matrix ({opt cluster()} makes them robust to the clusters); without ROC covariates the AUC, Phi(i_cons / sqrt(1 + s_cons^2)), is reported with a delta-method standard error. In the bootstrap, {opt cluster()} draws whole clusters.
{p_end}

{title:Examples}

{phang2}{cmd:. rocreg disease test1 test2, bseed(1)}{p_end}
{phang2}{cmd:. rocreg disease test1, ctrlcov(age) roc(.2) pauc(.3) bseed(1)}{p_end}
{phang2}{cmd:. rocreg disease test1, probit roccov(male) bseed(1)}{p_end}

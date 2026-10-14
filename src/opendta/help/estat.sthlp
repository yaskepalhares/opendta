{smcl}
{title:Title}

{p 4 4 2}
{bf:estat} {hline 2} Postestimation statistics

{title:Syntax}

{p 8 16 2}{cmd:estat ic} [{cmd:,} {cmd:n(}{it:#}{cmd:)}]{p_end}
{p 8 16 2}{cmd:estat vce} [{cmd:,} {cmd:correlation} {cmd:format(}{it:%fmt}{cmd:)}]{p_end}
{p 8 16 2}{cmd:estat summarize} [{it:varlist}]{p_end}
{p 8 16 2}{cmd:estat hettest} [{it:varlist}] [{cmd:,} {cmd:rhs} {cmd:iid} {cmd:fstat}]  (after {cmd:regress}){p_end}
{p 8 16 2}{cmd:estat ovtest} [{cmd:,} {cmd:rhs}]  (after {cmd:regress}){p_end}
{p 8 16 2}{cmd:estat classification} [{cmd:,} {cmd:cutoff(}{it:#}{cmd:)}]  (after {cmd:logit}, {cmd:logistic}, {cmd:probit}){p_end}
{p 8 16 2}{cmd:estat gof} [{cmd:,} {cmd:group(}{it:#}{cmd:)}]  (after {cmd:logit}, {cmd:logistic}, {cmd:probit}){p_end}

{title:Description}

{phang}{cmd:estat ic}: AIC = -2 ll + 2 df and BIC = -2 ll + ln(N) df.{p_end}
{phang}{cmd:estat vce}: covariance (or correlation) matrix of the coefficients.{p_end}
{phang}{cmd:estat summarize}: mean, standard deviation, minimum and maximum of
the variables in the estimation sample.{p_end}
{phang}{cmd:estat hettest}: Breusch-Pagan / Cook-Weisberg test, by default with
the fitted values and the normal-errors version; {cmd:iid} uses N R2 from the
regression of the squared residuals.{p_end}
{phang}{cmd:estat ovtest}: Ramsey RESET test with the 2nd to 4th powers of the
standardized fitted values.{p_end}
{phang}{cmd:estat classification}: classification table and sensitivity,
specificity and predictive values at the cutoff (default .5).{p_end}
{phang}{cmd:estat gof}: Pearson chi-squared over covariate patterns, or the
Hosmer-Lemeshow test with {cmd:group()}.{p_end}

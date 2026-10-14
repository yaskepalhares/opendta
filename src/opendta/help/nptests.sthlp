{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:ranksum} {hline 2}}Nonparametric and normality tests{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:ranksum} {it:varname} [{it:if}] [{it:in}] {cmd:,} {opt by(groupvar)} [{opt porder}]{p_end}
{p 8 16 2}{cmd:signrank} {it:varname} {cmd:=} {it:exp} [{it:if}] [{it:in}]{p_end}
{p 8 16 2}{cmd:kwallis} {it:varname} [{it:if}] [{it:in}] {cmd:,} {opt by(groupvar)}{p_end}
{p 8 16 2}{cmd:spearman} [{it:varlist}] [{it:if}] [{it:in}] [{cmd:,} {opt stats(list)} {opt pw}]{p_end}
{p 8 16 2}{cmd:swilk} {it:varlist} [{it:if}] [{it:in}]{p_end}
{p 8 16 2}{cmd:sktest} {it:varlist} [{it:if}] [{it:in}] [{cmd:,} {opt noadj:ust}]{p_end}
{p 8 16 2}{cmd:sdtest} {it:varname} {cmd:==} {it:#} | {cmd:sdtest} {it:varname} {cmd:,} {opt by(groupvar)}{p_end}

{title:Description}

{pstd}
{cmd:ranksum}: Wilcoxon rank-sum (Mann-Whitney) test. {cmd:signrank}: Wilcoxon matched-pairs signed-rank test. {cmd:kwallis}: Kruskal-Wallis test. {cmd:spearman}: Spearman's rank correlation. {cmd:swilk}: Shapiro-Wilk normality test. {cmd:sktest}: skewness and kurtosis test for normality, with Royston's adjustment. {cmd:sdtest}: variance-comparison tests.
{p_end}

{title:Remarks}

{pstd}
swilk with tied values can differ slightly from Stata in the W statistic.
{p_end}

{title:Examples}

{phang2}{cmd:. ranksum mpg, by(foreign)}{p_end}
{phang2}{cmd:. kwallis bmi, by(region)}{p_end}
{phang2}{cmd:. swilk bmi}{p_end}
{phang2}{cmd:. sktest bmi}{p_end}

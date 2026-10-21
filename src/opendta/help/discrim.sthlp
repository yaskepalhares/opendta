{smcl}
{title:Title}

{p 4 4 2}
{bf:discrim}, {bf:candisc} {hline 2} Discriminant analysis

{title:Syntax}

{p 8 16 2}
{cmd:discrim} {cmd:lda} | {cmd:qda} | {cmd:logistic} varlist [if] [in] [weight] {cmd:,} {cmd:group(}varname{cmd:)} [{cmd:priors(}equal|proportional|numlist{cmd:)} {cmd:lootable} {cmd:notable} {cmd:nopriors}]

{p 8 16 2}
{cmd:discrim} {cmd:knn} varlist [if] [in] {cmd:,} {cmd:group(}varname{cmd:)} {cmd:k(}#{cmd:)} [{cmd:measure(}measure{cmd:)} {cmd:priors()} {cmd:lootable} {cmd:notable}]

{p 8 16 2}
{cmd:candisc} varlist [if] [in] [weight] {cmd:,} {cmd:group(}varname{cmd:)}

{title:Description}

{pstd}
{cmd:discrim lda} and {cmd:qda} classify by the largest posterior probability under multivariate normal groups with a pooled or a group-specific covariance matrix; {cmd:knn} uses the k nearest neighbors; {cmd:logistic} uses the multinomial logit fitted by {cmd:mlogit}. The leave-one-out table ({cmd:lootable}) is available after lda, qda and knn. {cmd:candisc} reports the canonical linear discriminant functions with Wilks' lambda tests.

{pstd}
Not implemented yet: {cmd:ties()}, {cmd:priors(}matname{cmd:)}, weights and {cmd:lootable} with {cmd:logistic}, {cmd:estat} options other than those listed below.

{title:Postestimation}

{phang2}{cmd:estat} {cmd:classtable} [{cmd:,} {cmd:loo} {cmd:nopriors}] | {cmd:errorrate} [{cmd:,} {cmd:pp} {cmd:nopriors}] | {cmd:list} [{cmd:,} {cmd:misclassified} {cmd:probabilities}] | {cmd:grdistances} (lda) | {cmd:covariance}{p_end}
{phang2}{cmd:predict} newvar [{cmd:,} {cmd:classification}]   {cmd:predict} newvarlist {cmd:,} {cmd:pr} | {cmd:mahalanobis}{p_end}

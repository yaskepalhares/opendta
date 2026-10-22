{smcl}
{title:Title}

{p 4 4 2}
{bf:teffects}, {bf:eteffects}, {bf:stteffects}, {bf:etregress}, {bf:etpoisson} {hline 2} Treatment effects

{title:Syntax}

{p 8 16 2}
{cmd:teffects} {c -(}{cmd:ra}|{cmd:ipw}|{cmd:ipwra}|{cmd:aipw}|{cmd:nnmatch}|{cmd:psmatch}{c )-} {cmd:(}ovar [omvarlist] [{cmd:,} omodel]{cmd:)} {cmd:(}tvar [tmvarlist] [{cmd:,} tmodel]{cmd:)} [if] [in] [{cmd:,} {cmd:ate}|{cmd:atet}|{cmd:pomeans} {cmd:nneighbor(}#{cmd:)} {cmd:metric(}{cmd:mahalanobis}|{cmd:euclidean}|{cmd:ivariance}{cmd:)} {cmd:vce(}{cmd:robust}[{cmd:, nn(}#{cmd:)}]{cmd:)}]

{p 8 16 2}
{cmd:eteffects} {cmd:(}ovar omvarlist{cmd:)} {cmd:(}tvar tmvarlist{cmd:)} [if] [in] [{cmd:,} {cmd:ate}|{cmd:atet}|{cmd:pomeans}]

{p 8 16 2}
{cmd:stteffects} {c -(}{cmd:ra}|{cmd:ipw}|{cmd:ipwra}{c )-} ... [{cmd:,} {cmd:ate}|{cmd:atet}|{cmd:pomeans}]

{p 8 16 2}
{cmd:etregress}|{cmd:etpoisson} depvar [indepvars] [if] [in] {cmd:, treat(}tvar {cmd:=} varlist{cmd:)} [{cmd:twostep}]

{title:Description}

{pstd}
Binary treatments. The outcome model is {cmd:linear}, {cmd:logit}, {cmd:probit} or {cmd:poisson} and the treatment model {cmd:logit} or {cmd:probit}; the variance is robust (stacked estimating equations, matching: Abadie and Imbens). {cmd:eteffects} takes a linear outcome and a probit treatment (control function). {cmd:stteffects} fits {cmd:exponential}, {cmd:weibull} or {cmd:lognormal} survival models and reports mean survival times. {cmd:etregress} is maximum likelihood or two-step, and {cmd:etpoisson} integrates the unobserved heterogeneity by Gauss-Hermite quadrature.

{pstd}
Not available: multivalued treatments, other outcome models, {cmd:stteffects wra}, {cmd:predict} after {cmd:teffects}.

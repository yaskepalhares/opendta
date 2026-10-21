{smcl}
{title:Title}

{p 4 4 2}
{bf:ivprobit}, {bf:ivtobit}, {bf:ivpoisson}, {bf:ivregress liml}, {bf:ivregress gmm} {hline 2} Models with endogenous regressors

{title:Syntax}

{p 8 16 2}
{cmd:ivprobit} depvar [varlist1] {cmd:(}varlist2 {cmd:=} varlist_iv{cmd:)} [if] [in] [weight] [{cmd:,} {cmd:vce(}vcetype{cmd:)} ...]

{p 8 16 2}
{cmd:ivtobit} depvar [varlist1] {cmd:(}varlist2 {cmd:=} varlist_iv{cmd:)} [if] [in] [weight] {cmd:,} {cmd:ll(}#{cmd:)}|{cmd:ul(}#{cmd:)} [...]

{p 8 16 2}
{cmd:ivpoisson} {c -(}{cmd:gmm}|{cmd:cfunction}{c )-} depvar [varlist1] {cmd:(}varlist2 {cmd:=} varlist_iv{cmd:)} [if] [in] [{cmd:,} {cmd:multiplicative} {cmd:onestep} {cmd:vce(}vcetype{cmd:)} ...]

{p 8 16 2}
{cmd:ivregress} {c -(}{cmd:2sls}|{cmd:liml}|{cmd:gmm}{c )-} depvar [varlist1] {cmd:(}varlist2 {cmd:=} varlist_iv{cmd:)} [if] [in] [weight] [{cmd:,} {cmd:wmatrix(}{cmd:robust}|{cmd:unadjusted}|{cmd:cluster} clustvar{cmd:)} {cmd:small} {cmd:vce(}vcetype{cmd:)}]

{title:Description}

{pstd}
{cmd:ivprobit} and {cmd:ivtobit} are fitted by joint maximum likelihood with one endogenous regressor, and report {cmd:/athrho}, {cmd:/lnsigma} and a Wald test of exogeneity; the {cmd:twostep} option is not available yet. {cmd:ivpoisson gmm} minimizes the GMM criterion in two steps (robust weight matrix, or one step with {cmd:onestep}); {cmd:cfunction} adds the first-stage residual to the exponential mean. {cmd:ivregress liml} uses the k-class estimator with the smallest eigenvalue as kappa; {cmd:ivregress gmm} is two-step efficient GMM.

{title:Postestimation}

{phang2}{cmd:predict} newvar [{cmd:,} {cmd:xb} {cmd:residuals}]{p_end}
{phang2}{cmd:test}, {cmd:lincom}, {cmd:estat ic}{p_end}

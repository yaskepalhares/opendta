{smcl}
{title:Title}

{p 4 4 2}
{bf:pca}, {bf:factor}, {bf:rotate}, {bf:rotatemat}, {bf:alpha}, {bf:screeplot}, {bf:loadingplot}, {bf:scoreplot} {hline 2} Principal components, factor analysis, rotation and scale reliability

{title:Syntax}

{p 8 16 2}
{cmd:pca} varlist [if] [in] [weight] [{cmd:,} {cmd:components(}#{cmd:)} {cmd:mineigen(}#{cmd:)} {cmd:correlation} {cmd:covariance} {cmd:vce(none}|{cmd:normal)} {cmd:level(}#{cmd:)} {cmd:novce} {cmd:means} {cmd:blanks(}#{cmd:)} {cmd:tol(}#{cmd:)} {cmd:ignore} {cmd:norotated}]

{p 8 16 2}
{cmd:factor} varlist [if] [in] [weight] [{cmd:,} {cmd:pf} {cmd:pcf} {cmd:ipf} {cmd:ml} {cmd:factors(}#{cmd:)} {cmd:mineigen(}#{cmd:)} {cmd:blanks(}#{cmd:)} {cmd:heywood} {cmd:iterate(}#{cmd:)} {cmd:tolerance(}#{cmd:)} {cmd:norotated}]

{p 8 16 2}
{cmd:rotate} [{cmd:,} {it:method} {cmd:oblique} {cmd:orthogonal} {cmd:normalize} {cmd:factors(}#{cmd:)} {cmd:components(}#{cmd:)} {cmd:blanks(}#{cmd:)} {it:optimize_options}]

{p 8 16 2}
{cmd:rotatemat} {it:matname} [{cmd:,} {it:method} {cmd:oblique} {cmd:orthogonal} {cmd:normalize} {cmd:format(}%fmt{cmd:)} {cmd:blanks(}#{cmd:)} {cmd:nodisplay} {cmd:noloading} {cmd:norotation} {cmd:matname(}string{cmd:)} {cmd:colnames(}string{cmd:)} {it:optimize_options}]

{p 8 16 2}
{it:method}: {cmd:varimax} {cmd:vgpf} {cmd:quartimax} {cmd:equamax} {cmd:parsimax} {cmd:entropy} {cmd:tandem1} {cmd:tandem2} {cmd:promax}[{cmd:(}#{cmd:)}] {cmd:oblimin}[{cmd:(}#{cmd:)}] {cmd:cf(}#{cmd:)} {cmd:bentler} {cmd:oblimax} {cmd:quartimin} {cmd:target(}Tg{cmd:)} {cmd:partial(}Tg W{cmd:)}; {it:optimize_options}: {cmd:iterate(}#{cmd:)} {cmd:tolerance(}#{cmd:)} {cmd:gtolerance(}#{cmd:)} {cmd:ltolerance(}#{cmd:)} {cmd:protect(}#{cmd:)} {cmd:maxstep(}#{cmd:)} {cmd:init(}matname{cmd:)} {cmd:random}

{p 8 16 2}
{cmd:alpha} varlist [if] [in] [{cmd:,} {cmd:asis} {cmd:casewise} {cmd:detail} {cmd:item} {cmd:label} {cmd:min(}#{cmd:)} {cmd:reverse(}varlist{cmd:)} {cmd:std} {cmd:generate(}newvar{cmd:)}]

{p 8 16 2}
{cmd:screeplot} [{cmd:,} {cmd:mean} {cmd:ci} {cmd:yline(}#{cmd:)}]   {cmd:loadingplot}   {cmd:scoreplot}

{title:Description}

{pstd}
{cmd:pca} computes the eigenvalues and eigenvectors of the correlation (default) or covariance matrix. {cmd:factor} fits principal factors ({cmd:pf}), principal-component factors ({cmd:pcf}), iterated principal factors ({cmd:ipf}) or maximum likelihood ({cmd:ml}, with the likelihood-ratio tests against independence and against the saturated model). {cmd:rotate} rotates the last {cmd:pca} or {cmd:factor} solution and {cmd:rotatemat} rotates a matrix; the criteria other than varimax and promax are minimized by the gradient projection algorithm. Rotations are orthogonal unless {cmd:oblique} is given or the method is {cmd:promax}. {cmd:pca, vce(normal)} adds the standard errors of the eigenvalues, eigenvectors and explained variance under multivariate normality and the likelihood-ratio tests of independence and sphericity. {cmd:alpha} reports Cronbach's alpha. The sign of each eigenvector is chosen so that its column sum is positive (VERIFY).

{title:Postestimation}

{phang2}{cmd:estat} {cmd:kmo} | {cmd:smc} | {cmd:anti} | {cmd:loadings} | {cmd:correlation} | {cmd:covariance} | {cmd:residuals} | {cmd:summarize} | {cmd:structure} | {cmd:common} | {cmd:factors} (after {cmd:factor}) | {cmd:rotatecompare}{p_end}
{phang2}{cmd:predict} newvarlist [{cmd:,} {cmd:regression} {cmd:bartlett} {cmd:fitted} {cmd:residuals} {cmd:comp(}numlist{cmd:)}]{p_end}

{title:Stored results}

{pstd}
{cmd:pca} and {cmd:factor} store {cmd:e(N)}, {cmd:e(f)} (number retained), {cmd:e(trace)}, the matrices {cmd:e(Ev)} and {cmd:e(L)} and the macros {cmd:e(cmd)}, {cmd:e(properties)} and {cmd:e(estat_cmd)}. {cmd:alpha} stores {cmd:r(alpha)}, {cmd:r(rho)}, {cmd:r(k)} and the average covariance and correlation.

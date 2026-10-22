{smcl}
{title:Title}

{p 4 4 2}
{bf:clogit}, {bf:rologit}, {bf:asclogit}, {bf:slogit}, {bf:mprobit}, {bf:asmprobit}, {bf:asroprobit} {hline 2} Choice models

{title:Syntax}

{p 8 16 2}
{cmd:clogit} depvar indepvars [if] [in] [weight] {cmd:, group(}varname{cmd:)} [{cmd:or} {cmd:vce(}vcetype{cmd:)}]

{p 8 16 2}
{cmd:rologit} depvar indepvars [if] [in] {cmd:, group(}varname{cmd:)} [{cmd:ties(}{cmd:exactm}|{cmd:breslow}|{cmd:efron}|{cmd:none}{cmd:)}]

{p 8 16 2}
{cmd:asclogit} depvar [indepvars] [if] [in] {cmd:, case(}varname{cmd:)} {cmd:alternatives(}varname{cmd:)} [{cmd:casevars(}varlist{cmd:)} {cmd:basealternative(}#{cmd:)} {cmd:noconstant}]

{p 8 16 2}
{cmd:slogit} depvar [indepvars] [if] [in] [{cmd:,} {cmd:dim(1)} {cmd:base(}#{cmd:)}]

{p 8 16 2}
{cmd:mprobit} depvar indepvars [if] [in] [{cmd:,} {cmd:baseoutcome(}#{cmd:)} {cmd:intpoints(}#{cmd:)}]

{p 8 16 2}
{cmd:asmprobit}|{cmd:asroprobit} depvar [indepvars] [if] [in] {cmd:, case(}varname{cmd:)} {cmd:alternatives(}varname{cmd:)} [{cmd:casevars(}varlist{cmd:)} {cmd:correlation(}{cmd:independent}|{cmd:unstructured}{cmd:)} {cmd:stddev(}{cmd:homoskedastic}|{cmd:heteroskedastic}{cmd:)} {cmd:intpoints(}#{cmd:)}]

{title:Description}

{pstd}
{cmd:clogit} drops the groups in which the outcome does not vary and accepts any number of chosen alternatives per group. {cmd:rologit} is the exploded logit for rankings. {cmd:mprobit}, {cmd:asmprobit} and {cmd:asroprobit} integrate the normal probabilities by the GHK simulator with Hammersley draws, so their results follow the draws and are compared with Stata only up to simulation error.

{pstd}
Not available: {cmd:slogit} with {cmd:dim()} above 1; {cmd:asmprobit} and {cmd:asroprobit} with {cmd:altwise}, unbalanced alternatives or other covariance structures.

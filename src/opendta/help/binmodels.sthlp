{smcl}
{title:Title}

{p 4 4 2}
{bf:cloglog}, {bf:scobit}, {bf:hetprobit}, {bf:fracreg}, {bf:glogit}, {bf:binreg}, {bf:zip}, {bf:zinb}, {bf:tpoisson}, {bf:tnbreg}, {bf:gnbreg} {hline 2} Binary, fractional and count-data models

{title:Syntax}

{p 8 16 2}
{cmd:cloglog} depvar [indepvars] [if] [in] [weight] [{cmd:,} {cmd:eform} {cmd:offset(}varname{cmd:)} {cmd:noconstant} {cmd:vce(}vcetype{cmd:)} {cmd:level(}#{cmd:)}]

{p 8 16 2}
{cmd:scobit} depvar [indepvars] [if] [in] [weight] [{cmd:,} {cmd:or} {cmd:offset(}varname{cmd:)} ...]

{p 8 16 2}
{cmd:hetprobit} depvar [indepvars] [if] [in] [weight] {cmd:,} {cmd:het(}varlist{cmd:)} [...]

{p 8 16 2}
{cmd:fracreg} {c -(}{cmd:logit}|{cmd:probit}|{cmd:cloglog}{c )-} depvar [indepvars] [if] [in] [weight] [{cmd:,} ...]

{p 8 16 2}
{cmd:binreg} depvar [indepvars] [if] [in] [weight] [{cmd:,} {cmd:or}|{cmd:rr}|{cmd:hr}|{cmd:rd} {cmd:n(}#|varname{cmd:)} ...]

{p 8 16 2}
{cmd:glogit} pos pop [indepvars] [if] [in] [{cmd:,} {cmd:or} ...]

{p 8 16 2}
{cmd:zip} depvar [indepvars] [if] [in] [weight] {cmd:,} {cmd:inflate(}varlist{cmd:)} [{cmd:probit} {cmd:vuong} {cmd:irr} ...]

{p 8 16 2}
{cmd:zinb} depvar [indepvars] [if] [in] [weight] {cmd:,} {cmd:inflate(}varlist{cmd:)} [{cmd:probit} {cmd:vuong} {cmd:irr} ...]

{p 8 16 2}
{cmd:tpoisson}|{cmd:tnbreg} depvar [indepvars] [if] [in] [weight] [{cmd:,} {cmd:ll(}#|varname{cmd:)} {cmd:ul(}#|varname{cmd:)} {cmd:irr} {cmd:dispersion(mean}|{cmd:constant)} ...]

{p 8 16 2}
{cmd:gnbreg} depvar [indepvars] [if] [in] [weight] [{cmd:,} {cmd:lnalpha(}varlist{cmd:)} {cmd:irr} ...]

{title:Description}

{pstd}
Maximum-likelihood estimators with Stata-style Newton iterations and numerical Hessians of analytic gradients. {cmd:binreg} is {cmd:glm} with a binomial family and the link chosen by the risk option; {cmd:glogit} is weighted least squares on the empirical logit. The layouts of the headers are marked VERIFY in the code until compared with logs.

{title:Postestimation}

{phang2}{cmd:predict} newvar [{cmd:,} {cmd:xb} {cmd:pr} {cmd:n} {cmd:stdp}]{p_end}
{phang2}{cmd:test}, {cmd:lincom}, {cmd:margins}, {cmd:estat ic}{p_end}

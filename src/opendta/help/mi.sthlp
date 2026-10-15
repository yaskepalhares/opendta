{smcl}
{title:Title}

{p 4 4 2}
{bf:mi} {hline 2} Multiple imputation

{title:Syntax}

{p 8 16 2}{cmd:mi set} {cmd:mlong}|{cmd:flong}{p_end}
{p 8 16 2}{cmd:mi register} {cmd:imputed}|{cmd:passive}|{cmd:regular} varlist; {cmd:mi unregister} varlist{p_end}
{p 8 16 2}{cmd:mi describe}; {cmd:mi query}; {cmd:mi misstable} ...{p_end}
{p 8 16 2}{cmd:mi impute} {cmd:regress}|{cmd:pmm}|{cmd:logit}|{cmd:ologit}|{cmd:mlogit}|{cmd:poisson} ivar [{cmd:=} indepvars]{cmd:,} {cmd:add(}#{cmd:)}|{cmd:replace} [{cmd:rseed(}#{cmd:)} {cmd:knn(}#{cmd:)}]{p_end}
{p 8 16 2}{cmd:mi impute chained} ({it:method}) varlist ({it:method}) varlist ... [{cmd:=} indepvars]{cmd:,} {cmd:add(}#{cmd:)}|{cmd:replace} [{cmd:burnin(}#{cmd:)} {cmd:rseed(}#{cmd:)}]{p_end}
{p 8 16 2}{cmd:mi estimate} [{cmd:,} {cmd:level(}#{cmd:)} {cmd:imputations(}numlist{cmd:)} {cmd:nosmall} {cmd:noisily}]{cmd::} estimation_command{p_end}
{p 8 16 2}{cmd:mi passive:} {cmd:generate}|{cmd:replace}|{cmd:egen} ...{p_end}
{p 8 16 2}{cmd:mi xeq} [numlist]{cmd::} command [{cmd:;} command ...]{p_end}
{p 8 16 2}{cmd:mi extract} #{cmd:, clear}{p_end}

{title:Description}

{pstd}
The imputations are stored in the mlong or flong style, with the system variables {cmd:_mi_m}, {cmd:_mi_id} and {cmd:_mi_miss}. {cmd:mi impute} draws the model parameters from their posterior and then the missing values; {cmd:chained} imputes several variables by fully conditional specification. {cmd:mi estimate} runs any estimation command on each imputation and combines the results with Rubin's rules (small-sample degrees of freedom when the command has residual degrees of freedom; model F test assuming equal fractions of missing information). The imputed values come from OpenDTA's random-number generator and differ from Stata's. The wide and flongsep styles, {cmd:mi impute mvn} and {cmd:monotone} are not available.

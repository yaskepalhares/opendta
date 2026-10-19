{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:unab} {hline 2}}Unabbreviate and expand varlists, for programmers{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:unab} {it:lmacname} {cmd::} {it:varlist} [{cmd:,} {it:options}]{break}
{cmd:tsunab} {it:lmacname} {cmd::} {it:varlist} [{cmd:,} {it:options}]{break}
{cmd:fvunab} {it:lmacname} {cmd::} {it:varlist} [{cmd:,} {it:options}]

{p 8 16 2}
{cmd:unabcmd} {it:name}

{p 8 16 2}
{cmd:fvexpand} {it:varlist} [{cmd:if}] [{cmd:in}]

{p 8 16 2}
{cmd:fvrevar} {it:varlist} [{cmd:if}] [{cmd:in}] [{cmd:,} {opt subst:itute} {opt tsonly}
{opt list} {opt stub(stub)}]

{p 8 16 2}
{cmd:tsrevar} {it:varlist} [{cmd:if}] [{cmd:in}] [{cmd:,} {opt subst:itute} | {opt list}]

{p 8 16 2}
{cmd:_rmcoll} {it:varlist} [{cmd:if}] [{cmd:in}] [{it:weight}] [{cmd:,} {opt nocons:tant}
{opt coll:inear} {opt exp:and} {opt forced:rop}]

{synoptset 16}{...}
{synopthdr:unab options}
{synoptline}
{synopt :{opt min(#)}}fewest variables allowed; default {cmd:1}{p_end}
{synopt :{opt max(#)}}most variables allowed; default {cmd:32000}{p_end}
{synopt :{opt name(string)}}label used in the error message{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:unab} puts the unabbreviated variable names (wildcards and ranges
expanded) in the local macro {it:lmacname}. {cmd:tsunab} also accepts
time-series operators ({cmd:L(1/3).mpg} becomes {cmd:L.mpg L2.mpg L3.mpg}) and
{cmd:fvunab} also factor variables. With {cmd:set varabbrev off}, an
abbreviation is an error.
{p_end}

{pstd}
{cmd:unabcmd} stores in {cmd:r(cmd)} the full name of a command, built in or a program
found on the adopath; if there is none, the error is {cmd:r(111)}.
{p_end}

{pstd}
{cmd:fvexpand} returns in {cmd:r(varlist)} the varlist with the factor variables expanded
to one term for each level present in the sample, the base level marked
({cmd:1b.rep78} {cmd:2.rep78} ...).
{p_end}

{pstd}
{cmd:fvrevar} and {cmd:tsrevar} return in {cmd:r(varlist)} a varlist without operators: a
temporary variable replaces each factor-variable term (one for each level, the base level
included), interaction or time-series-operated variable, and a variable that needs none
stays as it is. A temporary variable is used again when the same one is asked for. The options:
{cmd:tsonly} replaces only the time-series-operated variables and returns the other terms
expanded; {cmd:list} creates nothing and returns the base variable names, once each, in the
order of the dataset; {cmd:stub(}{it:stub}{cmd:)} (of {cmd:fvrevar}) creates {it:stub}{cmd:1},
{it:stub}{cmd:2}, ... instead of temporary variables.
{p_end}

{pstd}
{cmd:_rmcoll} flags the collinear variables of the list with {cmd:o.} ({cmd:o.tt}) and
writes a note for each one; a level that identifies no observation is noted as well.
{cmd:r(varlist)} is the list specific to the sample ({cmd:i(1 2 3 4 5)b1.rep78}); with {cmd:expand}
it has one term for each level; {cmd:r(k_omitted)} is the number of omitted terms.
{cmd:forcedrop} drops the collinear variables from the list instead of flagging them.
{p_end}

{title:Example}

{phang2}{cmd:. tsunab lags : L(1/2).(price turn)}{p_end}
{phang2}{cmd:. display "`lags'"}{p_end}
{phang2}{cmd:L.price L2.price L.turn L2.turn}{p_end}

{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:set seed} {hline 2}}Random numbers{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:set seed} {it:#}  |  {cmd:set seed} {it:statecode}{p_end}
{p 8 16 2}{cmd:set rng} {c -(}{cmd:default}|{cmd:mt64}{c )-}{p_end}

{title:Functions}

{p2colset 9 36 38 2}{...}
{p2col :{cmd:runiform()}, {cmd:runiform(}{it:a,b}{cmd:)}}uniform on (0,1) or ({it:a},{it:b}){p_end}
{p2col :{cmd:runiformint(}{it:a,b}{cmd:)}}integers from {it:a} to {it:b}{p_end}
{p2col :{cmd:rnormal()}, {cmd:rnormal(}{it:m,s}{cmd:)}}normal{p_end}
{p2col :{cmd:rbinomial(}{it:n,p}{cmd:)}, {cmd:rpoisson(}{it:m}{cmd:)}}{p_end}
{p2col :{cmd:rchi2(}{it:df}{cmd:)}, {cmd:rt(}{it:df}{cmd:)}, {cmd:rbeta(}{it:a,b}{cmd:)}, {cmd:rgamma(}{it:a,b}{cmd:)}}{p_end}
{p2col :{cmd:rexponential(}{it:b}{cmd:)}, {cmd:rlogistic()}, {cmd:rweibull(}{it:a,b}{cmd:)}}{p_end}
{p2col :{cmd:rnbinomial(}{it:n,p}{cmd:)}, {cmd:rhypergeometric(}{it:N,K,n}{cmd:)}}{p_end}
{p2colreset}{...}

{title:Description}

{pstd}
OpenDTA uses the 64-bit Mersenne Twister. The same seed always gives the
same numbers in OpenDTA, but the sequences are not guaranteed to match
other programs. {cmd:c(seed)} holds the current state; {cmd:set seed}
accepts it back, so a sequence can be resumed exactly.
{p_end}

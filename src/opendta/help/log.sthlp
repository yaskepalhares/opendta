{smcl}
{title:Title}

{p2colset 5 16 18 2}{...}
{p2col :{cmd:log} {hline 2}}Record the session in a file{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:log using} {it:filename} [{cmd:,} {cmd:append} {cmd:replace} {cmd:text} {cmd:smcl} {cmd:name(}{it:name}{cmd:)}]

{p 8 16 2}
{cmd:log} {c -(}{cmd:close}|{cmd:off}|{cmd:on}|{cmd:query}{c )-} [{it:name}|{cmd:_all}]

{p 8 16 2}
{cmd:cmdlog using} {it:filename} [{cmd:,} {cmd:append} {cmd:replace}]{space 4}{cmd:cmdlog close}

{title:Description}

{pstd}
{cmd:log} copies everything shown in Results to a file: plain text
({cmd:.log}) or SMCL ({cmd:.smcl}, the default), which keeps the colors and
opens in the Viewer. Several logs can be open at once with {cmd:name()}.
{cmd:cmdlog} records only the commands that were typed.
{p_end}

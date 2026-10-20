{smcl}
{title:Title}

{p 4 4 2}
{bf:cumul}, {bf:stem}, {bf:lv}, {bf:smooth} {hline 2} Cumulative distribution, stem-and-leaf plot, letter values and resistant smoothers

{title:Syntax}

{p 8 16 2}
{cmd:cumul} {it:varname} [{it:weight}] {cmd:,} {cmd:generate(}{it:newvar}{cmd:)} [{cmd:freq} {cmd:equal}]

{p 8 16 2}
{cmd:stem} {it:varname} [{cmd:,} {cmd:lines(}#{cmd:)} {cmd:round(}#{cmd:)} {cmd:digits(}#{cmd:)} {cmd:prune}]

{p 8 16 2}
{cmd:lv} [{it:varlist}]

{p 8 16 2}
{cmd:smooth} {it:smoother}[{cmd:, twice}] {it:varname} {cmd:,} {cmd:generate(}{it:newvar}{cmd:)}

{title:Description}

{pstd}
{cmd:cumul} creates the empirical cumulative distribution function. {cmd:stem} displays a stem-and-leaf plot. {cmd:lv} displays letter values, the spreads, pseudosigmas and the inner and outer fences. {cmd:smooth} applies Tukey's smoothers: digits (running medians; an even span is followed by 2), {cmd:R} (repeat), {cmd:S} (split), {cmd:E} (end-point rule) and {cmd:H} (hanning).

{title:Limitations}

{pstd}
{cmd:stem} does not take {cmd:width()}; {cmd:lv} does not take {cmd:generate()} or {cmd:tail()}. Details that the manual does not show are marked VERIFY.

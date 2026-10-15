{smcl}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:set hints} {hline 2}}Explain errors{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:set hints} {c -(}{cmd:on}|{cmd:off}{c )-} [{cmd:,} {cmdab:perm:anently}]

{title:Description}

{pstd}
When an error stops a command, OpenDTA shows the usual error message and
return code {cmd:r(#)}. With {cmd:set hints on} (the default), an explanation
written for people learning Stata comes in between:
{p_end}

{p 8 12 2}
{bf:what went wrong} in this command, with the names and numbers of the data
in memory: which variable is missing in how many observations, which variable
is text and which is numeric, which values the outcome takes, similar variable,
option, command and file names;
{p_end}
{p 8 12 2}
{bf:why}: the idea behind the rule (= and ==, missing values, generate and
egen, sorting for {cmd:by}...);
{p_end}
{p 8 12 2}
{bf:How to fix}: the corrected command, when it can be written, or what to do;
{p_end}
{p 8 12 2}
{bf:Learn more}: a {cmd:help} topic. In the Results window it is a link: click
it to open the page.
{p_end}

{pstd}
The explanation follows the interface language ({help language:set language}).
{p_end}

{pstd}
The return codes never change, so do-files that test {cmd:_rc} behave the
same with hints on or off. Errors raised on purpose with {cmd:error #}
are not explained.
{p_end}

{pstd}
{cmd:permanently} also remembers the choice the next time OpenDTA opens.
The same switch is in the menu {bf:Help > Explain Errors}.
{p_end}

{title:Example}

        {cmd:. list if sex = 1}
        {err:invalid syntax}
          {txt:-> The condition after "if" (sex = 1) uses a single "=".}
            {txt:In Stata, a single "=" gives a value (as in "generate y = 1"); to ask}
            {txt:whether two things are equal, use a double "==". ...}
            {txt:How to fix: list if sex == 1}
            {txt:Learn more: help operators}
        {err:r(198);}

{title:Also see}

{pstd}
{help operators}, {help missing}, {help functions}, {help by}
{p_end}

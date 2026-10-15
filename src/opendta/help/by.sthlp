{smcl}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:by} {hline 2}}Repeat a command for each group{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:by} {varlist}{cmd::} {it:command}

{p 8 16 2}
{cmd:bysort} {varlist}{cmd::} {it:command}

{title:Description}

{pstd}
{cmd:by} runs the command once for each group of observations with the same
values of the variables in {varlist}. Inside each group, {cmd:_n} is the number
of the observation in the group and {cmd:_N} the size of the group.
{p_end}

{pstd}
{cmd:by} needs the data sorted by those variables, because a group is a run of
consecutive observations. {cmd:bysort} (or {cmd:by ..., sort}) sorts first and
then runs the command.
{p_end}

{title:Examples}

        {cmd:. bysort city: summarize income}
        {cmd:. bysort family: generate size = _N}
        {cmd:. bysort family (age): generate oldest = (_n == _N)}
        {cmd:. bysort id: egen total = total(spending)}

{pstd}
In the third example, the variable in parentheses (age) sorts the observations
within each family without defining new groups.
{p_end}

{title:Also see}

{pstd}
{help datamgmt:sort}, {help egen}, {help collapse}
{p_end}

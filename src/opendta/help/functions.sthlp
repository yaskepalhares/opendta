{smcl}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:functions} {hline 2}}Functions for expressions{p_end}
{p2colreset}{...}

{title:Description}

{pstd}
Functions are used in expressions ({cmd:generate}, {cmd:replace}, {cmd:if},
{cmd:display}...). Their names are lowercase and are followed directly by
parentheses. They compute one value per observation, from that observation's
values; statistics over many observations (means, totals, ranks) are made with
{help egen}.
{p_end}

{title:Mathematical}

{p2colset 5 28 30 2}{...}
{p2col :{cmd:abs(x)} {cmd:sqrt(x)}}absolute value, square root{p_end}
{p2col :{cmd:exp(x)} {cmd:ln(x)} {cmd:log10(x)}}exponential and logarithms{p_end}
{p2col :{cmd:round(x, 0.1)}}round to a multiple ({cmd:floor()}, {cmd:ceil()}, {cmd:int()}){p_end}
{p2col :{cmd:mod(x, y)}}remainder of x / y{p_end}
{p2col :{cmd:min(a, b, ...)} {cmd:max()}}smallest and largest of the arguments{p_end}
{p2colreset}{...}

{title:Missing values and choices}

{p2colset 5 28 30 2}{...}
{p2col :{cmd:missing(x)}}1 if x is missing ({cmd:!missing(x)}: not missing){p_end}
{p2col :{cmd:cond(c, a, b)}}a when the condition c is true, b otherwise{p_end}
{p2col :{cmd:inlist(x, 1, 3, 5)}}1 if x is one of the values{p_end}
{p2col :{cmd:inrange(x, 18, 65)}}1 if 18 <= x <= 65{p_end}
{p2colreset}{...}

{title:Text (string)}

{p2colset 5 28 30 2}{...}
{p2col :{cmd:strlen(s)}}number of characters{p_end}
{p2col :{cmd:substr(s, 1, 3)}}part of the text: from position 1, 3 characters{p_end}
{p2col :{cmd:strupper(s)} {cmd:strlower(s)}}upper and lower case{p_end}
{p2col :{cmd:strtrim(s)}}remove spaces at the ends{p_end}
{p2col :{cmd:strpos(s, "x")}}position of "x" in s (0 if absent){p_end}
{p2col :{cmd:subinstr(s, "a", "b", .)}}replace "a" by "b"{p_end}
{p2col :{cmd:real(s)} {cmd:string(x)}}text to number and number to text{p_end}
{p2colreset}{...}

{title:Statistical}

{p2colset 5 28 30 2}{...}
{p2col :{cmd:normal(z)} {cmd:invnormal(p)}}standard normal distribution{p_end}
{p2col :{cmd:ttail(df, t)} {cmd:chi2tail(df, x)}}tail probabilities{p_end}
{p2col :{cmd:runiform()} {cmd:rnormal()}}random numbers (see {help set_seed:set seed}){p_end}
{p2colreset}{...}

{title:Dates}

{p2colset 5 28 30 2}{...}
{p2col :{cmd:mdy(m, d, y)}}date from month, day and year{p_end}
{p2col :{cmd:date(s, "DMY")}}date from text{p_end}
{p2col :{cmd:year(d)} {cmd:month(d)} {cmd:day(d)}}parts of a date{p_end}
{p2colreset}{...}

{title:Examples}

        {cmd:. generate lninc = ln(income)}
        {cmd:. generate first3 = substr(name, 1, 3)}
        {cmd:. generate group = cond(age < 18, 1, 2) if !missing(age)}
        {cmd:. list if inlist(state, "MG", "SP")}

{title:Also see}

{pstd}
{help operators}, {help egen}, {help generate}
{p_end}

{smcl}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:missing values} {hline 2}}How missing data are stored and used{p_end}
{p2colreset}{...}

{title:Description}

{pstd}
A numeric variable that has no value in an observation holds a missing value,
shown as a dot: {cmd:.}. There are also extended missing values, {cmd:.a} to
{cmd:.z}, to record why a value is missing (refused, not applicable...). A
missing string is the empty text {cmd:""}.
{p_end}

{title:Remarks}

{pstd}
{bf:Missing is larger than any number.} In comparisons, {cmd:.} counts as larger
than every number, and {cmd:.a} < {cmd:.b} < ... < {cmd:.z}. So
{cmd:if age > 60} also keeps the observations with age missing; write
{cmd:if age > 60 & !missing(age)}.
{p_end}

{pstd}
{bf:Calculations.} Any calculation with a missing value gives a missing value:
if x is missing, {cmd:x + 1} is missing too.
{p_end}

{pstd}
{bf:Statistics and models.} {cmd:summarize}, {cmd:regress}, {cmd:logit} and the
other estimation commands use only the observations with a value in every
variable of the command (complete cases). A variable with many missing values
shrinks the sample; when no observation is complete, the command stops with
{err:no observations}.
{p_end}

{title:Finding missing values}

{p2colset 5 36 38 2}{...}
{p2col :{cmd:misstable summarize}}how many values are missing in each variable{p_end}
{p2col :{cmd:count if missing(x)}}how many observations lack x{p_end}
{p2col :{cmd:list id x if missing(x)}}which observations lack x{p_end}
{p2col :{cmd:tabulate x, missing}}frequencies including the missing ones{p_end}
{p2colreset}{...}

{title:Examples}

        {cmd:. misstable summarize income age schooling}
        {cmd:. generate rich = (income > 10000) if !missing(income)}
        {cmd:. replace income = . if income == -99}
        {cmd:. mvdecode income, mv(-99 = .a)}

{title:Also see}

{pstd}
{help operators}, {help codebook} ({cmd:misstable}), {help recode} ({cmd:mvdecode})
{p_end}

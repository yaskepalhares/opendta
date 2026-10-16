{smcl}
{title:Title}

{p 4 4 2}
{bf:codebook} {hline 2} Describe data contents{break}
{bf:misstable} {hline 2} Tabulate missing values

{title:Syntax}

{p 8 16 2}
{cmd:codebook} [{it:varlist}] [{it:if}] [{it:in}] [{cmd:,} {cmd:compact} {cmd:problems}
{cmd:header} {cmd:tabulate(}{it:#}{cmd:)}]

{p 8 16 2}
{cmd:misstable} {cmd:summarize} [{it:varlist}] [{it:if}] [{it:in}] [{cmd:,} {cmd:all} {cmd:showzeros}]

{p 8 16 2}
{cmd:misstable} {cmd:patterns} [{it:varlist}] [{it:if}] [{it:in}] [{cmd:,} {cmd:frequency}]

{p 8 16 2}
{cmd:misstable} {cmd:tree} [{it:varlist}] [{it:if}] [{it:in}] [{cmd:,} {cmd:frequency}]

{p 8 16 2}
{cmd:misstable} {cmd:nested} [{it:varlist}] [{it:if}] [{it:in}]

{title:Description}

{pstd}
{cmd:codebook} shows, for each variable, its type, value label, range,
units, number of unique values and of missing values; variables with at most
{cmd:tabulate()} unique values (default 9) get a frequency table, the others
the mean, standard deviation and percentiles. String variables show examples
and a warning about embedded blanks. {cmd:compact} gives one line per
variable; {cmd:problems} lists constant variables and variables with
embedded blanks.

{pstd}
{cmd:misstable summarize} counts, for each numeric variable with missing
values, the observations equal to {cmd:.}, greater than {cmd:.} (extended
missing) and nonmissing. {cmd:misstable patterns} shows the patterns of
missing values (1 means complete). {cmd:misstable tree} shows the nested
pattern: the variable with the most missing values splits the observations
into missing and not missing, the next one splits each part, and so on
(percentages, or counts with {cmd:frequency}). {cmd:misstable nested} lists
chains such as {cmd:a(10) -> b(20)}: whenever {cmd:a} is missing, {cmd:b} is
missing too.

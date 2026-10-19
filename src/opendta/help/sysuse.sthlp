{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:sysuse} {hline 2}}Example datasets{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:sysuse} {it:filename} [{cmd:,} {opt clear}]{space 4}{cmd:sysuse dir} [{cmd:,} {opt all}]{p_end}
{p 8 16 2}{cmd:webuse} {it:filename} [{cmd:,} {opt clear}]{space 4}{cmd:webuse query}{space 4}{cmd:webuse set} [{it:url}]{p_end}

{title:Description}

{pstd}
{cmd:sysuse} loads one of OpenDTA's own example datasets ({cmd:sysuse dir} lists them). They are synthetic data built by {cmd:tools/make_datasets.py}; Stata's example datasets are StataCorp files and are not included.
{p_end}

{pstd}
{cmd:webuse} downloads a dataset from a web address: by default the public address of Stata's example datasets, so that the examples of the manuals can be followed; the file is downloaded by you when you ask for it and is not part of OpenDTA. {cmd:webuse set} changes the address and {cmd:webuse query} shows it.
{p_end}

{title:Examples}

{phang2}{cmd:. sysuse dir}{p_end}
{phang2}{cmd:. sysuse diet, clear}{p_end}
{phang2}{cmd:. webuse auto, clear}{p_end}

{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:save} {hline 2}}Save the data in memory{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:save} [{it:filename}] [{cmd:,} {cmd:replace} {cmd:emptyok}]

{p 8 16 2}
{cmd:saveold} {it:filename} [{cmd:,} {cmd:replace} {cmd:version(}{it:#}{cmd:)}]

{title:Description}

{pstd}
{cmd:save} writes the data, labels, notes and characteristics to a
{cmd:.dta} file (format 118, or 119 with more than 32,767 variables).
Without a file name, the data go back to the file they came from.
The file is written to a temporary name first, so a failure never
destroys the previous version.
{p_end}

{pstd}
{cmd:saveold} writes formats read by older programs: {cmd:version(13)}
(the default), {cmd:version(12)} or {cmd:version(11)}. Versions 11 and 12
cannot hold long strings; they are cut to 244 characters with a note.
{p_end}

{title:Options}

{phang}
{cmd:replace} overwrites an existing file.

{phang}
{cmd:emptyok} allows saving a dataset with no variables.

{title:Also see}

{pstd}
{help use}, {help import}
{p_end}

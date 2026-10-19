{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:filesys} {hline 2}}Files, folders and the operating system{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:copy} {it:from} {it:to} [{cmd:,} {opt text} {opt pub:lic} {opt replace}]{p_end}
{p 8 16 2}{cmd:mkdir} {it:directory} [{cmd:,} {opt pub:lic}]{p_end}
{p 8 16 2}{cmd:rmdir} {it:directory}{p_end}
{p 8 16 2}{cmd:dir} [{it:filespec}] [{cmd:,} {opt w:ide}]{space 4}({cmd:ls} is a synonym){p_end}
{p 8 16 2}{cmd:shell} [{it:operating_system_command}]{space 4}({cmd:!} is a synonym){p_end}
{p 8 16 2}{cmd:filefilter} {it:oldfile} {it:newfile} {cmd:,} {opt from(oldpattern)} {opt to(newpattern)} [{opt replace}]{p_end}
{p 8 16 2}{cmd:zipfile} {it:file|directory} [...] {cmd:,} {opt saving(zipfile [, replace])}{p_end}
{p 8 16 2}{cmd:unzipfile} {it:zipfile} [{cmd:,} {opt replace}]{p_end}
{p 8 16 2}{cmd:checksum} {it:filename} [{cmd:,} {opt save} {opt saving(filename)} {opt replace}]{p_end}

{title:Description}

{pstd}
{cmd:copy} copies a file (also from a URL); {cmd:mkdir} and {cmd:rmdir} create and remove a folder (the folder must be empty); {cmd:dir} lists files; {cmd:shell} runs a command of the operating system and shows its output; {cmd:filefilter} copies a text file replacing a pattern (with the codes {cmd:\n}, {cmd:\r}, {cmd:\t} and {cmd:\BS} for backslash); {cmd:zipfile} and {cmd:unzipfile} create and extract zip archives; {cmd:checksum} reports a file's checksum and size.
{p_end}

{title:Examples}

{phang2}{cmd:. mkdir results}{p_end}
{phang2}{cmd:. copy mydata.dta results/backup.dta, replace}{p_end}
{phang2}{cmd:. dir *.dta}{p_end}
{phang2}{cmd:. filefilter old.csv new.csv, from(";") to(",") replace}{p_end}
{phang2}{cmd:. zipfile results, saving(results.zip, replace)}{p_end}

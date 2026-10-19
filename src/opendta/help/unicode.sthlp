{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:unicode} {hline 2}}Converting old files to Unicode{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:unicode analyze} {it:filespec}{p_end}
{p 8 16 2}{cmd:unicode encoding set} {it:encoding}{p_end}
{p 8 16 2}{cmd:unicode translate} {it:filespec} [{cmd:,} {opt redo} {opt invalid} {opt transutf8} {opt nodata}]{p_end}
{p 8 16 2}{cmd:unicode retranslate} {it:filespec} [{cmd:,} {opt invalid} {opt transutf8} {opt replace}]{p_end}
{p 8 16 2}{cmd:unicode restore} {it:filespec}{p_end}
{p 8 16 2}{cmd:unicode convertfile} {it:srcfile} {it:dstfile} [{cmd:,} {opt srcenc:oding(enc)} {opt dstenc:oding(enc)} {opt replace}]{p_end}

{title:Description}

{pstd}
{cmd:unicode} converts the text of datasets and do-files written by older versions (in an extended ASCII encoding such as Latin-1) to UTF-8: {cmd:analyze} reports which files need it, {cmd:encoding set} chooses the old encoding, {cmd:translate} converts the files and keeps backups, which {cmd:restore} brings back; {cmd:convertfile} converts any text file between two encodings.
{p_end}

{title:Examples}

{phang2}{cmd:. unicode analyze *.dta}{p_end}
{phang2}{cmd:. unicode encoding set latin1}{p_end}
{phang2}{cmd:. unicode translate *.dta}{p_end}

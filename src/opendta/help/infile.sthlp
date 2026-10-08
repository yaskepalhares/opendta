{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:infile} {hline 2}}Read text data in free format or fixed columns{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:infile} [{it:type}] {it:newvar} [{it:newvar} ...] {cmd:using} {it:filename}
[{it:if}] [{it:in}] [{cmd:,} {cmd:clear} {cmdab:a:utomatic}]

{p 8 16 2}
{cmd:infile using} {it:dictionary} [{cmd:,} {cmd:using(}{it:datafile}{cmd:)} {cmd:clear}]

{p 8 16 2}
{cmd:infix} [{it:#} {cmd:lines}] [{it:#}{cmd::}] [{it:type}] {it:newvar} {it:start}{cmd:-}{it:end} ...
{cmd:using} {it:filename} [{cmd:,} {cmd:clear}]

{title:Description}

{pstd}
In free format, values are separated by spaces or commas and read in
order regardless of line breaks; string variables must be declared, as
in {cmd:str20 nome}. A type applies only to the name that follows it;
{cmd:int(a b c)} applies it to several.
{p_end}

{pstd}
Fixed-column files, such as census microdata, are read with {cmd:infix}
or with a dictionary file ({cmd:.dct}) that lists {cmd:_column(}{it:#}{cmd:)},
the type, the name, a format such as {cmd:%5.2f} (two implied decimals)
and an optional label for each variable.
{p_end}

{title:Example dictionary}

        dictionary using pessoas.raw {
            _column(1)  long  id     %6f   "Identifier"
                        str10 nome   %10s  "Name"
                        float renda  %7.2f "Income"
        }

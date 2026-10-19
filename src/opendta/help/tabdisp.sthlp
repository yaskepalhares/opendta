{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:tabdisp} {hline 2}}Display tables{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:tabdisp} {it:rowvar} [{it:colvar} [{it:supercolvar}]] [{cmd:if}] [{cmd:in}]{cmd:,}
{opth c:ellvar(varnames)} [{it:options}]

{synoptset 24}{...}
{synopthdr}
{synoptline}
{synopt :{opth by(varlist)}}superrow variables (up to four){p_end}
{synopt :{opt format(%fmt)}}format of the numbers; default {cmd:%9.0g}{p_end}
{synopt :{opt center}}center the contents of the cells{p_end}
{synopt :{opt left}}left-align the column labels{p_end}
{synopt :{opt concise}}omit rows whose entries are all missing{p_end}
{synopt :{opt missing}}show the missing-value code in the cells that have one{p_end}
{synopt :{opt totals}}label a missing value of a classifying variable as {cmd:Total}{p_end}
{synopt :{opt dotz}}{cmd:.z} instead of {cmd:.} marks the totals{p_end}
{synopt :{opt cellwidth(#)}}{opt csepwidth(#)} {opt scsepwidth(#)} {opt stubwidth(#)}: widths{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:tabdisp} displays the data in memory as a table, like {cmd:list} but with the
classifying variables as the headers: it calculates nothing, and a cell is the value of the
cell variables in the first observation with that row, column, supercolumn and superrows.
Value labels and variable labels are used. With no {it:colvar}, the cell variables are the columns.
{p_end}

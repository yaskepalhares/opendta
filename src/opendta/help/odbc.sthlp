{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:odbc} {hline 2}}Reading and writing databases through ODBC{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:odbc list}{p_end}
{p 8 16 2}{cmd:odbc query} [{it:"DataSourceName"}] [{cmd:,} {it:connect_options}]{p_end}
{p 8 16 2}{cmd:odbc describe} [{it:"TableName"}] [{cmd:,} {it:connect_options}]{p_end}
{p 8 16 2}{cmd:odbc load} [{it:extvarlist}] [{cmd:,} {opt table("TableName")} | {opt exec("SqlStmt")} {opt clear} {it:connect_options}]{p_end}
{p 8 16 2}{cmd:odbc insert} [{it:varlist}] {cmd:,} {opt table("TableName")} [{opt create} {it:connect_options}]{p_end}
{p 8 16 2}{cmd:odbc exec(}{it:"SqlStmt"}{cmd:)} [{cmd:,} {it:connect_options}]{p_end}

{p 8 16 2}{it:connect_options}: {opt dsn("DataSourceName")} | {opt conn:ectionstring("ConnectStr")} [{opt u:ser(UserID)} {opt p:assword(Password)}]{p_end}

{title:Description}

{pstd}
{cmd:odbc} lists the data sources, describes their tables, loads a table or the result of a query into memory, inserts the data in memory into a table and runs SQL statements.
{p_end}

{title:Remarks}

{pstd}
{cmd:odbc} needs the Python package {cmd:pyodbc} ({cmd:pip install pyodbc}) and an ODBC driver for the database on your computer; without them it explains what is missing.
{p_end}

{title:Examples}

{phang2}{cmd:. odbc list}{p_end}
{phang2}{cmd:. odbc load, table("patients") dsn("clinic") clear}{p_end}

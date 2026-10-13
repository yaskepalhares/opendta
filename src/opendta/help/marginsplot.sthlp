{smcl}
{title:Title}

{p 4 4 2}
{bf:marginsplot} {hline 2} Graph results from margins

{title:Syntax}

{p 8 16 2}
{cmd:marginsplot} [, noci twoway_options]

{title:Description}

{pstd}
Plots the margins and confidence intervals of the last {cmd:margins}. With {cmd:at()} varying one variable, it goes on the x axis; otherwise the levels of the factors. Interactions such as {cmd:a#b} give one line per level of {cmd:b}.

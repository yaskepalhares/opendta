{smcl}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:operators} {hline 2}}Arithmetic, comparison and logical operators{p_end}
{p2colreset}{...}

{title:Description}

{pstd}
Expressions, such as the right side of {cmd:generate y = ...} or the condition
after {cmd:if}, combine numbers, variables and functions with operators.
{p_end}

{title:Operators}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:+ - * /}}add, subtract, multiply, divide ({cmd:+} also joins text: {cmd:"a" + "b"}){p_end}
{p2col :{cmd:^}}power: {cmd:x^2}{p_end}
{p2col :{cmd:==}}equal to (two signs: one {cmd:=} gives a value){p_end}
{p2col :{cmd:!=} or {cmd:~=}}different from{p_end}
{p2col :{cmd:< > <= >=}}less than, greater than, less or equal, greater or equal{p_end}
{p2col :{cmd:&}}and: both conditions are true{p_end}
{p2col :{cmd:|}}or: at least one condition is true{p_end}
{p2col :{cmd:!}}not: {cmd:!missing(x)}{p_end}
{p2colreset}{...}

{pstd}
A comparison is 1 when true and 0 when false, so {cmd:generate adult = (age >= 18)}
creates a 0/1 variable. Parentheses set the order: {cmd:(a + b) / 2}.
{p_end}

{title:Remarks}

{pstd}
{bf:= and ==.} A single {cmd:=} gives a value ({cmd:generate y = 1}); a double
{cmd:==} asks whether two values are equal ({cmd:list if y == 1}).
{p_end}

{pstd}
{bf:Text.} Text values go in double quotes and are compared with text:
{cmd:if city == "Belo Horizonte"}. Upper case, accents and spaces count.
{p_end}

{pstd}
{bf:Missing values.} A missing value ({cmd:.}) is larger than any number, so
{cmd:if income > 1000} also selects the observations where income is missing.
Write {cmd:if income > 1000 & !missing(income)}; see {help missing}.
{p_end}

{title:Examples}

        {cmd:. generate bmi = weight / (height/100)^2}
        {cmd:. list name if age >= 18 & sex == 1}
        {cmd:. count if city != "Belo Horizonte"}
        {cmd:. generate obese = (bmi >= 30) if !missing(bmi)}

{title:Also see}

{pstd}
{help missing}, {help functions}, {help generate}
{p_end}

{smcl}
{title:Title}

{p 4 4 2}
{bf:stset} {hline 2} Declare data to be survival-time data

{title:Syntax}

{p 8 16 2}
{cmd:stset} timevar [weight] [, {cmd:failure(}failvar[{cmd:==}numlist]{cmd:)} {cmd:id(}idvar{cmd:)} {cmd:time0(}varname{cmd:)} {cmd:origin(time} exp{cmd:)} {cmd:enter(time} exp{cmd:)} {cmd:exit(failure}|{cmd:time} exp{cmd:)} {cmd:scale(}#{cmd:)} {cmd:if(}exp{cmd:)} {cmd:noshow}]

{p 8 16 2}
{cmd:stset}{space 3}(redisplay and re-run the current settings)

{p 8 16 2}
{cmd:stset, clear}

{p 8 16 2}
{cmd:st}{space 6}(show the settings)

{p 8 16 2}
{cmd:streset} [, options]{space 3}(stset again, changing some options)

{title:Description}

{pstd}
{cmd:stset} declares the data to be survival-time data. It creates {cmd:_st} (1 for observations used), {cmd:_d} (1 for failure), {cmd:_t} (analysis time at the end of the record) and {cmd:_t0} (analysis time at its beginning), and reports the observations excluded and why. With {cmd:id()}, records of the same subject are spells; without {cmd:time0()}, each spell begins where the previous one ended. Analysis time is (timevar - origin)/scale.

{pstd}
The settings are kept in the data characteristics and used by {cmd:stdescribe}, {cmd:stsum}, {cmd:sts}, {cmd:stcox} and {cmd:streg}.

{title:Example}

{phang2}{cmd:. stset studytime, failure(died)}{p_end}
{phang2}{cmd:. stset t, id(patid) failure(event==1) origin(time dob) scale(365.25)}{p_end}

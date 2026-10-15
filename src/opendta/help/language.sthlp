{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:set language} {hline 2}}Language of the interface (OpenDTA addition){p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:set language} {cmd:en}|{cmd:pt} [{cmd:,} {opt perm:anently}]{p_end}
{p 8 16 2}{cmd:set messages} {cmd:stata}|{cmd:translated} [{cmd:,} {opt perm:anently}]{p_end}

{title:Description}

{pstd}
OpenDTA addition, not part of Stata. {cmd:set language pt} shows menus, windows, dialogs and the error explanations ({help hints}) in Portuguese; {cmd:set language en} goes back to English, the default.
{p_end}

{pstd}
Messages in Stata's format (error messages such as {cmd:variable x not found}) and all command output stay in English, so that results can be compared with Stata logs and manuals. {cmd:set messages translated} translates the error messages too (only with {cmd:set language pt}); return codes do not change.
{p_end}

{pstd}
{bf:Settings} (Ctrl+, or Cmd+,) has the same choices, saved for the next session; so does {opt permanently}.
{p_end}

{title:Examples}

{phang2}{cmd:. set language pt}{p_end}
{phang2}{cmd:. set messages translated, permanently}{p_end}

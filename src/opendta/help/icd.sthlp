{smcl}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:icd} {hline 2}}ICD-9 and ICD-10 diagnosis codes{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:icd9}|{cmd:icd10} {cmd:check} {it:varname} [{it:if}] [{it:in}] [{cmd:,} {opt any} {opt l:ist} {opt gen:erate(newvar)} {opt sum:mary}]{p_end}
{p 8 16 2}{cmd:icd9}|{cmd:icd10} {cmd:clean} {it:varname} [{it:if}] [{it:in}] [{cmd:,} {opt dots} {opt gen:erate(newvar)}]{p_end}
{p 8 16 2}{cmd:icd9}|{cmd:icd10} {cmd:generate} {it:newvar} {cmd:=} {it:varname} {cmd:,} {opt cat:egory} | {opt sh:ort} | {opt range(codelist)} | {opt desc:ription} [{opt l:ong} {opt add:code(begin|end)} {opt nodots}]{p_end}
{p 8 16 2}{cmd:icd9}|{cmd:icd10} {cmd:lookup} {it:codelist}{p_end}
{p 8 16 2}{cmd:icd9}|{cmd:icd10} {cmd:search} {it:"text"} [...] [{cmd:,} {opt or}]{p_end}
{p 8 16 2}{cmd:icd9}|{cmd:icd10} {cmd:query}{p_end}
{p 8 16 2}{cmd:icd9p} {cmd:check}|{cmd:clean}|{cmd:generate} ...{space 4}(ICD-9 procedure codes){p_end}

{title:Description}

{pstd}
{cmd:check} reports codes that are badly formed or not defined in the code list; {cmd:clean} standardizes them (upper case, without spaces; {opt dots} adds the period); {cmd:generate} creates the category, the first four characters, an indicator for code ranges, or the description; {cmd:lookup} lists codes and descriptions (ranges {it:a}{cmd:/}{it:b} and wildcards {cmd:*}); {cmd:search} finds codes by the words in their descriptions.
{p_end}

{title:Remarks}

{pstd}
OpenDTA uses the public-domain code lists of the U.S. Centers for Medicare & Medicaid Services: ICD-9-CM diagnosis codes (version 29) for {cmd:icd9} and ICD-10-CM (2018) for {cmd:icd10}. ICD-10-CM extends the WHO ICD-10 with codes of up to seven characters, so a few codes differ from the WHO list. {cmd:icd9p} checks procedure codes by their format only, since the procedure list is not included. {opt any} checks the format only.
{p_end}

{title:Examples}

{phang2}{cmd:. icd10 check diag, generate(problem)}{p_end}
{phang2}{cmd:. icd10 generate desc = diag, description}{p_end}
{phang2}{cmd:. icd10 search "type 2 diabetes"}{p_end}
{phang2}{cmd:. icd9 lookup 250/250.93}{p_end}

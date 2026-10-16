{smcl}
{title:OpenDTA help}

{pstd}
OpenDTA runs analysis scripts (do-files), reads and writes {cmd:.dta}
datasets and keeps every change as a command. Type {cmd:help} followed by
a command name, or click a topic below.
{p_end}

{title:Basics}

{p2colset 5 24 26 2}{...}
{p2col :{help operators}}=, ==, &, | and the other operators{p_end}
{p2col :{help functions}}functions for expressions (ln, substr, cond...){p_end}
{p2col :{help missing}}missing values and how commands use them{p_end}
{p2col :{help by}}repeat a command for each group ({cmd:bysort}){p_end}
{p2colreset}{...}

{title:Data files}

{p2colset 5 24 26 2}{...}
{p2col :{help use}}load a .dta dataset{p_end}
{p2col :{help save}}save the data in memory (also {cmd:saveold}){p_end}
{p2col :{help import}}import delimited text and Excel files{p_end}
{p2col :{help infile}}read free-format and fixed-column text ({cmd:infile}, {cmd:infix}){p_end}
{p2colreset}{...}

{title:Working with data}

{p2colset 5 24 26 2}{...}
{p2col :{help describe}}describe the dataset{p_end}
{p2col :{help list}}list values{p_end}
{p2col :{help generate}}create and change variables ({cmd:generate}, {cmd:replace}){p_end}
{p2col :{help preserve}}save a copy of the data and bring it back{p_end}
{p2col :{help assert}}check that a condition holds{p_end}
{p2colreset}{...}

{title:Statistics and data management}

{p2colset 5 24 26 2}{...}
{p2col :{help summarize}}summary statistics ({cmd:tabstat}, {cmd:centile}, {cmd:correlate}...){p_end}
{p2col :{help ttest}}tests for means and proportions ({cmd:prtest}, {cmd:ci}){p_end}
{p2col :{help tabulate}}frequency tables ({cmd:tab1}, {cmd:tab2}, {cmd:table}){p_end}
{p2col :{help egen}}extensions to generate{p_end}
{p2col :{help collapse}}group statistics ({cmd:contract}, {cmd:expand}, {cmd:fillin}){p_end}
{p2col :{help merge}}combine datasets ({cmd:append}, {cmd:joinby}, {cmd:cross}){p_end}
{p2col :{help reshape}}wide and long forms{p_end}
{p2col :{help recode}}recode and convert ({cmd:encode}, {cmd:destring}, {cmd:split}...){p_end}
{p2col :{help duplicates}}duplicates, {cmd:isid}, {cmd:levelsof}{p_end}
{p2col :{help set_seed:set seed}}random numbers{p_end}
{p2col :{help power}}power and sample size{p_end}
{p2colreset}{...}

{title:More data management}

{p2colset 5 24 26 2}{...}
{p2col :{help datamgmt:drop}}{cmd:drop}, {cmd:keep}, {cmd:sort}, {cmd:gsort}, {cmd:rename}, {cmd:order}, {cmd:count}{p_end}
{p2col :{help label}}labels, {cmd:notes} and {cmd:char}{p_end}
{p2col :{help format}}display formats, {cmd:recast}, {cmd:compress}{p_end}
{p2col :{help codebook}}{cmd:codebook} and {cmd:misstable}{p_end}
{p2col :{help putexcel}}write results to Excel{p_end}
{p2colreset}{...}

{title:Estimation}

{p2colset 5 24 26 2}{...}
{p2col :{help regress}}linear regression ({cmd:areg}, {cmd:ivregress}){p_end}
{p2col :{help anova}}{cmd:anova} and {help oneway}{p_end}
{p2col :{help logit}}binary outcomes ({cmd:logistic}, {cmd:probit}){p_end}
{p2col :{help ologit}}ordered and multinomial outcomes ({cmd:oprobit}, {cmd:mlogit}){p_end}
{p2col :{help poisson}}count outcomes ({cmd:nbreg}){p_end}
{p2col :{help glm}}generalized linear models{p_end}
{p2col :{help tobit}}censored outcomes{p_end}
{p2col :{help qreg}}quantile regression{p_end}
{p2col :{help sem}}structural equation models ({cmd:estat gof}){p_end}
{p2col :{help xtreg}}panel data ({cmd:xtset}, {help xtsum}, {help xtlogit}, {help xtpoisson}, {help xtgee}){p_end}
{p2col :{help mean}}{cmd:mean}, {cmd:proportion}, {cmd:ratio}, {cmd:total}{p_end}
{p2col :{help nptests:ranksum}}nonparametric and normality tests{p_end}
{p2col :{help immediate:cii}}immediate commands ({cmd:cii}, {cmd:tabi}, {cmd:ttesti}){p_end}
{p2col :{help bootstrap}}bootstrap standard errors{p_end}
{p2colreset}{...}

{title:After estimation}

{p2colset 5 24 26 2}{...}
{p2col :{help predict}}predictions and residuals{p_end}
{p2col :{help margins}}margins and marginal effects; {help marginsplot}{p_end}
{p2col :{help test}}{cmd:test}, {cmd:testparm}, {cmd:lincom}, {cmd:nlcom}{p_end}
{p2col :{help estat}}{cmd:estat ic}, {cmd:vce}, {cmd:hettest}, {cmd:gof}, {cmd:classification}...{p_end}
{p2col :{help estimates}}store and compare results; {help lrtest}, {help hausman}{p_end}
{p2colreset}{...}

{title:Survey data}

{p2colset 5 24 26 2}{...}
{p2col :{help svy}}{cmd:svyset}, {cmd:svydescribe} and the {cmd:svy:} prefix{p_end}
{p2colreset}{...}

{title:Time series, survival, multilevel models and imputation}

{p2colset 5 24 26 2}{...}
{p2col :{help arima}}ARIMA models; {help newey}, {help dfuller}, {help corrgram}{p_end}
{p2col :{help var}}vector autoregressions ({cmd:varsoc}, {cmd:vargranger}); {help irf}{p_end}
{p2col :{help stset}}survival data ({help sts}, {help stsum}, {help stdescribe}){p_end}
{p2col :{help stcox}}Cox regression; {help streg} parametric models{p_end}
{p2col :{help mixed}}multilevel linear models; {help melogit}{p_end}
{p2col :{help mi}}multiple imputation{p_end}
{p2colreset}{...}

{title:Graphics}

{p2colset 5 24 26 2}{...}
{p2col :{help graph}}managing graphs (save, export, combine, schemes){p_end}
{p2col :{help twoway}}scatter, line, fits and other twoway plots{p_end}
{p2col :{help graph_bar:graph bar}}bar, box, dot, pie and matrix graphs{p_end}
{p2col :{help histogram}}histograms and {cmd:kdensity}{p_end}
{p2col :{help rvfplot}}residual plots; {help qnorm}, {help pnorm}, {help lowess}{p_end}
{p2colreset}{...}

{title:Interface and files}

{p2colset 5 24 26 2}{...}
{p2col :{help interface:browse}}Data Editor, Viewer, dialogs and Graph Editor{p_end}
{p2col :{help doedit}}the Do-file Editor{p_end}
{p2col :{help do}}run do-files; {cmd:capture}, {cmd:confirm}, {cmd:timer}{p_end}
{p2col :{help files:cd}}working directory and files{p_end}
{p2col :{help display}}display text and results{p_end}
{p2col :{help set}}settings{p_end}
{p2colreset}{...}

{title:Programming}

{p2colset 5 24 26 2}{...}
{p2col :{help program}}define programs; {cmd:return}, {cmd:ereturn}{p_end}
{p2col :{help macro:local}}macros, scalars and {cmd:tokenize}{p_end}
{p2col :{help filecmd:file}}read and write files{p_end}
{p2col :{help syntax}}parse a program's arguments{p_end}
{p2col :{help matrix}}matrices and matrix expressions{p_end}
{p2col :{help log}}record the session in a file{p_end}
{p2col :{help mata}}the matrix programming language{p_end}
{p2colreset}{...}

{title:OpenDTA additions}

{p2colset 5 24 26 2}{...}
{p2col :{help hints}}error explanations ({cmd:set hints}){p_end}
{p2col :{help superscript}}superscript exponents, 1.0×10¹⁰ ({cmd:set superscript}){p_end}
{p2col :{help numerics}}Stata arithmetic or more precise arithmetic ({cmd:set numerics}){p_end}
{p2col :{help language}}interface in English or Portuguese ({cmd:set language}){p_end}
{p2colreset}{...}

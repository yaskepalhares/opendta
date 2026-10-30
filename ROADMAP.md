# Roadmap

**English** | [Português](ROADMAP.pt-BR.md)

Goal: cover the commands of Stata 14, including graphics, with an interface as close as possible to the original.

**Version 1.0** is released when every command of Stata 14 is covered. The commands still missing are listed in [Missing for 1.0](#missing-for-10); what comes afterwards is in [After 1.0](#after-10).

The phases are in dependency order. Each one ends with something that can be opened, tested and compared with Stata (see [docs/testing.md](docs/testing.md)). Deadlines are not estimated in advance. Each phase is split into reviewable pull requests.

[![Commands implemented per area](https://yaskepalhares.github.io/opendta/roadmap-coverage.svg)](https://yaskepalhares.github.io/opendta/)

One tile per area of [Missing for 1.0](#missing-for-10): its size is the number of commands and its colour the share already implemented (options of existing commands are not counted). It is redrawn on every merge into `main`. **[Open the interactive version](https://yaskepalhares.github.io/opendta/)** to see the commands of each area and search for one.

## Phase 0: foundation ✅

- [x] Do-file lexer: comments, `///`, `#delimit`
- [x] Local and global macros, nested expansion, `=exp`, increment, basic extended functions
- [x] Expressions: precedence, missing-value rules, strings, a first set of functions
- [x] `display` with formats and directives
- [x] Control flow: `forvalues`, `foreach`, `while`, `if`/`else`, `continue`
- [x] `quietly`, `noisily` and `capture` prefixes; `do`, `run` and `include`
- [x] `r(#)` error codes and echo of commands and blocks as in Stata
- [x] Interface: Results, Command, Review, Variables, Properties, menus and toolbar
- [x] Console, batch mode and a tool to compare with Stata logs
- [x] Automated tests and continuous integration

## Phase 1: data in memory and files ✅

**1a ✅**
- [x] Dataset with the `byte`, `int`, `long`, `float`, `double`, `str#` and `strL` types, reproducing `float` precision
- [x] varlists: abbreviations, wildcards (`*`, `?`, `~`), ranges (`a-c`), `_all`
- [x] `if` and `in` qualifiers and the `by`/`bysort` prefix (weights are recognized by the syntax; their use arrives with the statistical commands)
- [x] `_n`, `_N`, `x[_n-1]`, `sum()` and vectorized expressions; sequential `replace` with subscripts
- [x] `clear`, `set obs`, `describe`, `list`, `generate`, `replace`, `drop`, `keep`, `rename`, `order`, `sort`, `gsort`, `label` (variable/define/values/list/dir/drop/data), `format`, `count`, `compress`, `input`
- [x] Variables and Properties windows linked to the data

**1b ✅**
- [x] OpenDTA's own `.dta` reader and writer for formats 117 (13), 118 (14) and 119, with value labels, notes and characteristics; cross-checked against ReadStat in the tests
- [x] Reading old `.dta` files (up to 115, Stata 12 and earlier) through ReadStat (`pip install pyreadstat`)
- [x] `use`, `save`, `saveold`, `recast`, `notes`, `char`, `type`, `erase`, `browse`/`edit`
- [x] `insheet`, `outsheet`, `import delimited`, `export delimited`
- [x] File > Open/Save, drag and drop of `.dta` files, read-only Data Browser

**1c ✅**
- [x] Formats 113–115 (Stata 8–12): OpenDTA's own reader; writing 114/115 with `saveold, version(11|12)`
- [x] `import excel`/`export excel` (.xlsx; .xls with the xlrd package)
- [x] `infile` (free format and with a dictionary), `infix`, `file open/write/read/close`
- [x] Data Editor with cell editing (generates `replace`/`set obs`/`generate`) and Variables Manager

## Phase 2: complete programming language ✅

- [x] `program define`/`end`, `syntax`, `args`, `gettoken`, `marksample`, `markout`
- [x] `return`, `ereturn`, `sreturn`; `tempvar`, `tempname`, `tempfile`; `preserve`/`restore`
- [x] Extended macro functions (`: variable label`, `: type`, `: format`, `: label`, `: char`, `: dir`, `: list`, `: subinstr`, `: permname`, `: r(scalars)`…)
- [x] `matrix` (definition, operators, functions, submatrices, `matrix list`, row and column names)
- [x] Search for `.ado` files on the adopath (runs SSC packages); `sysdir`, `adopath`, `which`, `findfile`, `discard`
- [x] `log using`/`log close` (text and SMCL), `cmdlog`; Viewer with an SMCL renderer, `help`, `view`, `search`
- [x] `assert`, `set trace`, `creturn list`, `timer`
- [x] Help pages for all commands, written from scratch

## Phase 3: descriptive statistics and data management ✅

**3a ✅**
- [x] `summarize` (`detail`, `meanonly`, weights), `tabstat`, `ci`, `ttest`, `prtest`, `correlate`, `pwcorr`, `centile`, `pctile`, `_pctile`, `xtile`
- [x] `tabulate` (one- and two-way, `row`/`col`/`cell`/`expected`, `chi2`, `lrchi2`, `V`, `gamma`, `taub`, `exact` by enumeration), `tab1`, `tab2`, `table`

**3b ✅**
- [x] `egen`: group functions (`count`, `mean`, `sd`, `total`, `min`, `max`, `median`, `mode`, `pctile`, `iqr`, `skew`, `kurt`, `mad`, `mdev`, `std`, `rank`), `group`, `tag`, `seq`, `fill`, `cut`, `concat`, `ends`, `diff` and row functions (`row*`, `any*`)
- [x] `collapse` (all statistics, weights, `cw`), `contract`, `expand`, `fillin`
- [x] `merge` (1:1, m:1, 1:m, m:m, `_n`, `update`/`replace`, `keep()`, `assert()`, `keepusing()`), `append`, `joinby`, `cross`
- [x] `reshape long`/`wide` (with `@`, `string` and repetition of the last specification)

**3c ✅**
- [x] `duplicates` (report/examples/list/tag/drop), `isid`, `levelsof`
- [x] `sample` (percentage or `count`, `if`/`in`, `by()` and the by prefix)
- [x] `recode`, `encode`/`decode`, `destring`/`tostring`, `split`, `mvencode`/`mvdecode`
- `distinct` is not an official command: it is an SSC package and runs through the adopath (phase 2).

**3d ✅**
- [x] 64-bit Mersenne Twister generator (mt64), `set seed` (number or state), `c(seed)`, `set rng`
- [x] `runiform()`, `runiformint()`, `rnormal()`, `rbinomial()`, `rpoisson()`, `rchi2()`, `rt()`, `rbeta()`, `rgamma()`, `rexponential()`, `rlogistic()`, `rweibull()`, `rnbinomial()`, `rhypergeometric()`. The seed conversion and the algorithm of each distribution are not documented: results are statistically equivalent but not identical to Stata's (cases in `compat/do/0307_aleatorios.do`).

**Checked against Stata**: the 28 cases in `compat/do/` have Stata references; 25 match line by line. The known differences are listed in `compat/known_differences.txt`: the progress of the network algorithm of `tabulate, exact` (tables larger than 2×2) and, for random numbers, `rnormal()`, `runiformint()`, `rpoisson()` and `rt()` (`runiform()`, `rbinomial()`, `rexponential()` and `rchi2(2)` already reproduce Stata's exact sequence).

## Phase 4: Mata ✅

**4a–4c ✅**
- [x] Lexer and parser with Mata's precedence; real, string, complex, pointer and struct types; scalars, vectors and matrices
- [x] Operators (arithmetic, colon, transpose, Kronecker, `,` and `\` joins, `..` and `::` ranges, relational, logical, `?:`), `[i,j]` and `[|...|]` subscripts, missing values propagated as in Stata
- [x] Control flow (`if`, `for`, `while`, `do`, `break`, `continue`, `return`), typed user functions, optional arguments and pass by reference, `struct`, pointers
- [x] `mata`/`mata:` ... `end` blocks with Stata's echo, `mata: statement`, `mata clear/describe/drop`; errors in Mata's format (`<istmt>:  3499  x not found`)
- [x] Library: ~150 functions (mathematical, matrix, statistical, missing values, strings, `printf`/`sprintf`)
- [x] Link to the data: `st_data`, `st_sdata`, `st_view` (copy), `st_store`, `st_addvar`, `st_local`, `st_global`, `st_numscalar`, `st_matrix`, `stata()`

**4d ✅**
- [x] Display, messages and errors checked against Stata (cases `compat/do/0401`–`0405`): numbers in `%12.0g`, single-width columns, block header, compile and run-time errors, skipped lines in `mata:`
- [x] `st_view` that writes to the data when elements are assigned; associative arrays (`asarray`); `optimize()` with d0/d1/d2 evaluators and nr/bfgs techniques
- [x] Check `optimize()`, `st_view` and `asarray` (cases `0406` and `0407`; in d0 `optimize()` only the noise of numerical derivatives differs)
- `moptimize()` moved to phase 5, together with `ml`; Mata classes, to phase 8

## Phase 5: estimation ✅

- [x] Factor variables (`i.`, `c.`, `#`, `##`) and time-series operators (`L.`, `F.`, `D.`, `S.`); `tsset`, `xtset`
- [x] `regress`, `predict`, `test`, `testparm`, `lincom`, `nlcom`, `margins`, `estimates store/table/restore`
- [x] `logit`, `logistic`, `probit`, `ologit`, `oprobit`, `mlogit`, `poisson`, `nbreg`, `glm`, `tobit`, `ivregress`, `xtreg`, `areg`
- [x] Robust, cluster and bootstrap standard errors (`vce(bootstrap)` and the `bootstrap` prefix); `ml` engine (maximum likelihood)
- [x] `xtlogit` and `xtpoisson` (fe and re), `qreg`, `lrtest`, `hausman`, `xtsum`, `xtdescribe`
- [x] Mata's `moptimize()` (lf, d0–d2, lf0–lf2 and gf0–gf2 evaluators; oim, opg and robust variance; case `0408` checked against the Stata 14 log)
- Checked against the Stata logs: cases `0501`–`0513`

## Phase 6: complex surveys ✅

- [x] `svyset` (strata, PSUs, weights, FPC), `svydescribe`
- [x] `svy:` prefix with `mean`, `proportion`, `ratio`, `total`, `tabulate`, `regress`, `logit`, `logistic`, `poisson`; `subpop()` and `over()`
- [x] Taylor linearization; Rao–Scott correction in two-way `svy: tabulate`
- [x] Several stages, BRR, jackknife, bootstrap and SDR; `estat effects`; `svy:` also with `ologit`, `oprobit`, `mlogit`, `nbreg`, `glm` and `tobit`
- Checked against the Stata logs: cases `0601`, `0602` and `0603`

## Phase 7: graphics ✅

- [x] Graphics engine with OpenDTA's own schemes (`s2color`, `s1mono`… are accepted and mapped to them) and Graph window
- [x] `twoway` (`scatter`, `line`, `connected`, `lfit`, `qfit`, `lfitci`, `qfitci`, `area`, `bar`, `spike`, `dropline`, `rcap`, `rarea`, `rbar`, `rline`, `function`, `histogram`, `kdensity`, `scatteri`), `tsline`, `graph bar`/`hbar`, `graph box`/`hbox`, `graph pie`, `histogram`, `kdensity`, `graph matrix`
- [x] Title, axis, legend and marker options; `by()`
- [x] `graph combine`, `graph export` (PNG, PDF, SVG, EPS, PS, TIFF, JPEG), `graph save`/`use` (OpenDTA's own format), `graph dir`/`drop`/`rename`/`copy`/`display`/`describe`
- [x] `lowess`, `graph dot`, `marginsplot`, diagnostic plots (`rvfplot`, `rvpplot`, `qnorm`, `pnorm`)
- [x] Simple Graph Editor (titles, legend, scheme, size, grid, log scales), with undo
- Checked against the Stata logs: case `0701` (text messages only; the drawing is OpenDTA's own)

## Phase 8: the remaining core commands and interface tools ✅

- [x] Do-file Editor with tabs, syntax highlighting, find/replace and execution of the selection (`doedit`); Data, Graphics and Statistics menus; printing
- [x] Nonparametric and normality tests (`ranksum`, `signrank`, `kwallis`, `spearman`, `swilk`, `sktest`, `sdtest`) and immediate commands (`cii`, `ttesti`, `prtesti`, `tabi`)
- [x] `oneway`, `anova`, `tabulate, summarize()`, `estat` (`ic`, `vce`, `summarize`, `hettest`, `ovtest`, `classification`, `gof`, `bootstrap`)
- [x] `codebook`, `misstable`, `putexcel`
- [x] Mata classes (`class`, `extends`, virtual methods, `super`, `this`, static members, constructors; case `0409` checked against the Stata 14 log)
- [x] Survival analysis: `stset` (with `id()`, `time0()`, `origin()`, `enter()`, `exit()`, `scale()`), `st`, `streset`, `stdescribe`, `stsum`, `sts list`/`generate`/`test`/`graph` (Kaplan-Meier, Nelson-Aalen, log-rank, Wilcoxon, Tarone-Ware, Peto-Peto), `stcox` (Breslow and Efron ties, strata, robust), `estat phtest`, `streg` (exponential, Weibull, Gompertz, lognormal, loglogistic) and their `predict` options; cases `0901`–`0902` checked against the Stata 14 log
- [x] Multilevel models: `mixed` (ML and REML, nested levels, random slopes, identity/independent/exchangeable/unstructured covariances; `estat icc`, `group`, `recovariance`; `predict`) and `melogit`, `meprobit`, `mepoisson` (adaptive Gauss-Hermite quadrature); cases `1001`–`1002` checked against the Stata 14 log
- [x] Time series: `arima` (exact ML with the Kalman filter, OPG variance), `var`, `varsoc`, `vargranger`, `dfuller` (Fuller critical values, MacKinnon p-values), `corrgram`, `newey`; case `1101` checked against the Stata 14 log
- [x] Multiple imputation: `mi set` (mlong, flong), `register`, `describe`, `impute` (`regress`, `pmm`, `logit`, `ologit`, `mlogit`, `poisson`, `chained`), `estimate` (Rubin's rules, Barnard-Rubin df), `passive`, `xeq`, `extract`, `misstable`; case `1201` (imputed values differ from Stata's: known difference)
- [x] `sem` (path models, CFA and latent variables by ML, standardized solution, `estat gof`, `predict`), `power` (`onemean`, `twomeans`, `pairedmeans`, `oneproportion`, `twoproportions`, `onecorrelation`, `twocorrelations`) and `irf` (`create`, `table`, `graph`, `describe`, `rename`, `drop` after `var`); cases `1301`–`1303` checked against the Stata 14 log
- [x] Completed commands: `anova` (sequential sums of squares, nested terms, error terms, repeated measures), `xtgee` and the `pa` estimator of `xtreg`/`xtlogit`/`xtpoisson`, `xtreg, be` and `mle`, two-way `tabulate, summarize()`, `misstable tree`/`nested`, `ci`/`cii variances`, two-sample `prtest`/`prtesti`, `reshape i/j/xij/error`, the old `merge` syntax, `sort ... in`, `list, sepby()`, `display _request()`, `export excel` to `.xls`, `set rng kiss32`; sequential `replace` with time-series operators (`replace y = 2*L.y`); case `0805` checked against the Stata 14 log
- [x] Menu dialogs (38 Statistics and Graphics commands), with a preview of the command and OK, Submit and Copy
- [x] Help pages for all commands
- [x] macOS package (`tools/build_macos.sh`: `OpenDTA.app` and `.dmg`; "Open With" for `.do` and `.dta`)
- Checked against the Stata logs: cases `0801`–`0804`

## Phase 9: interface language ✅

- [x] Settings window with a sidebar and search: theme (light, dark or system), app icon, font with preview, language, error explanations, superscript exponents, numerics and updates, applied at once
- [x] Language selector in Settings: English (default) or Portuguese, for menus, windows, dialogs, error explanations and OpenDTA's own messages
- [x] With Portuguese on, Stata-style messages can stay in English (for those who compare with Stata logs and manuals) or be translated, at the user's choice (`set language`, `set messages`)

## Phase 10: epidemiology ✅

- [x] Tables for epidemiologists: `cs`/`csi` (risk difference and ratio, attributable or prevented fractions), `cc`/`cci` (odds ratio with exact, Cornfield, Woolf or test-based limits; Fisher's exact test), `ir`/`iri` (rate difference and ratio, exact limits, mid-p), `mcc`/`mcci` (McNemar's test and paired comparisons), with `by()` stratification: crude and Mantel-Haenszel estimates and tests of homogeneity
- [x] `tabodds` (odds by exposure level, homogeneity and trend tests, odds ratios adjusted by Mantel-Haenszel, graph) and `mhodds` (Mantel-Haenszel odds ratio, score test for trend)
- [x] `symmetry` (Bowker, Stuart-Maxwell, Bickenboller, exact significance), `dstdize` (direct standardization) and `istdize` (indirect standardization, SMR with exact limits)
- [x] Interrater agreement: `kap` (Cohen's kappa, weighted kappa; Fleiss's kappa for several raters) and `kappa` (frequency form)
- [x] ROC analysis: `roctab` (area with DeLong or Hanley-McNeil errors, detailed report, curve), `roccomp` (equality of areas, correlated or independent), `rocfit` (binormal model), `lroc` and `lsens` after `logit`, `logistic` and `probit`
- Checked against the Stata 14 log: cases `1401` and `1402` (the differences that remain, and why, are in `compat/known_differences.txt`)

## Phase 11: data management ✅

- [x] Variables and observations: `clonevar`, `insobs`, `range`, `separate`, `ipolate`, `expandcl`, `stack`, `xpose`, `drawnorm`, `corr2data`
- [x] Describing and comparing data: `ds`, `lookfor`, `inspect`, `labelbook`, `compare`, `cf`, `datasignature`, `snapshot`
- [x] Files and the system: `copy`, `mkdir`, `rmdir`, `dir`/`ls`, `shell` (and `!`), `filefilter`, `zipfile`, `unzipfile`, `checksum`, `outfile`, `unicode`
- [x] `statsby`, `putmata`; `sysuse` (OpenDTA's own example datasets) and `webuse` (downloads on demand)
- [x] `icd9` and `icd10` with the public-domain ICD-9-CM and ICD-10-CM code lists of the U.S. CMS (check, clean, generate, lookup, search); `icd9p` by format
- [x] `odbc` through the optional `pyodbc` package; business calendars (`bcal`, the `%tb` format, `bofd()` and `dofb()`)
- Checked against the Stata 14 log: case `1501` (matches)

## Phase 12: epidemiology, finished ✅

- [x] `cs`, `cc` and `ir` with `by()`: directly pooled estimates (`pool`) with their test of homogeneity, internally and externally standardized estimates (`istandard`, `estandard`) or with weights from a variable (`standard()`), standardized risk and rate differences (`rd`, `ird`), binomial frequency records (`binomial()`), and the Breslow-Day and Tarone tests of `cc` (`bd`, `tarone`)
- [x] `mhodds` with `by()`: the odds ratio within each group, the estimate controlling for the groups and a test of homogeneity
- [x] `roctab` with Bamber's standard error (`bamber`), exact binomial limits (`binomial`) and the Lorenz curve with the Gini and Pietra indices (`lorenz`); `roccomp` with binormal areas (`binormal`)
- [x] `rocreg`: nonparametric ROC estimation with covariate-adjusted controls (stratified or linear, empirical or normal percentile values) and its summaries (AUC, ROC(t), inverse ROC, partial AUC), and the parametric probit GLM model with ROC covariates; bootstrap standard errors with normal, percentile and bias-corrected intervals
- Checked against the Stata 14 log: case `1601` (the differences that remain, and why, are in `compat/known_differences.txt`)

## Phase 13: data management, options ✅

- [x] `file` with binary read and write (`%1b`, `%2b`, `%4b`, their signed and unsigned forms, `%4z`, `%8z`, `%#s`), `file seek`, `file set byteorder`, `file query`, `read write` and the `byteorder()` function
- [x] `infile` with `byvariable()`; `putexcel` to `.xls`
- [x] `import sas` (with `bcat()`), `import spss` and `import sasxport`, and a basic `export sasxport`, through the optional `pyreadstat` package
- [x] `describe` with `numbers`, `varlist`, `replace` and `using`; `list` with `header()`, `nodotz`, `mean`, `sum`, `N`, `labvar()` and `nocompress`
- Checked against the Stata 14 log: case `1701` (matches)

## Phase 14: programming ✅

- [x] `unab`, `tsunab`, `fvunab`, `unabcmd`, `fvexpand`, `fvrevar`, `tsrevar` and `_rmcoll`, with factor-variable and time-series operators (specific forms such as `i(2,3)bn.rep78` included)
- [x] `mark`, `markin` and `markout` (with `sysmissok`); `marksample` and `syntax varlist(fv ts)` with factor variables and time-series operators
- [x] `postfile`, `post`, `postclose` and `postutil`
- [x] `matlist` (column and row specifications, borders, formats, titles) and `tabdisp` (one to three variables, `by()`, totals, `concise`, `dotz`)
- [x] `program, sortpreserve`; `program name, plugin` (explains that Stata plugins cannot be loaded)
- [x] `pause`, `sleep`, `window` (in batch mode the dialogs act as cancelled) and `viewsource` (shows OpenDTA's own source)
- [x] `signestimationsample`, `checkestimationsample`, `makecns` with a constraint matrix and `matcproc`
- Checked against the Stata 14 log: case `1801` (the differences that remain, and why, are in `compat/known_differences.txt`)

## Phase 15: survival analysis ✅

- [x] `strate`, `stptime`, `stci` (percentiles, restricted and extended means), `stir`, `stmh`, `stmc` and `ltable`
- [x] `stsplit`, `stjoin`, `stfill`, `stgen`, `stvary`, `stbase` and `snapspan`
- [x] Count-time data: `ctset`, `ct`, `cttost`, `sttoct`; case-control sampling with `sttocc`
- [x] `stcrreg` (Fine-Gray competing risks) and `stcurve` after `stcox`, `streg` and `stcrreg`
- [x] `streg, distribution(ggamma)` and `frailty(gamma|invgaussian)`, with `shared()`; the `exactm` and `exactp` ties of `stcox`
- Checked against the Stata log: case `1901`

## Phase 16: multilevel models ✅

- [x] `meglm` (families gaussian, bernoulli, binomial, gamma, nbinomial, ordinal and poisson with their links), `mecloglog`, `menbreg` (mean and constant dispersion), `meologit`, `meoprobit`, `meqrlogit` and `meqrpoisson`; `predict` and `estat icc`/`recovariance`/`group` after them
- [x] `mestreg` (exponential, Weibull, lognormal, loglogistic and gamma with a random effect per group)
- [x] Crossed random effects in `mixed` (`|| _all: R.factor`)
- Checked against the Stata 14 log: case `1902` (the differences that remain, and why, are in `compat/known_differences.txt`)

## Phase 17: base statistics ✅

- [x] Resampling prefixes `jackknife` (observations or clusters, `mse`, `dof()`), `permute` and `simulate`
- [x] `ml` with the `lf`, `d0`, `d1` and `d2` evaluators (`ml model`, `init`, `search`, `maximize`, `display`, `clear`, `score`; `oim`, `opg`, `robust` and `cluster()` variances; numerical derivatives), `mleval`, `mlsum`, `mlvecsum`, `mlmatsum` and `mlexp`
- [x] Model building: `mkspline` (linear, `pctile`, `cubic`, `marginal`), `xi`, `nestreg`, `stepwise` and `fp`/`mfp`
- [x] After estimation: `testnl`, `predictnl`, `linktest`, `contrast`, `pwcompare`, `pwmean` and `suest`; `margins` after `tobit`, `ologit` and `oprobit`
- [x] Transformations: `boxcox` (left-hand-side model), `ladder`, `gladder` and `lnskew0`
- [x] Tests and correlations: `bitest`/`bitesti`, `ksmirnov`, `median`, `nptrend`, `runtest`, `tetrachoric`, `ktau` and `pcorr`
- [x] Agreement and effect sizes: `icc`, `loneway`, `esize`/`esizei` and `ameans`
- [x] Displays and graphs: `cumul`, `stem`, `lv`, `smooth`, `quantile`, `symplot`, `serrbar`, `dotplot`, `sunflower`, `lpoly`, `cusum`
- Checked against the Stata 14 log: case `1903` matches except for the differences in `compat/known_differences.txt` (starting steps of `ml`/`ologit`, the `lnskew0` search, the depth-1 pseudosigma of `lv`, `permute`/`simulate` draws and the `sunflower` bin table)

## Phase 18: time series ✅

- [x] Data tools: `tsappend`, `tsfill`, `tsreport` and the `rolling` prefix
- [x] Filters and smoothers: `tsfilter` (`hp`, `bk`, `cf`, `bw`) and `tssmooth` (`ma`, `exponential`, `dexponential`, `hwinters`, `shwinters`, `nl`)
- [x] Spectral analysis and white-noise tests: `pergram`, `cumsp`, `xcorr`, `wntestb`, `wntestq` and `psdensity` after `arima` and `arfima`
- [x] Unit roots and regression with AR(1) errors: `dfgls`, `pperron` and `prais` (Prais-Winsten, Cochrane-Orcutt, `rhotype()`, `ssesearch`, `twostep`)
- [x] VAR postestimation: `varlmar`, `varnorm`, `varstable`, `varwle`, `varbasic` and `fcast` (after `var`, `svar` and `vec`)
- [x] Cointegration and structural VAR: `vecrank`, `vec` and `svar` (short-run `aeq()`/`beq()` matrices); `irf` after `svar` and `vec`, with bootstrap standard errors (`bs`, `bsp`), and `irf ctable`, `cgraph`, `ograph` and `add`
- [x] Conditional variance and long memory: `arch` (ARCH, GARCH and the other variance terms, ARCH-M, ARMA errors, t and GED distributions, `het()`), `arfima` (exact ML) and `mgarch` (`dvech`, `ccc` and `dcc`)
- [x] `arima` with `sarima()`, `mar()`, `mma()`, `condition` and `diffuse`
- [x] State-space models by the Kalman filter with diffuse initial states: `ucm`, `sspace` and `dfactor`
- Checked against the Stata 14 log: case `1904` matches except for the differences in `compat/known_differences.txt` (the `dfgls` critical values by lag, the iteration logs of `prais, ssesearch` and `svar`, the `fcast` standard errors from the fourth step on, the bootstrap tables of `irf`, and the estimates of `ucm` with `cycle()` and of `mgarch dcc`, where Stata reaches a different point)

## Phase 19: panel data ✅

- [x] Random-effects models by Gauss-Hermite quadrature (adaptive `mvaghermite` and `ghermite`, `intpoints()`): `xtprobit`, `xtcloglog`, `xtologit`, `xtoprobit`, `xttobit`, `xtintreg` and `xtnbreg` (random effects, conditional fixed effects and `pa`), with the LR test of the random effect and `vce(robust)`/`vce(cluster)`; `xtlogit, re` moved to the same engine; `quadchk`
- [x] Linear panel estimators: `xtgls` (`panels()`, `corr()`, `igls`), `xtpcse`, `xtregar` (`fe`, `re`, `rhotype()`, `rhof()`, `twostep`), `xtrc`, `xthtaylor` (`endog()`, `constant()`, `varying()`, `amacurdy`) and `xtivreg` (G2SLS and EC2SLS random effects, `be`, `fe`, `fd`)
- [x] `xtfrontier`: time-invariant (`ti`) and time-varying decay (`tvd`) inefficiency, production and `cost` frontiers
- [x] `xtunitroot`: `ht`, `breitung`, `hadri` and `fisher`
- [x] Dynamic panel data by GMM: `xtabond`, `xtdpdsys` and `xtdpd` (one-step and two-step, Windmeijer-corrected robust errors, `pre()`, `endogenous()`, `inst()`, `maxldep()`), with `estat abond` and `estat sargan`
- [x] Description: `xtline` (with `overlay`), `xttab` and `xtdata` (`fe`, `be`)
- Checked against the Stata 14 log: case `1905` matches except for the differences in `compat/known_differences.txt` (the iteration paths of the numerical maximizers and of `igls`, `xtgls` with correlated panels and `igls` with `psar1`, the variance components of `xtregar, re`, the F test of `u_i` in `xtivreg, fe`, the between R-sq of `xtivreg, fd`, the flat likelihood of `xtfrontier, tvd` and the Breitung statistic); `xtunitroot llc` and `ips`, which Stata 14 runs, are not implemented yet

## Phase 20: multivariate statistics ✅

- [x] Principal components and factor analysis: `pca`, `factor` (`pf`, `pcf`, `ipf`, `ml`), `rotate` (`varimax`, `quartimax`, `equamax`, `parsimax`, `promax`, `oblimin`, `quartimin`), `alpha`, `screeplot`, `loadingplot` and `scoreplot`, with `estat` (`kmo`, `smc`, `anti`, `loadings`, `correlation`, `residuals`, `summarize`, `structure`, `common`) and `predict`
- [x] Multivariate tests and models: `hotelling`, `mvreg`, `manova` (Wilks, Pillai, Lawley-Hotelling and Roy), `canon` and `mvtest` (`means`, `covariances`, `correlations`, `normality`)
- [x] Cluster analysis: `cluster` (`kmeans`, `kmedians`, the seven hierarchical linkages, `generate`, `dendrogram`, `stop` with `calinski` and `duda`, `list`, `dir`, `drop`, `rename`, `notes`) and `matrix dissimilarity`
- [x] Discriminant analysis: `discrim` (`lda`, `qda`, `knn`, `logistic`) and `candisc`, with `estat` (`classtable`, `errorrate`, `list`, `grdistances`, `covariance`) and `predict`
- [x] Correspondence analysis and scaling: `ca`, `mca`, `mds` and `mdsmat` (classical and modern), orthogonal `procrustes` and `biplot`, with `estat` and `predict` where they apply
- Checked against the Stata 14 log: case `1906` matches except for the differences in `compat/known_differences.txt` (the iteration paths of `factor, ml`, `discrim logistic` and `mds, method(modern)`, and the maximum-likelihood factor numbers that follow from where Stata's `ml` stops)

## Phase 21: models, part 1 (binary and count, censored and selection, instrumental variables) ✅

- [x] Binary and fractional response: `cloglog`, `scobit`, `hetprobit` (`het()`), `fracreg` (`logit`, `probit`, `cloglog`), `binreg` (`or`, `rr`, `hr`, `rd`) and `glogit`
- [x] Count data: `zip` and `zinb` (`inflate()`, `probit`, `vuong`, `irr`), `tpoisson`, `tnbreg` and `gnbreg` (`lnalpha()`)
- [x] Censored, truncated and selection models: `intreg`, `truncreg` (`ll()`, `ul()`), `heckman` (maximum likelihood and `twostep`), `heckprobit`, `heckoprobit` and `biprobit` (shared varlist and `(y1 = x) (y2 = z)` forms), with the bivariate normal distribution function in `core/bvnorm.py` (Plackett's identity, exact derivatives)
- [x] Instrumental variables: `ivregress liml` (k-class) and `gmm` (two-step, `wmatrix()`), `ivprobit` and `ivtobit` (maximum likelihood, Wald test of exogeneity), `ivpoisson gmm` (additive and `multiplicative` moments, `onestep`) and `cfunction`
- [x] A multi-equation maximum-likelihood engine (`estimation/mleq.py`) with the log likelihoods apart from the commands (`estimation/likelihoods.py`), reused by the choice models of the next part
- Checked against the Stata 14 log: case `1907` matches, except for the differences in `compat/known_differences.txt` (the iteration paths and the stopping point of each optimizer; blocks whose data put a parameter at the boundary, so that Stata's estimates depend on its path: `scobit`, `zinb`, `tnbreg`, `gnbreg`; `binreg, hr`, where Stata's IRLS does not converge on the case data; `glogit`, run on other grouped counts than Stata's; the robust variance of `ivregress liml` and of `ivpoisson cfunction`). Found and fixed on the way: `binreg` is IRLS with the expected-information variance, `ivpoisson cfunction` uses the multiplicative residual, `ivregress gmm` with a robust weight matrix reports the sandwich variance, and the constant-only models of `tpoisson` and `gnbreg` follow Stata's reference points
- Not verified against Stata yet: `glogit` (data of the case), `fracreg cloglog`, `heckoprobit` and the options not exercised by case `1907`

## Phase 22: models, part 2 (choice, exact, systems and treatment effects) ✅

- [x] Systems and constraints: `constraint`, `cnsreg`, `sureg` (including `isure`), `reg3` (`2sls`, `3sls`, `ols`, `sure`), `eivreg`, `vwls` and `rreg` (Huber then biweight), with the algebra in `core/syseq.py`
- [x] Choice models: `clogit` (any number of chosen alternatives per group), `rologit` (`exactp`, `breslow`, `efron` and `none` for ties), `asclogit`, `slogit`, `mprobit`, `asmprobit` and `asroprobit` (GHK simulator with Hammersley draws), with the likelihoods in `core/choice.py`
- [x] Exact methods: `exlogistic` and `expoisson` (conditional maximum likelihood, median unbiased estimates, exact confidence intervals and p-values, from the exact distribution of the sufficient statistics in `core/exact.py`)
- [x] Nonlinear models: `nl`, `nlsur` (NLS, FGNLS and IFGNLS), `gmm` (one- and two-step, robust and cluster weight matrices) with substitutable expressions, and `frontier` (`hnormal`, `exponential`, `tnormal`)
- [x] Treatment effects: `teffects ra`, `ipw`, `ipwra`, `aipw`, `nnmatch` and `psmatch` (Abadie-Imbens variance), `eteffects` (control function), `stteffects ra`, `ipw` and `ipwra`, `etregress` (maximum likelihood and `twostep`) and `etpoisson` (Gauss-Hermite), with the estimators in `core/treatfx.py`
- Checked against the Stata 14 log: case `1908` (see `compat/known_differences.txt` for what does not match and why)
- Not verified against Stata yet: the table layouts of every command of this part (marked `VERIFY` in the code), the variance of `eivreg`, the weights of `stteffects` and the options the case does not exercise

## Phase 23: other suites (sem, gsem, irt, Bayesian estimation, power, mi impute and svy) ✅

- [x] `sem` with `method(mlmv)` (normal likelihood by missing-data pattern, EM starting values in `core/mvnmiss.py`) and `method(adf)`, `group()` and `ginvariant()`, path labels (`@name`) and the `estat` subcommands `eqgof`, `eqtest`, `framework`, `ggof`, `ginvariant`, `mindices`, `residuals`, `scoretests`, `stable`, `stdize:`, `summarize` and `teffects`
- [x] `gsem` with the `gaussian`, `bernoulli`, `binomial`, `poisson`, `nbinomial`, `gamma` and `ordinal` families, observation-level and cluster-level (`M[id]`) latent variables, adaptive Gauss-Hermite quadrature and `predict`
- [x] `irt` (`1pl`, `2pl`, `3pl`, `grm`, `nrm`, `pcm`, `gpcm`, `rsm` and `hybrid`) on top of `gsem`, with `estat report` and `irtgraph` (`icc`, `iif`, `tcc`, `tif`)
- [x] `bayesmh` (adaptive random-walk Metropolis-Hastings in `core/mcmc.py`, eight likelihoods and ten priors), the `bayes` prefix for `regress`, `logit`, `probit` and `poisson`, `bayesstats` (`summary`, `ess`, `ic`), `bayestest` (`interval`, `model`) and `bayesgraph`
- [x] `power` with `onevariance`, `twovariances`, `pairedproportions` and `oneway`, the `lrchi2` and `fisher` tests of `twoproportions` and the `graph` option
- [x] `mi impute mvn` (data augmentation), `monotone`, `truncreg`, `intreg` and `nbreg`
- [x] `svy:` with `cloglog`, `intreg`, `truncreg`, `heckman`, `heckprobit`, `zip`, `zinb`, `gnbreg`, `hetprobit`, `scobit`, `biprobit`, `fracreg`, `slogit`, `mprobit`, `clogit` and `ivregress`
- Checked against the Stata 14 log: case `1909` (the differences and the reason for each are in `compat/known_differences.txt`). Found and fixed on the way: the variance of `sem, method(adf)` is the numerical Hessian of the discrepancy, not its Gauss-Newton approximation; `ivregress` ignored weights; the `svy:` confidence intervals of `sigma`, `rho`, `alpha` and `lambda` use the design degrees of freedom; `irt 3pl` has one common pseudoguessing parameter unless `sepguessing` is given. The draws of `bayesmh` and `mi impute` come from OpenDTA's own generator and do not reproduce Stata's

## Phase 24: programming, part 2 ✅

- [x] `makecns` with a list of constraint numbers (constraints that cannot be read are left out with a note), the `displaycns` and `nocnsnotes` options, and `matrix dispCns`
- [x] `constraint #` without the word `define`
- [x] `matlist` with equation names: `showcoleq()` (`first`, `each`, `combined`, `lcombined`, `rcombined`), `roweqonly`, `coleqonly`, `keepcoleq`, `colorcoleq()` and `lines(eq)`; `matrix rownames` and `colnames` set the equation of names written `eq:name`
- [x] `window` in the graphical interface: `fopen` and `fsave` file dialogs, `stopbox` message boxes, `push` to the Review window, the User menu built with `window menu`, and `manage` (`minimize`, `restore`, `maintitle`, `forward`, `close`)
- Checked against the Stata 14 log: case `1910`. Found and fixed on the way: `makecns` skips numbers with no constraint, keeps redundant constraints and shows nothing with `r`; `matrix dispCns` needs results posted by an estimation command; `matlist` shows the column equations as spans, a label row for each row equation and, in later panels, leaves out the empty upper half of a symmetric matrix

## Phase 25: data management, part 2 ✅

- [x] `table` with a supercolumn variable (`table rowvar colvar supercolvar`), the `scol` option and `by()` with up to four variables
- [x] `list` with `string()` and the `display` format (also chosen when the table does not fit the line)
- [x] `import sasxport` with `describe` and `member()` (OpenDTA's own reader of every member of an XPORT file, `io/xport.py`) and the value labels of `formats.xpf`; `export sasxport` with `rename` and `vallabfile()` (`xpf`, `sascode`, `both`, `none`)
- Checked against the Stata 14 log: case `1911`. Found and fixed on the way: `table` keeps every row and every combination of the `by()` variables, with the supercolumns 4 blanks apart; `list, string()` ends cut values with `..`; `import sasxport, describe` has Stata's layout and stored results

## Phase 26: epidemiology, part 2 ✅

- [x] `rocreg` with `probit ml`: the binormal ROC model by maximum likelihood, with ROC covariates (`roccov()`), control covariates (`ctrlcov()`, linear) and the AUC by the delta method
- [x] `rocreg` with `cluster()`: cluster-robust standard errors after `ml` and resampling of whole clusters in the bootstrap
- Checked against the Stata 14 log: case `1912` (the differences and the reason for each are in `compat/known_differences.txt`). Found and fixed on the way: `probit ml` fits, for each classifier, the covariate control adjustment model (case and control means and standard deviations) and takes the ROC parameters and the AUC from it by the delta method; several classifiers go in one model, with the test of equal AUCs

## Phase 27: base statistics, part 2 ✅

- [x] `ml` with the `lf0`, `lf1`, `lf2`, `gf0`, `gf1` and `gf2` methods, `constraints()` (maximization over the constrained subspace, `e(Cns)`) and `svy` (design weights and linearized variance)
- [x] `boxcox` with `model(rhsonly)`, `model(lambda)` and `model(theta)` and `notrans()`
- [x] `icc` with `testvalue()` and unbalanced one-way data; `loneway` with `mean`, `median` and `exact`
- [x] `esize` and `esizei` with `unequal` and `welch`; `pwmean` with `unequal` (Satterthwaite degrees of freedom per comparison); `ameans` with `fweight` and `aweight`
- [x] `lv` with `generate` and `tail()`; `stem` with `width()` (and a fix for decimal data stored as `float`)
- Checked against the Stata 14 log: case `1913` (the differences and the reason for each are in `compat/known_differences.txt`). Found and fixed on the way: `ml` with constraints posts e(Cns) and the survey results of Stata; `boxcox` iterates by Newton–Raphson from 1 and transforms every variable of the varlist; `icc` leaves out the targets with missing ratings; `loneway` stores Stata's r(); `esize, unequal/welch` report their degrees of freedom

## Phase 28: panel data, part 2 ✅

- [x] `xtivreg` with `first` (first-stage regressions with the same estimator) and `regress` (every regressor exogenous)
- [x] `xtdata` with `re` and `ratio()` (the random-effects GLS transform, with a `constant` variable)
- [x] `xtrc` with `betas` (best linear predictors of the panel-specific coefficients, `e(beta_ps)` and `e(V_ps)`)
- [x] `xtunitroot fisher` with `drift`, `demean` and `pperron`; `xtunitroot hadri` with `kernel()` (Bartlett, Parzen and quadratic-spectral long-run variances with a given number of lags)
- [x] `xtabond`, `xtdpdsys` and `xtdpd` with `artests()`; `xtdpd` with `hascons`
- Checked against the Stata 14 log: case `1914` (the differences and the reason for each are in `compat/known_differences.txt`). Found and fixed on the way: the first-stage regressions of `xtivreg, first` follow the estimator (G2SLS has its own header); `xtdata, re` keeps the panel variable and adds the constant; `xtunitroot fisher` and `hadri` take Stata's layout and r(); `xtabond` drops duplicated instruments and the constant of `hascons` with a note

## Phase 29: time series, part 2 ✅

- [x] `prais` with `vce(robust)`, `vce(hc2)`, `vce(hc3)` and `vce(cluster)` (semirobust variance of the transformed regression)
- [x] `arch` with `constraints()` and `from()`; `arfima` with `vce(robust)`
- [x] `svar` with `aconstraints()` and `bconstraints()`; `sspace` and `dfactor` with `constraints()` on the coefficients
- [x] `predict` after `ucm`, `sspace` and `dfactor` (one-step predictions, residuals and standardized residuals); the `frequency()` suboption of `ucm`'s `cycle()`; `psdensity` after `ucm`
- Checked against the Stata 14 log: case `1915` (the differences and the reason for each are in `compat/known_differences.txt`). Found and fixed on the way: the log likelihood of `prais` leaves out the Jacobian term; `arfima` refines its estimates in a second maximization; `svar` posts every entry of A and B in e(b), with zero variance for the fixed ones, and e(Cns); `ucm` checks the cycle frequency and its Wald test covers the frequency and the damping; `sspace` and `dfactor` standardize the residuals jointly

## Phase 30: multilevel models, part 2 ✅

- [x] `meglm` and the other `me` commands with nested levels of random effects (`|| school: || class:`), by product adaptive quadrature over the levels
- [x] Weights in the `me` commands and `mestreg`: `fweight`, `pweight` and `iweight` at the observation level and the `fweight()` and `pweight()` options of each random-effects equation
- [x] `vce(robust)` and `vce(cluster clustvar)` (the default with sampling weights): sandwich variance with the scores of the highest-level groups
- [x] `predict` after `mestreg` with `median`, `mean` and `hazard` for every distribution
- Checked against the Stata 14 log: case `1916` (the differences and the reason for each are in `compat/known_differences.txt`). Found and fixed on the way: with nested levels the quadrature points of a lower level move with the effects of the level above and the derivatives follow the points, which gives Stata's estimates and standard errors; a variance going to zero no longer inflates the likelihood; `meqrlogit` and `meqrpoisson` take nested levels; frequency weights count in the group table; the `mestreg` median uses the fixed effects only

## Phase 31: multivariate statistics, part 2 ✅

- [x] `clustermat` (hierarchical clustering of a dissimilarity matrix, with `shape()`, `add`, `clear` and `labelvar()`); the binary similarity measures of `cluster` and `matrix dissimilarity` (`matching`, `Jaccard`, `Russell`, `Hamann`, `Dice`, `antiDice`, `Sneath`, `Rogers`, `Ochiai`, `Yule`, `Anderberg`, `Kulczynski`, `Pearson`, `Gower2`); `cluster kmeans` and `kmedians` with `keepcenters`
- [x] `manovatest` (tests of terms after `manova`, with error terms after `/` and `ytransform()`)
- [x] `camat`, `cabiplot`, `caprojection`, `mcaplot` and `mcaprojection`
- [x] `mdslong`, `mdsconfig` and `mdsshepard`; `procoverlay`; `procrustes` with `transform(unrestricted)` and `transform(oblique)`
- Checked against the Stata 14 log: case `1917` matches except for the iteration log of `mdslong, method(modern)` (Stata's optimizer is undocumented; see `compat/known_differences.txt`). Found and fixed on the way: the `_ord` and `_hgt` variables of a hierarchical cluster analysis are by dendrogram position (the id of the k-th leaf and the height of the merge that joins it to the next one), the cluster holding the smallest observation goes first, binary similarity measures keep their heights as similarities and `cluster generate` refuses a number of groups that ties make ambiguous; `clustermat` reports the new observations and the cluster name and puts the label variable last; `manovatest` with several terms tests them together in the base-level coding, with an error term of fewer degrees of freedom than dependent variables it stops, and it has `showorder`; `list` abbreviates a name so that it matches only that variable (`_clus~id`)

## Phase 32: models, part 3 ✅

- [x] `tpoisson` and `tnbreg` with `ul()` and limits given by variables; `tnbreg` and `nbreg` with `dispersion(constant)` (`nbreg` took the option and fitted the mean-dispersion model)
- [x] `glogit` with `or`; `glm` and `binreg` with `mu()`
- [x] `offset()` in `intreg`, `truncreg`, `heckman`, `heckprobit` and `heckoprobit` (accepted and ignored before), `offset1()` and `offset2()` in `biprobit`; the offset row in the coefficient tables of these models and of the count models; `heckman` with `mills()` (two-step and maximum likelihood)
- [x] `sureg` and `reg3` with `constraints()`
- [x] `ivregress gmm` with `wmatrix(hac kernel #)`, `center` and `igmm`
- [x] `frontier` with `uhet()`, `vhet()` and `cm()`
- Checked against the Stata 14 log: case `1918` (the differences and the reason for each are in `compat/known_differences.txt`). Stata 14 has no `ul()` in `tpoisson` and `tnbreg` (it came later): OpenDTA keeps it as an extension. Found and fixed on the way: the negative binomial likelihoods sum ln Gamma(y + r) - ln Gamma(r) term by term when the dispersion goes to zero (the difference of the two lost digits and put the NB1 likelihood above the Poisson one); the window probability of a count with an upper limit is summed term by term; `tnbreg` and `nbreg` with `dispersion(constant)` start the full model from the constant-only fit and `tnbreg` stores Stata's e(); `sureg` and `reg3` with constraints take the first residual covariance from the constrained equation-by-equation fit and store Stata's e() (`Sigma`, `Cns`, the statistics of each equation), and constrained `2sls` reports z statistics; `ivregress gmm` with a HAC weight matrix: the quadratic-spectral bandwidth is # + 1, the footnote of the HAC VCE, e(J), e(hac_lag), e(S) and e(W), and the log of `igmm`; `frontier` with `uhet()` or `vhet()` reports both variances as equations, and with `cm()` the parameters `/lnsigma2` and `/ilgtgamma`

## Phase 33: power and sample size, part 2 ✅

- [x] `power twoway` (noncentral F, row, column and interaction effects, `cellweights()`, `showcellsizes`), `power repeated` (between, within and between-within effects, compound symmetry or `covmatrix()`, Greenhouse-Geisser correction), `power cmh`, `power mcc` (with `m()` and `compare`), `power trend`, `power cox`, `power exponential` and `power logrank`; sample size, power and effect size, lists of values in tables
- Statistics in `core/power_suites.py`; command in `commands/power3.py`. Case `1919` has no Stata reference yet: what the manual does not settle (return-list order, solver values, table layouts, the definition of `delta` in `cmh` and `mcc`, the `m` of the between-within test, rounding of events before subjects) is marked `VERIFY` in the code and runs under `capture noisily`.
- Not available yet in these methods: `graph()`, `saving()`, `show`, `st1()`, `n1()` and `n2()`

## Distribution and documentation

- [x] MIT License and license notices for third-party components
- [x] macOS download (Apple Silicon and Intel), with the Release built by GitHub Actions for every `v*` tag
- [x] All documentation published on GitHub in English, with a Portuguese copy (`*.pt-BR.md`) and an **English | Português** selector at the top of each page; code comments and commit messages in English
- [x] Update check and installation from inside the app (GitHub Releases, verified download, daily check that can be turned off)
- [x] Windows download (64-bit, Intel/AMD): per-user installer and portable `.zip`, built by the same workflow
- [x] Version 1.0: short README with a screenshot, documentation in `docs/` (getting started, features, testing, development), `CONTRIBUTING.md`, issue templates (difference from Stata, bug, command request) and a pull request template; the release workflow checks that the tag matches the version

## Improvements over Stata (throughout the phases)

- [x] Error explanations: after the message and before r(#), one line says what went wrong in the context of the command and another suggests how to fix it (similar names, similar files in the folder, number of observations…). In English or Portuguese, following the interface language; `set hints off` turns them off
- [x] Explanations for beginners: what went wrong with the data at hand (missing values per variable, string or numeric variables, values of the outcome, contradictions of assert, matrix sizes), why the rule exists, the corrected command (= for ==, egen functions in generate, quotes for text and value labels, bysort, names without spaces or accents...) and a clickable help link; new help pages for operators, functions, missing values and by
- [x] `set numerics precise`: QR in `regress`, `D.`/`S.` in double, a final Newton step in `ml`/`optimize()`, Richardson extrapolation for numerical derivatives and new variables in double (the default `stata` reproduces Stata's rounding)
- [ ] Statistical options beyond Stata: `vce(cr2)` with Satterthwaite degrees of freedom, *wild cluster bootstrap*, `logit, firth`
- [x] Superscript exponents in displayed scientific notation (`1.0000×10¹⁰`, `2.5×10⁻⁵`), keeping the format width; `set superscript off` (or Help → Superscript Exponents) restores Stata's `e+10`. `string()`, macros, `file write`, exports and Mata's `sprintf()`/`strofreal()` always use Stata's form
- [x] Review can list the commands run by do-files (`do`, `run`, `include` and the Do and Run buttons of the Do-file Editor) under the do line, nested do-files included and failed commands in red, with the editor's selections named after their file (Settings > General > Session; off by default, as in Stata, where Review shows only `do "<temporary file>"`)

## Missing for 1.0

Commands and options of Stata 14 that OpenDTA does not have yet, by area. All of them are needed for 1.0; the order is the one of the blocks in [docs/next-steps.md](docs/next-steps.md), adjusted by the requests opened with the [Command request](https://github.com/yaskepalhares/opendta/issues/new?template=command-request.md) template.

- **Data management**: `icd9p` with the ICD-9-CM procedure code list.
- **Models**: `ivprobit` and `ivtobit` with `twostep` and with more than one endogenous regressor; `ivpoisson cfunction` with more than one endogenous regressor; `ivregress gmm` with `wmatrix(hac kernel opt)`; `biprobit` with `partial` and `constraints()`; `asmprobit` and `asroprobit` with `altwise`, unbalanced alternatives and other covariance structures; `slogit` with `dim()` above 1; `exlogistic` and `expoisson` with `terms`, `group()`, `binomial()`, `condvars()`, `saving()`, `exposure()` and `offset()`; the program-based forms of `nl`, `nlsur` and `gmm`, and `gmm` with `derivatives()`, `parameters()` and `igmm`; `frontier` with `constraints()`; `teffects` with multivalued treatments, other outcome models and `predict`; `stteffects wra`; `eteffects` with other outcome and treatment models; `etregress` and `etpoisson` with `cfunc`.
- **Panel data**: `xtunitroot` with `llc` and `ips` (they need the tabulated moments of the original papers), `trend` for `ht`, `lags()` for `breitung` and the automatic lag choices of `kernel()`; `xtregar` with `lbi`; `xtdata` with `fa`; `xtabond`, `xtdpdsys` and `xtdpd` with `teffects`, `small`, `diffvars()`, `collapse`, `h()` and `vce(cluster)` (some of these may be options of the user-written `xtabond2` rather than of Stata 14: to check in the manual).
- **Time series**: `vec` with `bconstraints()`, `aconstraints()` and `sindicators()`, and the Osterwald-Lenum critical values for the trends other than `constant` and `rconstant`; `svar` with `lrconstraints()`; `fcast` with `bs`, `bsp`, `exog` and non-consecutive lags; `arfima` with `mpl` and `smemory`; `ucm` with `constraints()`; `sspace` with `covstate()` and `covobserved()` other than `diagonal`, constraints on the variances and states with lags above one; `dfactor` with more than one factor, `ma()` and error autoregressions; `mgarch vcc`, `distribution(t)`, `het()`, `constraints()` and `vce(robust)`; `predict` with the states, components and `smethod()` and `estat` after `ucm`, `sspace` and `dfactor`; `predict` and `estat` after `mgarch`.
- **Multivariate statistics**: `manovatest` with `test()`; `pca` with `vce(normal)`, `means` and supplementary variables; `factor` with `estat factors` and `estat rotatecompare`; `rotate` with `target()`, `matrix()` and the other criteria (`oblimax`, `bentler`, `cf()`…); `manova` with error terms; `canon` with `estat`; `mvtest` with `compound`, `block` and `pairwise` hypotheses and `normality` with `univariate`, `bivariate`, `by()` and Doornik-Hansen; `discrim` with `ties()`, `priors(matname)`, weights and `lootable` after `logistic`, `estat correlations`, `estat grsummarize` and `estat manova`; `ca` and `mca` with supplementary points, `method(indicator)` and `method(johnson)` and the other `estat` subcommands; `mds` with `method(nonmetric)`, `loss()`, `transform()` and `addconstant`; `mdslong` with `s2d()`.
- **Multilevel models**: `intmethod(laplace)` with more than one level; crossed random effects (`_all:`) in the `me` commands.
- **Other suites**: `power` with `graph()`, `saving()`, `show`, `st1()` and `n1()`/`n2()` in `twoway`, `repeated`, `cmh`, `mcc`, `trend`, `cox`, `exponential` and `logrank`; `bayes` with the estimation commands other than `regress`, `logit`, `probit` and `poisson`; `irt` with `group()`; `gsem` with observation-level and cluster-level latent variables in the same model; `svy:` with `ivprobit`, `ivtobit` and the other estimators that do not take weights yet; `e(b)` of `sem`, `gsem` and `irt` with the constrained parameters (and `e(Cns)`); `svy:` storing the `e()` results of the estimator.
- **Programming**: Stata's class system (Mata classes already work).

Left out of 1.0 on purpose, as known limitations: `import haver`, which reads Haver Analytics databases only through Haver's own software (Windows and a paid license), and loading compiled Stata plugins (`program name, plugin` and `plugin call`), which are machine code written against Stata's internal interface; the same computation can be written in a do-file, an ado-file or Mata.

## After 1.0

What comes once every command of Stata 14 is covered, in priority order. These are improvements over Stata for population-based health surveys and epidemiology, plus later Stata versions. Most of them have no Stata 14 counterpart, so their reference is the literature or R (`survey`) rather than a Stata log; each one states how it is validated. The statistics go in `core/`, apart from the syntax (see [docs/next-steps.md](docs/next-steps.md)).

**P0: complex surveys and epidemiology**

1. Weight calibration: post-stratification, raking, linear and logit calibration, trimming, and weight diagnostics (Kish effective sample size, variability of weights).
2. Prevalence confidence intervals for rare or extreme estimates (Clopper-Pearson and Wilson with effective sample size) and a reliability report (CV, degrees of freedom, effective sample size, PSUs per domain), as a warning in the `svy:` output.
3. Ready-made results tables in one command: a weighted Table 1 (`n`, `%`, CI, p-value, `by()`), regressions side by side (`esttab`-like) with prevalence ratios and CIs in English or Portuguese, exported to Word, Excel and LaTeX (reusing `putexcel`).
4. Prevalence ratios and adjusted prevalence: `binreg` (block 10) with `svy:`, and `margins`/`predict` with CI for adjusted prevalence after `poisson`, `glm` and `binreg`.
5. The options already listed under [Improvements over Stata](#improvements-over-stata-throughout-the-phases): `vce(cr2)` with Satterthwaite degrees of freedom, *wild cluster bootstrap* and `logit, firth`.
6. Pre-run checks and hints for `svy:`: single-PSU strata, `if` used instead of `subpop()`, missing `svyset` or weights (extends `core/hints.py`).

**P1: survey workflow**

7. Importer for the survey dictionaries of IBGE (PNS, POF, PNAD) and Vigitel layouts: generates the `infile` and the value labels.
8. Sample size for surveys: design effect, non-response, optimal allocation by stratum, clusters (ICC).
9. Health inequality indices (slope and relative index of inequality, concentration index) with `svy:`.
10. Non-response: inverse probability weighting, and `mi` consistent with `svyset`.
11. Scales and indices with weights: Cronbach's alpha and omega, weighted factor analysis and IRT (after blocks 8 and 11), and food-security and diet scores (EBIA, NOVA, diet adequacy).
12. Direct and indirect standardization with the survey design, and trends across surveys (joinpoint/APC), with tests of differences between domains and multiple-comparison corrections.

**P2: reporting, graphics and tools**

13. Reproducible reports: a do-file that writes a Markdown/HTML/docx document with commands, tables and graphs.
14. Graphs: export and import of Plotly JSON (a documented subset of graph types); integrity tests for the SVG, PNG, JPEG and PDF exports; SVG with editable text. OpenDTA keeps its own `.gph` format: reading Stata's `.gph` stays out of scope (see below).
15. Python inside do-files (`python:` blocks, with pandas and scikit-learn).
16. Parallel bootstrap, jackknife and replicate weights; the environment and the seed recorded with each result.

**P3: later**

17. Small-area estimation (Fay-Herriot area models).
18. The commands and changes of each later Stata version, one at a time, up to Stata 19.
19. A version that runs in the web browser (see Long term).

## Long term

For the record, the order after 1.0:

1. Finish every command of the target version (the blocks in [docs/next-steps.md](docs/next-steps.md)) and release 1.0. Nothing from [After 1.0](#after-10) comes first.
2. The priorities of [After 1.0](#after-10), then the commands and changes of each later Stata version, one version at a time, up to Stata 19.
3. A version that runs in the web browser, without installing anything: the statistical core (Python, NumPy and SciPy) runs in WebAssembly through Pyodide. It builds on the separation between the statistical core and the Stata layer that new commands already follow.

Out of scope: reading the `.gph` files saved by Stata, whose format has no public specification ([development policy](docs/clean-room-policy.md)); OpenDTA saves graphs in its own `.gph` format.

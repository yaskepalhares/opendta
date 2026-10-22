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
- [x] Mata's `moptimize()` (lf, d0–d2, lf0–lf2 and gf0–gf2 evaluators; oim, opg and robust variance; case `0408` awaits the Stata run)
- Checked against the Stata logs: cases `0501`–`0513`

## Phase 6: complex surveys ✅

- [x] `svyset` (strata, PSUs, weights, FPC), `svydescribe`
- [x] `svy:` prefix with `mean`, `proportion`, `ratio`, `total`, `tabulate`, `regress`, `logit`, `logistic`, `poisson`; `subpop()` and `over()`
- [x] Taylor linearization; Rao–Scott correction in two-way `svy: tabulate`
- [x] Several stages, BRR, jackknife, bootstrap and SDR; `estat effects`; `svy:` also with `ologit`, `oprobit`, `mlogit`, `nbreg`, `glm` and `tobit`
- Checked against the Stata logs: cases `0601` and `0602`; `0603` awaits the next run

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
- [x] Mata classes (`class`, `extends`, virtual methods, `super`, `this`, static members, constructors; case `0409` awaits the Stata run)
- [x] Survival analysis: `stset` (with `id()`, `time0()`, `origin()`, `enter()`, `exit()`, `scale()`), `st`, `streset`, `stdescribe`, `stsum`, `sts list`/`generate`/`test`/`graph` (Kaplan-Meier, Nelson-Aalen, log-rank, Wilcoxon, Tarone-Ware, Peto-Peto), `stcox` (Breslow and Efron ties, strata, robust), `estat phtest`, `streg` (exponential, Weibull, Gompertz, lognormal, loglogistic) and their `predict` options; cases `0901`–`0902` await the Stata run
- [x] Multilevel models: `mixed` (ML and REML, nested levels, random slopes, identity/independent/exchangeable/unstructured covariances; `estat icc`, `group`, `recovariance`; `predict`) and `melogit`, `meprobit`, `mepoisson` (adaptive Gauss-Hermite quadrature); cases `1001`–`1002` await the Stata run
- [x] Time series: `arima` (exact ML with the Kalman filter, OPG variance), `var`, `varsoc`, `vargranger`, `dfuller` (Fuller critical values, MacKinnon p-values), `corrgram`, `newey`; case `1101` awaits the Stata run
- [x] Multiple imputation: `mi set` (mlong, flong), `register`, `describe`, `impute` (`regress`, `pmm`, `logit`, `ologit`, `mlogit`, `poisson`, `chained`), `estimate` (Rubin's rules, Barnard-Rubin df), `passive`, `xeq`, `extract`, `misstable`; case `1201` (imputed values differ from Stata's: known difference)
- [x] `sem` (path models, CFA and latent variables by ML, standardized solution, `estat gof`, `predict`), `power` (`onemean`, `twomeans`, `pairedmeans`, `oneproportion`, `twoproportions`, `onecorrelation`, `twocorrelations`) and `irf` (`create`, `table`, `graph`, `describe`, `rename`, `drop` after `var`); cases `1301`–`1303` await the Stata run
- [x] Completed commands: `anova` (sequential sums of squares, nested terms, error terms, repeated measures), `xtgee` and the `pa` estimator of `xtreg`/`xtlogit`/`xtpoisson`, `xtreg, be` and `mle`, two-way `tabulate, summarize()`, `misstable tree`/`nested`, `ci`/`cii variances`, two-sample `prtest`/`prtesti`, `reshape i/j/xij/error`, the old `merge` syntax, `sort ... in`, `list, sepby()`, `display _request()`, `export excel` to `.xls`, `set rng kiss32`; sequential `replace` with time-series operators (`replace y = 2*L.y`); case `0805` awaits the Stata run
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
- Cases `1401` and `1402` await the Stata run

## Phase 11: data management ✅

- [x] Variables and observations: `clonevar`, `insobs`, `range`, `separate`, `ipolate`, `expandcl`, `stack`, `xpose`, `drawnorm`, `corr2data`
- [x] Describing and comparing data: `ds`, `lookfor`, `inspect`, `labelbook`, `compare`, `cf`, `datasignature`, `snapshot`
- [x] Files and the system: `copy`, `mkdir`, `rmdir`, `dir`/`ls`, `shell` (and `!`), `filefilter`, `zipfile`, `unzipfile`, `checksum`, `outfile`, `unicode`
- [x] `statsby`, `putmata`; `sysuse` (OpenDTA's own example datasets) and `webuse` (downloads on demand)
- [x] `icd9` and `icd10` with the public-domain ICD-9-CM and ICD-10-CM code lists of the U.S. CMS (check, clean, generate, lookup, search); `icd9p` by format
- [x] `odbc` through the optional `pyodbc` package; business calendars (`bcal`, the `%tb` format, `bofd()` and `dofb()`)
- Case `1501` awaits the Stata run

## Phase 12: epidemiology, finished ✅

- [x] `cs`, `cc` and `ir` with `by()`: directly pooled estimates (`pool`) with their test of homogeneity, internally and externally standardized estimates (`istandard`, `estandard`) or with weights from a variable (`standard()`), standardized risk and rate differences (`rd`, `ird`), binomial frequency records (`binomial()`), and the Breslow-Day and Tarone tests of `cc` (`bd`, `tarone`)
- [x] `mhodds` with `by()`: the odds ratio within each group, the estimate controlling for the groups and a test of homogeneity
- [x] `roctab` with Bamber's standard error (`bamber`), exact binomial limits (`binomial`) and the Lorenz curve with the Gini and Pietra indices (`lorenz`); `roccomp` with binormal areas (`binormal`)
- [x] `rocreg`: nonparametric ROC estimation with covariate-adjusted controls (stratified or linear, empirical or normal percentile values) and its summaries (AUC, ROC(t), inverse ROC, partial AUC), and the parametric probit GLM model with ROC covariates; bootstrap standard errors with normal, percentile and bias-corrected intervals
- Case `1601` awaits the Stata run

## Phase 13: data management, options ✅

- [x] `file` with binary read and write (`%1b`, `%2b`, `%4b`, their signed and unsigned forms, `%4z`, `%8z`, `%#s`), `file seek`, `file set byteorder`, `file query`, `read write` and the `byteorder()` function
- [x] `infile` with `byvariable()`; `putexcel` to `.xls`
- [x] `import sas` (with `bcat()`), `import spss` and `import sasxport`, and a basic `export sasxport`, through the optional `pyreadstat` package
- [x] `describe` with `numbers`, `varlist`, `replace` and `using`; `list` with `header()`, `nodotz`, `mean`, `sum`, `N`, `labvar()` and `nocompress`
- Case `1701` awaits the Stata run

## Phase 14: programming ✅

- [x] `unab`, `tsunab`, `fvunab`, `unabcmd`, `fvexpand`, `fvrevar`, `tsrevar` and `_rmcoll`, with factor-variable and time-series operators (specific forms such as `i(2,3)bn.rep78` included)
- [x] `mark`, `markin` and `markout` (with `sysmissok`); `marksample` and `syntax varlist(fv ts)` with factor variables and time-series operators
- [x] `postfile`, `post`, `postclose` and `postutil`
- [x] `matlist` (column and row specifications, borders, formats, titles) and `tabdisp` (one to three variables, `by()`, totals, `concise`, `dotz`)
- [x] `program, sortpreserve`; `program name, plugin` (explains that Stata plugins cannot be loaded)
- [x] `pause`, `sleep`, `window` (in batch mode the dialogs act as cancelled) and `viewsource` (shows OpenDTA's own source)
- [x] `signestimationsample`, `checkestimationsample`, `makecns` with a constraint matrix and `matcproc`
- Case `1801` awaits the Stata run

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
- Checked against the Stata log: case `1902` awaits the Stata run

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

## Phase 20: multivariate statistics

- [x] Principal components and factor analysis: `pca`, `factor` (`pf`, `pcf`, `ipf`, `ml`), `rotate` (`varimax`, `quartimax`, `equamax`, `parsimax`, `promax`, `oblimin`, `quartimin`), `alpha`, `screeplot`, `loadingplot` and `scoreplot`, with `estat` (`kmo`, `smc`, `anti`, `loadings`, `correlation`, `residuals`, `summarize`, `structure`, `common`) and `predict`
- [x] Multivariate tests and models: `hotelling`, `mvreg`, `manova` (Wilks, Pillai, Lawley-Hotelling and Roy), `canon` and `mvtest` (`means`, `covariances`, `correlations`, `normality`)
- [x] Cluster analysis: `cluster` (`kmeans`, `kmedians`, the seven hierarchical linkages, `generate`, `dendrogram`, `stop` with `calinski` and `duda`, `list`, `dir`, `drop`, `rename`, `notes`) and `matrix dissimilarity`
- [x] Discriminant analysis: `discrim` (`lda`, `qda`, `knn`, `logistic`) and `candisc`, with `estat` (`classtable`, `errorrate`, `list`, `grdistances`, `covariance`) and `predict`
- [x] Correspondence analysis and scaling: `ca`, `mca`, `mds` and `mdsmat` (classical and modern), orthogonal `procrustes` and `biplot`, with `estat` and `predict` where they apply
- Checked against the Stata 14 log: case `1906` matches except for the differences in `compat/known_differences.txt` (the iteration paths of `factor, ml`, `discrim logistic` and `mds, method(modern)`, and the maximum-likelihood factor numbers that follow from where Stata's `ml` stops)

## Phase 21: models, part 1 (binary and count, censored and selection, instrumental variables)

- [x] Binary and fractional response: `cloglog`, `scobit`, `hetprobit` (`het()`), `fracreg` (`logit`, `probit`, `cloglog`), `binreg` (`or`, `rr`, `hr`, `rd`) and `glogit`
- [x] Count data: `zip` and `zinb` (`inflate()`, `probit`, `vuong`, `irr`), `tpoisson`, `tnbreg` and `gnbreg` (`lnalpha()`)
- [x] Censored, truncated and selection models: `intreg`, `truncreg` (`ll()`, `ul()`), `heckman` (maximum likelihood and `twostep`), `heckprobit`, `heckoprobit` and `biprobit` (shared varlist and `(y1 = x) (y2 = z)` forms), with the bivariate normal distribution function in `core/bvnorm.py` (Plackett's identity, exact derivatives)
- [x] Instrumental variables: `ivregress liml` (k-class) and `gmm` (two-step, `wmatrix()`), `ivprobit` and `ivtobit` (maximum likelihood, Wald test of exogeneity), `ivpoisson gmm` (additive and `multiplicative` moments, `onestep`) and `cfunction`
- [x] A multi-equation maximum-likelihood engine (`estimation/mleq.py`) with the log likelihoods apart from the commands (`estimation/likelihoods.py`), reused by the choice models of the next part
- Checked against the Stata 14 log: case `1907` matches, except for the differences in `compat/known_differences.txt` (the iteration paths and the stopping point of each optimizer; blocks whose data put a parameter at the boundary, so that Stata's estimates depend on its path: `scobit`, `zinb`, `tnbreg`, `gnbreg`; `binreg, hr`, where Stata's IRLS does not converge on the case data; `glogit`, run on other grouped counts than Stata's; the robust variance of `ivregress liml` and of `ivpoisson cfunction`). Found and fixed on the way: `binreg` is IRLS with the expected-information variance, `ivpoisson cfunction` uses the multiplicative residual, `ivregress gmm` with a robust weight matrix reports the sandwich variance, and the constant-only models of `tpoisson` and `gnbreg` follow Stata's reference points
- Not verified against Stata yet: `glogit` (data of the case), `fracreg cloglog`, `heckoprobit` and the options not exercised by case `1907`

## Phase 22: models, part 2 (choice, exact, systems and treatment effects)

- [x] Systems and constraints: `constraint`, `cnsreg`, `sureg` (including `isure`), `reg3` (`2sls`, `3sls`, `ols`, `sure`), `eivreg`, `vwls` and `rreg` (Huber then biweight), with the algebra in `core/syseq.py`
- [x] Choice models: `clogit` (any number of chosen alternatives per group), `rologit` (`exactp`, `breslow`, `efron` and `none` for ties), `asclogit`, `slogit`, `mprobit`, `asmprobit` and `asroprobit` (GHK simulator with Hammersley draws), with the likelihoods in `core/choice.py`
- [x] Exact methods: `exlogistic` and `expoisson` (conditional maximum likelihood, median unbiased estimates, exact confidence intervals and p-values, from the exact distribution of the sufficient statistics in `core/exact.py`)
- [x] Nonlinear models: `nl`, `nlsur` (NLS, FGNLS and IFGNLS), `gmm` (one- and two-step, robust and cluster weight matrices) with substitutable expressions, and `frontier` (`hnormal`, `exponential`, `tnormal`)
- [x] Treatment effects: `teffects ra`, `ipw`, `ipwra`, `aipw`, `nnmatch` and `psmatch` (Abadie-Imbens variance), `eteffects` (control function), `stteffects ra`, `ipw` and `ipwra`, `etregress` (maximum likelihood and `twostep`) and `etpoisson` (Gauss-Hermite), with the estimators in `core/treatfx.py`
- Checked against the Stata 14 log: case `1908` (see `compat/known_differences.txt` for what does not match and why)
- Not verified against Stata yet: the table layouts of every command of this part (marked `VERIFY` in the code), the variance of `eivreg`, the weights of `stteffects` and the options the case does not exercise

## Phase 23: other suites (sem, gsem, irt, Bayesian estimation, power, mi impute and svy)

- [x] `sem` with `method(mlmv)` (normal likelihood by missing-data pattern, EM starting values in `core/mvnmiss.py`) and `method(adf)`, `group()` and `ginvariant()`, path labels (`@name`) and the `estat` subcommands `eqgof`, `eqtest`, `framework`, `ggof`, `ginvariant`, `mindices`, `residuals`, `scoretests`, `stable`, `stdize:`, `summarize` and `teffects`
- [x] `gsem` with the `gaussian`, `bernoulli`, `binomial`, `poisson`, `nbinomial`, `gamma` and `ordinal` families, observation-level and cluster-level (`M[id]`) latent variables, adaptive Gauss-Hermite quadrature and `predict`
- [x] `irt` (`1pl`, `2pl`, `3pl`, `grm`, `nrm`, `pcm`, `gpcm`, `rsm` and `hybrid`) on top of `gsem`, with `estat report` and `irtgraph` (`icc`, `iif`, `tcc`, `tif`)
- [x] `bayesmh` (adaptive random-walk Metropolis-Hastings in `core/mcmc.py`, eight likelihoods and ten priors), the `bayes` prefix for `regress`, `logit`, `probit` and `poisson`, `bayesstats` (`summary`, `ess`, `ic`), `bayestest` (`interval`, `model`) and `bayesgraph`
- [x] `power` with `onevariance`, `twovariances`, `pairedproportions` and `oneway`, the `lrchi2` and `fisher` tests of `twoproportions` and the `graph` option
- [x] `mi impute mvn` (data augmentation), `monotone`, `truncreg`, `intreg` and `nbreg`
- [x] `svy:` with `cloglog`, `intreg`, `truncreg`, `heckman`, `heckprobit`, `zip`, `zinb`, `gnbreg`, `hetprobit`, `scobit`, `biprobit`, `fracreg`, `slogit`, `mprobit`, `clogit` and `ivregress`
- Not verified against Stata yet: case `1909` waits for its reference log; the layouts of every command of this phase are marked `VERIFY` in the code. The draws of `bayesmh` and `mi impute` come from OpenDTA's own generator and do not reproduce Stata's

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

## Missing for 1.0

Commands and options of Stata 14 that OpenDTA does not have yet, by area. All of them are needed for 1.0; the order is the one of the blocks in [docs/next-steps.md](docs/next-steps.md), adjusted by the requests opened with the [Command request](https://github.com/yaskepalhares/opendta/issues/new?template=command-request.md) template.

- **Data management**: `table` with three or more variables; `import haver`; `import sasxport` with `describe`, `member()` and the value labels of `formats.xpf`; `export sasxport` with `rename` and `vallabfile()`; `list` with `string()` and `display`; `icd9p` with the ICD-9-CM procedure code list.
- **Epidemiology**: `rocreg` with the `ml` and `cluster()` options.
- **Base statistics**: `margins` after more commands; `ml` with the `lf0`, `lf1`, `lf2` and `gf` methods, `constraints()` and `svy`; `boxcox` with `model(lhsrhs)`, `rhsonly` and `notrans()`; `icc` with `testvalue()` and unbalanced data; `loneway` with `median` and `exact`; `esize` with `unequal` and `welch`; `pwmean` with `unequal`; `ameans` with weights; `lv` with `generate()` and `tail()`; `stem` with `width()`; `cusum` with the test statistics.
- **Models**: `ivprobit` and `ivtobit` with `twostep` and with more than one endogenous regressor; `ivpoisson cfunction` with more than one endogenous regressor; `ivregress gmm` with `wmatrix(hac)`, `center` and `igmm`; `glogit` with `or`; `tpoisson` and `tnbreg` with `ul()`; `tnbreg` with `dispersion(constant)`; `binreg` with `mu()`; `biprobit` with `partial`, `offset1()`, `offset2()` and `constraints()`; `heckman` and `heckprobit` with weights in `twostep` and with `mills()`; `intreg` and `truncreg` with `offset()`; `sureg` and `reg3` with `constraints()`; `asmprobit` and `asroprobit` with `altwise`, unbalanced alternatives and other covariance structures; `slogit` with `dim()` above 1; `exlogistic` and `expoisson` with `terms`, `group()`, `binomial()`, `condvars()`, `saving()`, `exposure()` and `offset()`; the program-based forms of `nl`, `nlsur` and `gmm`, and `gmm` with `derivatives()`, `parameters()` and `igmm`; `frontier` with `uhet()`, `vhet()`, `cm()` and `constraints()`; `teffects` with multivalued treatments, other outcome models and `predict`; `stteffects wra`; `eteffects` with other outcome and treatment models; `etregress` and `etpoisson` with `cfunc`.
- **Panel data**: `xtunitroot` with `llc` and `ips` (they need the tabulated moments of the original papers), `trend` for `ht`, `lags()` for `breitung`, `demean`, `drift` and `pperron` for `fisher` and the kernel options of `hadri`; `xtregar` with `lbi`; `xtrc` with `betas`; `xtivreg` with `first` and `regress`; `xtdata` with `re` and `fa`; `xtabond`, `xtdpdsys` and `xtdpd` with `teffects`, `small`, `diffvars()`, `collapse`, `hascons`, `h()`, `artests()` and `vce(cluster)`.
- **Time series**: `vec` with `bconstraints()`, `aconstraints()` and `sindicators()`, and the Osterwald-Lenum critical values for the trends other than `constant` and `rconstant`; `svar` with `aconstraints()`, `bconstraints()` and `lrconstraints()`; `fcast` with `bs`, `bsp`, `exog` and non-consecutive lags; `arch` with `constraints()` and `from()`; `arfima` with `mpl`, `smemory` and `vce(robust)`; `ucm` with `constraints()` and the `frequency()` suboption of `cycle()`; `sspace` with `covstate()` and `covobserved()` other than `diagonal`, `constraints()` and states with lags above one; `dfactor` with more than one factor, `ma()` and error autoregressions; `mgarch vcc`, `distribution(t)`, `het()`, `constraints()` and `vce(robust)`; `prais` with `vce(robust)`; `predict` and `estat` after `ucm`, `sspace`, `dfactor` and `mgarch`; `psdensity` after `ucm`.
- **Multivariate statistics**: `clustermat`, `manovatest`, `camat`, `cabiplot`, `cagraph`, `caprojection`, `mcaplot`, `mcaprojection`, `mdsconfig`, `mdsshepard`, `mdslong`, `procoverlay`; `pca` with `vce(normal)`, `means` and supplementary variables; `factor` with `estat factors` and `estat rotatecompare`; `rotate` with `target()`, `matrix()` and the other criteria (`oblimax`, `bentler`, `cf()`…); `manova` with error terms; `canon` with `estat`; `mvtest` with `compound`, `block` and `pairwise` hypotheses and `normality` with `univariate`, `bivariate`, `by()` and Doornik-Hansen; `cluster` with the binary similarity measures and `keepcenters`; `discrim` with `ties()`, `priors(matname)`, weights and `lootable` after `logistic`, `estat correlations`, `estat grsummarize` and `estat manova`; `ca` and `mca` with supplementary points, `method(indicator)` and `method(johnson)` and the other `estat` subcommands; `mds` with `method(nonmetric)`, `loss()`, `transform()` and `addconstant`; `procrustes` with `transform(oblique)` and `transform(unrestricted)`.
- **Multilevel models**: `meglm` and the other `me` commands with more than one level of random effects, weights or `vce(robust)`; `mestreg` with weights, `vce(robust)` and the `median` and `mean` options of `predict` (and `hazard` after the lognormal, loglogistic and gamma distributions); `predict, reffects` after `mixed` with `R.` equations.
- **Other suites**: `power` with the `twoway`, `repeated`, `cmh`, `mcc`, `trend`, `cox`, `exponential` and `logrank` methods; `bayes` with the estimation commands other than `regress`, `logit`, `probit` and `poisson`; `irt` with `group()`; `gsem` with observation-level and cluster-level latent variables in the same model; `svy:` with `ivprobit`, `ivtobit` and the other estimators that do not take weights yet.
- **Programming**: `makecns` with a constraint list and `displaycns` (they need `constraint`, from the Models area); `matlist` with the equation options; `window` dialogs and menus in the graphical interface; `plugin` loading (compiled Stata plugins cannot run in OpenDTA); Stata's class system (Mata classes already work).

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

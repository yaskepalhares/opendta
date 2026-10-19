# Next steps: closing the roadmap

**English** | [Português](next-steps.pt-BR.md)

A working plan for the commands still missing in [After 1.0](../ROADMAP.md#after-10), written so that each block can be done in a fresh session without rediscovering the code. Version 1.0 is released when every block is done; after that come improvements to what exists and newer Stata versions.

## How to run a block

Start a new session (a long conversation costs more on every step) with an instruction such as:

> Implement the **Programming** block as described in docs/next-steps.md. Open a pull request and merge it when CI passes.

Each block follows the same steps:

1. Branch from `main` (`block-<name>`).
2. For each command: the statistics in `src/opendta/core/` (see the next section), the command in `src/opendta/commands/` (new module registered in `commands/__init__.py`), its help page, its tests.
3. A compat case `compat/do/NNNN_<block>.do` with the commands and `return list`/`ereturn list`; it has no reference until the Stata run. Run under `capture noisily` every command whose behavior in Stata the manual does not document (and do the same in any old case you change), so that a single Stata run shows the whole case instead of stopping at the first error.
4. `ROADMAP.md` and `ROADMAP.pt-BR.md`: a new phase listing what was done; in *After 1.0*, leave only what is missing. Each clause of an *After 1.0* line that lists options must start with the command they belong to (`` `cs` with the `pool` option``), so the coverage map does not count options as commands.
5. `tools/roadmap_coverage.py`: add the new module to `MODULE_AREA`.
6. `docs/features.md` and `docs/features.pt-BR.md`: one line for the block.
7. The Stata run, once, for everything the block needs (cases without a reference included): ask the user to run `compat/gerar_esperados.do`, never a loose block of commands. That file already sets `dp period` and `linesize 255`, erases the `odta_*` leftovers and runs `mata clear`. Before asking, check with git that the branch is on GitHub and tell the user to `git pull` it. If a run fails for an environment reason (decimal comma, a file name, leftovers of another run), say so before asking for another.
8. `python tools/compare.py` with no new differences; open the PR; merge when the tests pass on Linux, macOS and Windows. The coverage page updates itself after the merge.

Notes for sessions that run in the cloud:

- The work is done in a clone in the cloud, which reaches GitHub; the user's computer does not. After each push, tell the user to `git pull` the branch. Use `git push -u origin <branch>` on the first push, so that git knows the remote branch.
- `pytest` may not install there (no access to PyPI): the tests that need optional packages (`pyreadstat`, `xlrd`, `xlwt`) run only in CI. Wait for CI and read the result before asking for the Stata run.
- The CI logs cannot be downloaded from the cloud, but the annotations of each job show the `tools/compare.py` differences (`gh api repos/<owner>/<repo>/check-runs/<id>/annotations`). `gh pr` fails because GraphQL is blocked: open and merge the pull request with `gh api` (REST), using a merge commit.
- A reference log from the Stata run is never edited by hand: when a case changes, the Stata run is repeated.

Next free numbers: compat cases from **1801**, roadmap phases from **14**.

## Rule: statistics in the core, Stata in the command

The computation goes in plain functions in `src/opendta/core/` (or `src/opendta/estimation/` for models): they take arrays and numbers and return numbers, arrays or small result objects, without a `Session`, without reading variables and without printing. The command module parses the syntax, reads the data, calls those functions, writes Stata's output and stores `r()`/`e()`. Example: `core/epi.py` (formulas) and `commands/epitab.py` (tables). Tests check the core functions against closed formulas or SciPy, and the command against the output.

## What to reuse

| Need | Where |
|---|---|
| syntax: varlist, if/in, weights, options | `lang/syntax.py` (`parse_standard`, `match_options`, `find_top`), `commands/_util.py` (`touse`, `eval_vector`, `plural`), `commands/summarize.py` (`weights`, `_level`, `g9`, `name12`) |
| groups of `by()` and strata | `commands/epitab.py` (`_strata`), `core/grouping.py` (`group_ids`) |
| `r()` and matrices | `s.r = {...}` (insertion order is the `return list` order, reversed where noted), `commands/matrix.py` (`Matrix`, `evaluate_matrix`) |
| estimation results and table | `estimation/results.py` (`Estimates`, `post`, `coef_table`, `CoefRow`) |
| maximum likelihood | `estimation/ml.py` (`maximize`, `print_log`): Stata's Newton steps and iteration log |
| predictions after a model | `estimation/postest.py` (`xb_all`, `esample`), `estimation/estat.py` (`_binary`) |
| prefixes (bootstrap-like) | `estimation/bootstrap.py` |
| graphs | `graphics/model.py` (`Graph`, `Layer`, `Panel`), `commands/graphcmd.py` (`finish`), `graphics/opts.py` (`graph_options`) |
| help pages | `src/opendta/help/*.sthlp`, `ALIASES` in `commands/helpcmd.py`; `tests/test_i18n.py` checks every command has one |
| errors | `core/errors.py` (`StataError(rc, message)`); an option not implemented yet raises `StataError(198, "option X is not available yet in OpenDTA")` instead of being ignored |

Tests use the `run` fixture of `tests/conftest.py` (`run("commands")` returns the output; `run.session.r`, `run.session.e`). The cloud workspace may not install pytest: then CI runs it.

## Blocks

Estimates from the work so far (±50%): session time and tokens generated; a session reads its context again at every step, so its total use is about 100 to 150 times the tokens generated. The order is a suggestion: cheap blocks that close an area first, then the ones most used in health research.

### 1. Epidemiology: finish ✅

Done in phase 12 (compat case 1601), except `rocreg` with `ml` and `cluster()`, which stay in *After 1.0*. The original plan:

`rocreg` (nonparametric and parametric ROC regression; bootstrap standard errors); `cs`, `cc` and `ir` with `estandard`, `istandard`, `standard()`, `binomial` and `pool`; `mhodds` with `by()`; `roctab` with `binomial`, `bamber` and `lorenz`; `roccomp` with `binormal`. Reuse `core/epi.py`, `commands/epitab.py`, `commands/roc.py`; remove the matching "not available yet" errors.

### 2. Data management: options ✅

Done in phase 13 (compat case 1701), except `import haver` (a proprietary format), `icd9p` (no copy of the CMS ICD-9-CM procedure list could be obtained), `import sasxport, describe` and `member()`, the value labels of `formats.xpf`, `export sasxport` with `rename` and `vallabfile()`, and `list` with `string()` and `display`, which stay in *After 1.0*. The original plan:

`import sas`, `import spss` and `import sasxport` (through `pyreadstat`, optional like for old `.dta`; `io/`), `file` with binary read/write (`commands/filecmd.py`), `infile` with `byvariable()` (`commands/infile.py`), `putexcel` to `.xls`, the missing options of `list` and `describe`. `import haver` needs a proprietary format: leave it out and say so. `icd9p` needs the public-domain ICD-9-CM procedure list (CMS); add it to `src/opendta/data/icd/` with `tools/make_icd_tables.py` if a copy can be obtained.

### 3. Programming ✅

Done in phase 14 (compat case 1801), except `makecns` with a constraint list and `displaycns` (they need `constraint`, from block 10), `matlist` with the equation options, the `window` dialogs and menus of the graphical interface, and loading plugins, which stay in *After 1.0*. The original plan:

`unab`, `tsunab`, `fvunab`, `unabcmd`, `fvexpand`, `fvrevar`, `tsrevar`, `mark`, `_rmcoll`, `matlist`, `tabdisp`, `sortpreserve` (program option), `signestimationsample`, `makecns`, `postfile`/`post`/`postclose`, `pause`, `sleep`, `window` (dialog commands, as far as the interface allows), `plugin` (cannot load Stata plugins: explain), `viewsource` (show OpenDTA's own source of a command, not Stata's). Mostly in `lang/` and `commands/program.py`, `commands/progcmd.py`, `commands/syntaxcmd.py`; factor variables in `estimation/fvars.py`.

### 4. Survival analysis (1.5–2.5 h, ~90k tokens)

`stsplit`, `stjoin`, `stfill`, `stgen`, `stbase`, `stvary`, `stci`, `strate`, `stptime`, `stir`, `stmh`, `stmc`, `ltable`, `sttocc`, `sttoct`, `snapspan`, `ct`, `stcurve` (graph after `stcox`/`streg`), `stcrreg` (Fine-Gray competing risks, ML); `streg, distribution(ggamma)` and `frailty()`; the `exactm` and `exactp` ties of `stcox`. Reuse `estimation/st.py`, `km.py`, `stcox.py`, `streg.py`, `stfit.py`.

### 5. Multilevel models (1–1.5 h, ~55k tokens)

`meglm` (the general form), `mecloglog`, `menbreg`, `meologit`, `meoprobit`, `mestreg`, `meqrlogit` and `meqrpoisson`; crossed random effects in `mixed`. Reuse `estimation/meglm.py` (adaptive quadrature) and `mixed.py`.

### 6. Panel data (2–3 h, ~100k tokens)

Random-effects models by quadrature: `xtprobit`, `xtcloglog`, `xtnbreg`, `xttobit`, `xtintreg`, `xtologit`, `xtoprobit`, `quadchk`; dynamic panels (GMM): `xtabond`, `xtdpd`, `xtdpdsys`; `xtivreg`, `xthtaylor`, `xtgls`, `xtpcse`, `xtregar`, `xtrc`, `xtfrontier`, `xtunitroot`; description: `xtline`, `xttab`, `xtdata`. Reuse `estimation/panel.py`, `panel2.py`, `xtgee.py`, the quadrature of `meglm.py`.

### 7. Base statistics (3–4 h, ~135k tokens)

Prefixes and resampling: `jackknife`, `permute`, `simulate` (reuse `estimation/bootstrap.py`); `ml` (the command around `estimation/ml.py`), `mlexp`; model building: `nestreg`, `stepwise`, `xi`, `mkspline`, `fp`, `mfp`; after estimation: `contrast`, `pwcompare`, `pwmean`, `margins` after more commands (`estimation/margins.py`), `suest`, `predictnl`, `testnl`, `linktest`; transformations: `boxcox`, `ladder`, `gladder`, `lnskew0`; tests: `bitest`, `ksmirnov`, `median`, `nptrend`, `runtest`, `tetrachoric`, `ktau`, `pcorr`; agreement and effect sizes: `icc`, `loneway`, `esize`, `ameans`; displays: `cumul`, `stem`, `dotplot`, `lv`, `lpoly`, `smooth`, `cusum`, `serrbar`, `sunflower`, `symplot`, `quantile`.

### 8. Multivariate statistics (2–3 h, ~90k tokens)

`pca`, `factor`, `rotate`, `alpha`, `cluster` (linkages, dendrogram, `cluster generate`), `discrim` (lda, qda, knn, logistic), `manova`, `mvreg`, `mvtest`, `canon`, `ca`, `mca`, `mds`, `biplot`, `procrustes`, `hotelling`, `screeplot`, `scoreplot`, each with its `estat` subcommands. New `core/multivariate.py`.

### 9. Time series (3–4 h, ~150k tokens)

`arch`, `arfima`, `mgarch`, `dfactor`, `sspace`, `ucm` (state space: share one Kalman filter with `arima` in `estimation/tsmodels.py`), `vec`, `vecrank`, `svar`, `fcast`, `varbasic`, `varlmar`, `varnorm`, `varstable`, `varwle`, `dfgls`, `pperron`, `prais`, `rolling`, `tsappend`, `tsfill`, `tsreport`, `tsfilter`, `tssmooth`, `pergram`, `cumsp`, `xcorr`, `wntestb`, `wntestq`, `psdensity`; `arima` with `sarima()`, `condition` and `diffuse`; `irf` after `svar` and `vec`, with bootstrap standard errors, and `irf ctable`, `cgraph`, `ograph`, `add` (`estimation/irf.py`).

### 10. Models (6–8 h, ~285k tokens; can be split)

Binary and count: `cloglog`, `scobit`, `hetprobit`, `glogit`, `binreg`, `fracreg`, `gnbreg`, `zip`, `zinb`, `tpoisson`, `tnbreg`; censored and selection: `intreg`, `truncreg`, `heckman`, `heckprobit`, `heckoprobit`, `biprobit`; instrumental variables: `ivprobit`, `ivtobit`, `ivpoisson`, `ivregress liml` and `gmm`; choice models: `clogit`, `mprobit`, `rologit`, `slogit`, `asclogit`, `asmprobit`, `asroprobit`; exact: `exlogistic`, `expoisson`; systems and others: `cnsreg` and `constraint`, `nl`, `nlsur`, `gmm`, `rreg`, `sureg`, `reg3`, `eivreg`, `vwls`, `frontier`; treatment effects: `teffects`, `stteffects`, `etregress`, `etpoisson`, `eteffects`. Reuse `estimation/models.py`, `models2.py`, `glm.py`, `ml.py`. Suggested split: 10a binary/count/censored, 10b IV/choice/exact, 10c systems and treatment effects.

### 11. Other suites (3–4 h, ~130k tokens)

`gsem`, `bayesmh` and the `bayes` prefix, `irt`; `sem` with `mlmv`, `adf`, `group()` and the other `estat` subcommands (`estimation/sem.py`); the other `power` methods, the `lrchi2` and `fisher` tests and graphs (`commands/power.py`); `mi impute mvn`, `monotone`, `truncreg`, `intreg` and `nbreg` (`estimation/mi.py`); `svy:` with more estimators (`estimation/svy.py`).

## Before 1.0

When every block is done: the Stata run of all compat cases without reference, the differences fixed, the version set to 1.0.0 in `src/opendta/__init__.py` and `pyproject.toml`, the release notes written, and the tag `v1.0.0` pushed (see [development.md](development.md#publishing-a-version)).

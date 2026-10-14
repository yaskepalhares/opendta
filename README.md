# OpenDTA

**English** | [Português](README.pt-BR.md)

**An open-source interpreter for statistical analysis built around do-files and `.dta` data files.**

OpenDTA runs analysis scripts (do-files), reads and writes `.dta` datasets and comes with a desktop interface: Results, Command, Review, Variables and Properties windows, plus a data browser. The goal is a free and transparent tool for reproducible statistical work, from data management to estimation, complex survey analysis and graphics, able to run existing scripts without changes.

It is a personal study project under active development, released under the [MIT License](LICENSE). The full plan is in the [ROADMAP](ROADMAP.md).

## Download

[![Download for macOS (Apple Silicon)](https://img.shields.io/badge/Download-macOS%20Apple%20Silicon-1f6feb?style=for-the-badge&logo=apple&logoColor=white)](https://github.com/yaskepalhares/opendta/releases/latest/download/OpenDTA-macos-apple-silicon.dmg)
[![Download for macOS (Intel)](https://img.shields.io/badge/Download-macOS%20Intel-57606a?style=for-the-badge&logo=apple&logoColor=white)](https://github.com/yaskepalhares/opendta/releases/latest/download/OpenDTA-macos-intel.dmg)
[![Download for Windows (64-bit)](https://img.shields.io/badge/Download-Windows%2064--bit-0078d4?style=for-the-badge&logo=windows&logoColor=white)](https://github.com/yaskepalhares/opendta/releases/latest/download/OpenDTA-windows-x64-setup.exe)

**macOS:** requires macOS 12 or later. Choose **Apple Silicon** for Macs with an M1, M2, M3 or M4 chip and **Intel** for the others (check  → About This Mac). Open the `.dmg` and drag OpenDTA to Applications. The app is not notarized by Apple, so on first launch open **System Settings → Privacy & Security** and click **Open Anyway** (only once).

**Windows:** requires 64-bit Windows 10 or 11 on an Intel or AMD processor. Run `OpenDTA-windows-x64-setup.exe`; it installs for the current user only, with no administrator rights, and adds OpenDTA to the **Open with** menu of `.do`, `.ado` and `.dta` files. The installer is not signed, so Windows SmartScreen may show "Windows protected your PC": click **More info → Run anyway**. If you prefer not to install anything, download the [portable version](https://github.com/yaskepalhares/opendta/releases/latest/download/OpenDTA-windows-x64-portable.zip), unzip it and run `OpenDTA.exe`.

All versions are on the [Releases](https://github.com/yaskepalhares/opendta/releases) page.

## Current state: phases 5 to 9 (estimation, complex surveys, graphics, interface and language)

What already works:

- do-file reading: `*`, `//` and `/* */` comments (nested), `///` continuation and `#delimit ;`;
- local and global macros, including nesting (``` ``i'' ```), `` `=exp' ``, `` `++i' ``, `local ++i` and extended functions (`word count`, `word # of`, `display`, `list`, `piece`, `env`);
- expressions with Stata's operator precedence and missing-value rules (`.`, `.a`–`.z`);
- 123 function names (counting synonyms such as `trim`/`strtrim`): mathematical, string, logical, distribution, date and random-number functions;
- `display` with formats (`%9.2f`, `%10.3e`, `%-12s`, `%td`…) and directives (`_col`, `_skip`, `_dup`, `_n`, `_continue`, `as`…);
- `forvalues`, `foreach`, `while`, `if`/`else if`/`else`, `continue` and `continue, break`;
- the `quietly`, `noisily` and `capture` prefixes (also on `{ }` blocks);
- `do`, `run`, `include`, `args`, `tokenize`, `confirm`, `scalar`, `macro`, `set`, `version`, `cd`, `pwd`, `exit` and `error`;
- errors with Stata's messages and `r(#)` return codes, followed by an OpenDTA explanation of what went wrong in the command and how to fix it (turn it off with `set hints off`, or with Help → Explain Errors);
- `set numerics precise`: the same estimates with more exact arithmetic (QR in `regress`, `D.`/`S.` differences in double, a final Newton step in `ml` and `optimize()`, numerical derivatives with Richardson extrapolation, new variables in double); the default `set numerics stata` reproduces Stata's rounding;
- scientific notation with superscript exponents on screen (`1.0000×10¹⁰` instead of `1.00000e+10`); `set superscript off` restores Stata's form, and text that can be read back (`string()`, macros, `file write`, exports) always uses Stata's form;
- graphical interface with the Results, Command, Review, Variables and Properties windows, menus, a toolbar with OpenDTA's own icons (light/dark) and Preferences (app icon and font);
- text console and batch mode (`opendta -b do file.do`, which writes `file.log` like `stata -b`);
- programming: `program`, `syntax`, `gettoken`, `marksample`, `return`/`ereturn`, `tempvar`/`tempfile`, `preserve`/`restore`, extended macro functions, `matrix`, third-party `.ado` files on the adopath, `log`/`cmdlog`, `assert`, `set trace` and `timer`;
- Viewer with SMCL `help` (OpenDTA's own pages and those of installed packages);
- data in memory with Stata's storage types (including `float` precision), varlists, `if`/`in`, `by`/`bysort`, `_n`/`_N` and subscripts;
- `set obs`, `input`, `generate`, `replace`, `drop`, `keep`, `list`, `describe`, `count`, `sort`, `gsort`, `rename`, `order`, `label`, `format`, `compress` and `recast`;
- `.dta` files: OpenDTA's own reader for formats 113 to 119 (Stata 8 onward) and writer for formats 114, 115, 117, 118 and 119, with labels, notes and characteristics; older formats are read through ReadStat when the `pyreadstat` package is installed;
- `use` (with varlist, `if` and `in`), `save`, `saveold`, `notes`, `char`, `type`, `erase`;
- delimited text: `import delimited`, `export delimited`, `insheet` and `outsheet`;
- spreadsheets: `import excel` and `export excel`;
- free-format and fixed-column text: `infile` (including `.dct` dictionaries) and `infix`; `file open/write/read/close`;
- File > Open/Save/Save as, dropping a `.dta` on the window, Data Editor (`browse` and `edit`, with every edit recorded as a command) and Variables Manager (`varmanage`);
- descriptive statistics: `summarize`, `tabstat`, `tabulate` (with `chi2`, `exact`, `gamma`, `taub`…), `tab1`, `tab2`, `table`, `ci`, `ttest`, `prtest`, `correlate`, `pwcorr`, `centile`, `pctile`, `xtile`, with `fweight`, `aweight` and `iweight`;
- data management: `egen` (group and row functions), `collapse`, `contract`, `expand`, `fillin`, `merge`, `append`, `joinby`, `cross`, `reshape`, `duplicates`, `isid`, `levelsof`, `recode`, `encode`/`decode`, `destring`/`tostring`, `split`, `mvencode`/`mvdecode`;
- random numbers: 64-bit Mersenne Twister, `set seed` and `runiform()`, `rnormal()`, `rbinomial()`, `rpoisson()` and other distributions;
- Mata: `mata ... end` blocks, user functions, structs, classes (inheritance, virtual methods) and pointers, `optimize()` and `moptimize()`, about 150 library functions and the link to the data (`st_data`, `st_store`, `st_local`, `st_matrix`…);
- estimation: `regress`, `logit`, `logistic`, `probit`, `ologit`, `oprobit`, `mlogit`, `poisson`, `nbreg`, `tobit`, `glm`, `xtreg`, `xtlogit`, `xtpoisson`, `areg`, `ivregress`, `qreg`, `anova` and `oneway`, with factor variables, time-series operators, weights and robust, cluster and bootstrap standard errors;
- postestimation: `predict`, `test`, `testparm`, `lincom`, `nlcom`, `margins`, `estimates`, `lrtest`, `hausman`, `estat` (`ic`, `vce`, `hettest`, `ovtest`, `classification`, `gof`, `bootstrap`…);
- complex surveys: `svyset` (several stages, replicate weights), `svydescribe`, `estat effects` and the `svy:` prefix with `mean`, `proportion`, `total`, `ratio`, `tabulate`, `regress`, `logit`, `logistic`, `probit`, `ologit`, `oprobit`, `mlogit`, `poisson`, `nbreg`, `glm` and `tobit` (Taylor linearization, BRR, jackknife, bootstrap and SDR);
- multilevel models: `mixed` (ML and REML, nested levels, random slopes and covariance structures) and `melogit`, `meprobit` and `mepoisson` (adaptive quadrature), with `estat icc`;
- survival analysis: `stset` (single and multiple records, delayed entry, origin and scale), `stdescribe`, `stsum`, `sts list`/`generate`/`test`/`graph` (Kaplan-Meier, Nelson-Aalen, log-rank and Wilcoxon tests), `stcox` (Breslow and Efron ties, strata, robust variance), `estat phtest` and `streg` (exponential, Weibull, Gompertz, lognormal and loglogistic);
- tests: `ranksum`, `signrank`, `kwallis`, `spearman`, `swilk`, `sktest`, `sdtest` and the immediate commands `cii`, `ttesti`, `prtesti` and `tabi`;
- describing data: `codebook` and `misstable`; results to Excel with `putexcel`;
- graphics with OpenDTA's own look: `twoway` (`scatter`, `line`, `lfit`, `qfit`, `function`…), `histogram`, `kdensity`, `graph bar`, `graph box`, `graph pie`, `graph matrix`, `graph combine`, `by()`, `graph export` (PNG, PDF, SVG…), `lowess`, `qnorm`, `rvfplot`, `marginsplot` and a Graph window with a simple Graph Editor;
- Do-file Editor with tabs, syntax highlighting and execution of the selection (Ctrl+D) or of the whole file;
- dialog boxes in the Statistics and Graphics menus, which show the command before running it;
- interface in English or Portuguese (**Edit → Preferences** or `set language pt`); output and Stata-style error messages stay in English unless translation is requested (`set messages translated`).

The outputs of cases `0603` (survey replicates), `0901` and `0902` (survival analysis) still await checking against the Stata log, as do the cases added in block 2 (`0408`, `0409`, `1001`, `1002`); all other cases in `compat/do` have been checked. The new estimators were checked against statsmodels (Cox, Kaplan-Meier, log-rank, Wilcoxon, `mixed`), against direct maximization of the likelihoods (`streg`) and against exact numerical integration (`melogit`, `mepoisson`).

## Installing from source

Requires Python 3.11 or later. In the Terminal:

```bash
git clone https://github.com/yaskepalhares/opendta.git
cd opendta
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

The commands are the same on Linux. On Windows, activate the environment with `.venv\Scripts\activate`.

### macOS application

To have OpenDTA in the Applications folder, with a Dock icon and "Open With" for `.do` and `.dta` files, build the package in the Terminal (inside the `opendta` folder):

```bash
tools/build_macos.sh
```

The script creates `dist/OpenDTA.app` and `dist/OpenDTA-<version>.dmg`. Open the `.dmg` and drag OpenDTA to Applications. The signature is local (ad hoc): on first launch macOS may ask for confirmation in **System Settings → Privacy & Security → Open Anyway**. On other systems, `pyinstaller packaging/opendta.spec` builds the `dist/OpenDTA` folder with the executable.

## Usage

| Command | What it does |
|---|---|
| `opendta` | opens the graphical interface |
| `opendta file.do` | opens the interface and runs the do-file |
| `opendta --console` | text console |
| `opendta -b do file.do` | batch mode: writes `file.log` in the current directory |

## Icon

The application icon ("Table") has a light (default) and a dark version, chosen in **Edit → Preferences**. The full sets (Windows `.ico`, macOS `.icns`, Linux, web) are in `assets/icons/opendta/` and `assets/icons/opendta-dark/`, and the SVG sources are in `assets/icons/icon-work/`.

## Tests

The automated tests run with `pytest`. The comparison with the Stata 14 reference logs runs with `python tools/compare.py`. The complete step-by-step guide for macOS is in [docs/testing.md](docs/testing.md).

## Code layout

```
src/opendta/
  core/       missing values, display formats, r(#) errors, output
  lang/       do-file lexer, macros, expressions, functions, interpreter
  commands/   commands (one module per group)
  io/         reading and writing .dta, text and spreadsheets
  estimation/ models, postestimation, svy and bootstrap
  graphics/   graphics (matplotlib, OpenDTA's own schemes)
  mata/       the Mata language
  gui/        graphical interface (PySide6), dialogs and Graph Editor
  help/       help pages (.sthlp)
  i18n.py     interface language (English or Portuguese)
  session.py  state of a session (macros, scalars, r(), e(), c())
  cli.py      command line: interface, console and batch
compat/
  do/         verification do-files
  expected/   reference logs produced by Stata 14
packaging/    PyInstaller spec (packaged application)
tools/
  compare.py         runs compat/do in OpenDTA and compares with compat/expected
  build_macos.sh     builds OpenDTA.app and the .dmg
  build_windows.ps1  builds the Windows installer and the portable .zip
```

## License

OpenDTA is released under the [MIT License](LICENSE): it may be used, copied, modified and redistributed freely, as long as the copyright notice and the license accompany the copies. The packaged application includes third-party components under their own licenses (Qt/PySide6 under LGPL-3.0, NumPy, SciPy, Matplotlib and others), listed in [packaging/THIRD_PARTY_NOTICES.md](packaging/THIRD_PARTY_NOTICES.md).

## Publishing a new version

The workflow [`.github/workflows/release.yml`](.github/workflows/release.yml) builds the `.dmg` files for Apple Silicon and Intel, the Windows installer and the portable `.zip`, and creates the Release when a `v*` tag is pushed:

```bash
git tag v0.2.1
git push origin v0.2.1
```

The download buttons always point to the latest Release. Remember to update the version in `src/opendta/__init__.py` and `pyproject.toml` before tagging. To build the Windows files locally, run `powershell -ExecutionPolicy Bypass -File tools\build_windows.ps1` (needs 64-bit Python 3.11+ and [Inno Setup 6](https://jrsoftware.org/isinfo.php)).

## Notice

Stata is a registered trademark of StataCorp LLC. OpenDTA has no relationship with StataCorp and uses none of its code; compatibility with Stata 14 do-files is built only from public documentation and observation of behavior. See the [development policy](docs/clean-room-policy.md).

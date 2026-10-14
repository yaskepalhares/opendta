# How to test and compare with Stata 14 (macOS)

**English** | [Português](testing.pt-BR.md)

These instructions are for macOS, with the Terminal (zsh) and Stata 14 for Mac. The repository is in `~/opendta`. If it is in another folder, change the path in the commands.

> **Careful when pasting commands into the Terminal.** In interactive zsh, `#` does not start a comment: a line such as `pytest  # run everything` fails. The blocks below have no comments and can be pasted as they are.

There are two levels of testing: the automated tests and the comparison with Stata.

## 0. Set up the environment (once)

1. Install Python 3.11 or later, with the installer from [python.org](https://www.python.org/downloads/macos/) or with Homebrew (`brew install python@3.12`).
2. In the Terminal, clone the repository and create the virtual environment:
   ```bash
   cd ~
   git clone https://github.com/yaskepalhares/opendta.git
   cd opendta
   python3 -m venv .venv
   source .venv/bin/activate
   pip install -e ".[dev]"
   ```

The `[dev]` extra installs pytest, ReadStat (`pyreadstat`) and pandas, which the tests use.

## Before each round of tests

Open the Terminal, go to the folder and activate the environment:

```bash
cd ~/opendta
source .venv/bin/activate
```

The prompt now shows `(.venv)`. To test a pull request before it enters `main`, switch to its branch (the name is shown on the PR):

```bash
git fetch
git checkout branch-name
git pull
pip install -e ".[dev]"
```

To go back to the main version:

```bash
git checkout main
git pull
```

Repeating `pip install -e ".[dev]"` is only needed when the dependencies change, but running it every time does no harm.

## 1. Automated tests (`pytest`)

They check the internal behavior and run on GitHub after every change (Linux, Windows and macOS). To run them on your Mac:

```bash
pytest
```

The expected result is `N passed`, with no `failed`. If something fails, copy the Terminal output and send it.

## 2. Open the program

```bash
opendta
```

Other ways to use it:

```bash
opendta --console
opendta -b do file.do
```

The first opens the text console (prompt `. `). The second runs a do-file in batch mode and writes `file.log` in the current folder, like `stata -b`.

## 3. Comparison with Stata 14

The `compat/` folder holds verification do-files (`compat/do/`) and the reference output produced by Stata 14 (`compat/expected/`). This comparison is what decides whether OpenDTA is faithful.

### Produce the references in Stata (for each new batch of do-files)

1. Open Stata 14.
2. In the Command window, go to the repository's `compat` folder and run the generator:
   ```stata
   cd "~/opendta/compat"
   do gerar_esperados.do
   ```
3. The `compat/expected/*.log` files appear. The generator forces `set dp period` while it runs and restores your setting at the end.
4. Send the logs: attach the new files to the conversation, or commit them from the Terminal:
   ```bash
   cd ~/opendta
   git add compat/expected
   git commit -m "Stata 14 reference logs"
   git push
   ```

### Compare

With the environment activated (`source .venv/bin/activate`), compare all do-files:

```bash
python tools/compare.py
```

To show the lines that differ:

```bash
python tools/compare.py --show-diff
```

To compare a single file, by the beginning of its name:

```bash
python tools/compare.py 0104
```

For each do-file, the output says whether the results are equal (`ok`), different (`XX`), different in a known and documented way (`~~`, see `compat/known_differences.txt`) or still without a reference (`?`). Log headers, blank lines, trailing spaces and the date stamp of `describe` are ignored.

### `.dta` files

The `.dta` reader and writer are OpenDTA's own. In the automated tests, every file written is read back by ReadStat, an independent implementation of the specification. Stata remains the final judge.

1. **Stata → OpenDTA.** In Stata, save a file with labels, notes and extended missing values:
   ```stata
   sysuse auto, clear
   notes: test note
   replace rep78 = .a in 1
   save "~/Desktop/auto14.dta", replace
   saveold "~/Desktop/auto13.dta", replace
   ```
   In OpenDTA, open the file and check it:
   ```stata
   use "~/Desktop/auto14.dta", clear
   describe
   notes
   list in 1/5
   ```
   Repeat with `auto13.dta`. The file can also be opened with File > Open or by dragging it onto the window.
2. **OpenDTA → Stata.** In OpenDTA:
   ```stata
   save "~/Desktop/from_opendta.dta", replace
   saveold "~/Desktop/from_opendta13.dta", replace
   ```
   In Stata, open each one and run `describe`, `notes`, `char list` and `list`. Any warning from Stata when opening the file is a difference.

The do-files `compat/do/0104_arquivos.do` to `0107_infile_infix.do` create and delete their own files (`.dta`, CSV, Excel and text) in the current folder, with names starting with `odta_`.

## Reporting a difference found in use

Open an *issue* with the **Difference from Stata** template: a minimal do-file, the Stata 14 output and the OpenDTA output. Each confirmed difference becomes a case in `compat/do/`, so that it does not happen again.

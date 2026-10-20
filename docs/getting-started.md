# Getting started

**English** | [Português](getting-started.pt-BR.md)

## Install

Download the installer for your system from the [latest release](https://github.com/yaskepalhares/opendta/releases/latest):

| System | File |
|---|---|
| macOS with Apple Silicon (M1, M2, M3, M4…) | `OpenDTA-macos-apple-silicon.dmg` |
| macOS with an Intel processor | `OpenDTA-macos-intel.dmg` |
| Windows 10 or 11, 64-bit | `OpenDTA-windows-x64-setup.exe` (installer) or `OpenDTA-windows-x64-portable.zip` |

**macOS** (12 or later): open the `.dmg` and drag OpenDTA to Applications. Not sure which file? Open  → About This Mac: "Chip Apple M…" is Apple Silicon. The app is not notarized by Apple, so on first launch open **System Settings → Privacy & Security** and click **Open Anyway** (only once).

**Windows**: the installer works for the current user only, with no administrator rights, and adds OpenDTA to the **Open with** menu of `.do`, `.ado` and `.dta` files. It is not signed, so SmartScreen may show "Windows protected your PC": click **More info → Run anyway**. The portable `.zip` runs from any folder.

**Updates**: OpenDTA looks for a new release when it opens (at most once a day; turn it off in **Settings**, Ctrl+, / ⌘,) and in **Help → Check for Updates**. **Install and Restart** downloads the right file, installs it and opens OpenDTA again.

### From source

Requires Python 3.11 or later:

```bash
git clone https://github.com/yaskepalhares/opendta.git
cd opendta
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

On Windows, activate the environment with `.venv\Scripts\activate`. Optional extras: `pip install -e ".[legacy]"` reads `.dta` files from Stata 7 and earlier (ReadStat), and `pip install -e ".[xls]"` reads and writes Excel 97-2003 `.xls` files.

## Run

| Command | What it does |
|---|---|
| `opendta` | opens the window |
| `opendta file.do` | opens the window and runs the do-file |
| `opendta --console` | text console (prompt `. `) |
| `opendta -b do file.do` | batch mode: writes `file.log` in the current folder, like `stata -b` |

## First steps

The window has the familiar parts: Results, Command, Review, Variables and Properties, plus the Data Editor (`browse`, `edit`), the Do-file Editor (`doedit`), the Viewer (`help`) and the Graph window. Type commands in the Command window, or open a do-file and run it with Ctrl+D (⌘+D).

```stata
use mydata, clear
describe
summarize age income
regress income age i.sex
```

When a command fails, OpenDTA shows Stata's message and code and, below them, an explanation for beginners: what went wrong with the data at hand, why, how to fix it, and a link to the help page. Turn it off with `set hints off`.

In the Command window, OpenDTA asks before a command that discards much of the data (`clear` with unsaved changes, `drop` or `keep` removing many observations, `save, replace` over a much larger file); each one has its own switch in **Settings → General → Safety**. The interface can be in English or Portuguese (**Settings**, or `set language pt`). Results and Stata-style messages stay in English unless you ask for the translation (`set messages translated`).

Type `help` for the list of commands, or see [what OpenDTA does](features.md).

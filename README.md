# OpenDTA

**English** | [Português](README.pt-BR.md)

**Do-files and `.dta` data, in a free and open program.**

OpenDTA runs do-files written for Stata 14, reads and writes `.dta` datasets and has the desktop interface people already know: the Results, Command, Review, Variables and Properties windows, the Data Editor, the Do-file Editor, the Viewer and the Graph window. It covers data management, descriptive statistics, estimation, complex surveys, survival analysis, multilevel models, time series, multiple imputation, structural equation models, Mata and graphics. When a command fails, it explains the error for beginners, looking at the data in memory, in English or Portuguese.

![OpenDTA with a regression and the explanation of an error](docs/screenshot.png)

- [Getting started](docs/getting-started.md)
- [What OpenDTA does](docs/features.md)
- [Roadmap and the commands still missing](ROADMAP.md)
- [Testing and comparing with Stata](docs/testing.md)
- [Development](docs/development.md)
- [Development policy (clean room)](docs/clean-room-policy.md)

## Install

<table>
<tr><th>macOS 12 or later</th><th>Windows 10 or 11, 64-bit</th><th>Linux and other systems</th></tr>
<tr>
<td align="center"><a href="https://github.com/yaskepalhares/opendta/releases/latest/download/OpenDTA-macos-apple-silicon.dmg"><img src="assets/download/macos-apple-silicon.svg" alt="Download for macOS, Apple Silicon (.dmg)" width="264"></a><br><a href="https://github.com/yaskepalhares/opendta/releases/latest/download/OpenDTA-macos-intel.dmg"><img src="assets/download/macos-intel.svg" alt="Download for macOS, Intel (.dmg)" width="264"></a></td>
<td align="center"><a href="https://github.com/yaskepalhares/opendta/releases/latest/download/OpenDTA-windows-x64-setup.exe"><img src="assets/download/windows-installer.svg" alt="Download the Windows installer (.exe)" width="264"></a><br><a href="https://github.com/yaskepalhares/opendta/releases/latest/download/OpenDTA-windows-x64-portable.zip"><img src="assets/download/windows-portable.svg" alt="Download the portable Windows version (.zip)" width="264"></a></td>
<td align="center" valign="top"><a href="docs/getting-started.md#from-source"><img src="assets/download/source.svg" alt="Install from source on Linux and other systems" width="264"></a></td>
</tr>
</table>

Not sure which Mac you have? Apple menu → **About This Mac**: "Apple M…" is Apple Silicon, "Intel" is Intel. The apps are not signed, so the first launch needs one confirmation; see [Getting started](docs/getting-started.md#install). Installed copies check for updates by themselves. All versions are on the [Releases](https://github.com/yaskepalhares/opendta/releases) page.

## Progress

[![Commands implemented per area](https://yaskepalhares.github.io/opendta/roadmap-coverage.svg)](https://yaskepalhares.github.io/opendta/)

Each tile is an area of commands, sized by how many it has and coloured by the share already implemented; it is redrawn on every change. **[Open the interactive map](https://yaskepalhares.github.io/opendta/)** to see every command and search for one, or read the [roadmap](ROADMAP.md).

## Contributing

Read [CONTRIBUTING.md](CONTRIBUTING.md) before opening an issue or a pull request. The most useful report is a difference from Stata: open an issue with the [Difference from Stata](https://github.com/yaskepalhares/opendta/issues/new?template=difference-from-stata.md) template and a short do-file. The interface translations live in `src/opendta/i18n.py`, and packaging is described in [docs/development.md](docs/development.md#publishing-a-version).

## Acknowledgements

Built with Python, NumPy, SciPy, Matplotlib, Qt (PySide6) and openpyxl; old `.dta` files are read through ReadStat. The third-party components and their licenses are listed in [packaging/THIRD_PARTY_NOTICES.md](packaging/THIRD_PARTY_NOTICES.md).

OpenDTA is a personal study project, independent and not affiliated with StataCorp. Stata is a registered trademark of StataCorp LLC; OpenDTA uses none of its code, and compatibility is built only from public documentation and observed behavior. Released under the [MIT License](LICENSE).

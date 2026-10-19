# Development

**English** | [Português](development.pt-BR.md)

How the code is organized, how it is tested and how a version is published. To propose a change, see [CONTRIBUTING.md](../CONTRIBUTING.md).

## Code layout

```
src/opendta/
  core/       missing values, display formats, r(#) errors, output, error explanations (hints/)
  lang/       do-file lexer, macros, expressions, functions, interpreter
  commands/   commands (one module per group)
  io/         reading and writing .dta, text and spreadsheets
  estimation/ models, postestimation, svy, sem, irf and bootstrap
  graphics/   graphics (matplotlib, OpenDTA's own schemes)
  mata/       the Mata language
  gui/        graphical interface (PySide6), dialogs and Graph Editor
  help/       help pages (.sthlp)
  i18n.py     interface language (English or Portuguese)
  session.py  state of a session (macros, scalars, r(), e(), c())
  cli.py      command line: interface, console and batch
compat/
  do/         verification do-files
  expected/   reference logs produced by Stata
packaging/    PyInstaller spec, Windows installer script and third-party notices
tools/
  compare.py         runs compat/do in OpenDTA and compares with compat/expected
  build_macos.sh     builds OpenDTA.app and the .dmg
  build_windows.ps1  builds the Windows installer and the portable .zip
```

The plan for the commands still missing, block by block, is in [next-steps.md](next-steps.md).

## Tests

`pip install -e ".[dev]"` installs pytest, ReadStat (`pyreadstat`) and pandas. Then:

```bash
pytest
python tools/compare.py
```

`pytest` checks the internal behavior and runs on GitHub after every change (Linux, Windows and macOS). `tools/compare.py` runs every do-file in `compat/do/` and compares its log with the Stata reference in `compat/expected/`; known differences are listed in `compat/known_differences.txt`. The step-by-step guide, including how to produce new references in Stata, is in [testing.md](testing.md).

## Icon

The application icon ("Table") has a light (default) and a dark version, chosen in **Settings**. The full sets (Windows `.ico`, macOS `.icns`, Linux, web) are in `assets/icons/opendta/` and `assets/icons/opendta-dark/`, and the SVG sources in `assets/icons/icon-work/`.

## Publishing a version

1. Update the version in `src/opendta/__init__.py` and `pyproject.toml` (they must match the tag).
2. Merge the pull request into `main`.
3. Push a tag:
   ```bash
   git tag v1.0.0
   git push origin v1.0.0
   ```

The workflow [`.github/workflows/release.yml`](../.github/workflows/release.yml) checks that the tag matches the version, builds the `.dmg` files for Apple Silicon and Intel, the Windows installer and the portable `.zip`, and creates the Release with the notes in `.github/release-notes.md`. The download links always point to the latest Release, and installed copies offer the update.

To build locally: `tools/build_macos.sh` on a Mac, or `powershell -ExecutionPolicy Bypass -File tools\build_windows.ps1` on Windows (needs 64-bit Python 3.11+ and [Inno Setup 6](https://jrsoftware.org/isinfo.php)). On other systems, `pyinstaller packaging/opendta.spec` builds the `dist/OpenDTA` folder.

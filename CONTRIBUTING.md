# Contributing

**English** | [Português](CONTRIBUTING.pt-BR.md)

Thank you for helping. OpenDTA is a personal study project, but reports and pull requests are welcome. Issues may be written in English or Portuguese.

## Reporting a problem

- **A result differs from Stata 14**: open an issue with the [Difference from Stata](https://github.com/yaskepalhares/opendta/issues/new?template=difference-from-stata.md) template. A short do-file that reproduces it, with the output of both programs, is the most useful report there is.
- **Something breaks or the window misbehaves**: use the [Bug report](https://github.com/yaskepalhares/opendta/issues/new?template=bug-report.md) template.
- **A Stata 14 command you need is missing**: check the list in the [ROADMAP](ROADMAP.md#after-10) and use the [Command request](https://github.com/yaskepalhares/opendta/issues/new?template=command-request.md) template. Requests help decide the order.

## Rules of the project

1. **Clean room.** OpenDTA is written from the Stata manuals, the public help, the `.dta` specification and the observed behavior of Stata 14. Never copy or adapt StataCorp code or files (`.ado`, `.mo`, `.sthlp`, `.dlg`, `.scheme`), and never decompile Stata. See the [development policy](docs/clean-room-policy.md).
2. **Language.** Code, comments, docstrings, names and commit messages are in English. Each documentation page (README, ROADMAP, CONTRIBUTING, `docs/`) has a Portuguese copy with the `.pt-BR.md` suffix and an **English | Português** selector at the top; when a page changes, update both. The interface is translated in `src/opendta/i18n.py`.
3. **Output like Stata.** Messages and results imitate Stata 14 (in English). After an error message, OpenDTA adds its own explanation (`src/opendta/core/hints/`).
4. **Undocumented details** get a `VERIFY` comment in the code and a case in `compat/do/`, so that a Stata run can settle them.
5. **Tests.** Before opening a pull request, `pytest` must pass and `python tools/compare.py` must show no new differences.
6. **Pull requests.** Work on a branch and open a pull request against `main`; nothing is pushed directly to `main`.

## Setting up

```bash
git clone https://github.com/yaskepalhares/opendta.git
cd opendta
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
pytest
```

The code layout, the tests and the release process are described in [docs/development.md](docs/development.md); the comparison with Stata, in [docs/testing.md](docs/testing.md).

# Instructions for development sessions

- Project: OpenDTA, a do-file interpreter with a Stata 14-style interface. The goal is in ROADMAP.md.
- Follow docs/clean-room-policy.md: no StataCorp code, files or graphical elements, and no decompiling. Use only manuals, public help, the .dta specification and observation of behavior.
- Language: everything in the repository is in English: documentation, code comments, docstrings, code names and commit messages. Each documentation page (README, ROADMAP, docs/) has a Portuguese copy with the `.pt-BR.md` suffix and an **English | Português** selector at the top; when a page changes, update both. The interface is bilingual through src/opendta/i18n.py.
- Error messages and output imitate Stata (in English); after an error message, OpenDTA adds its own explanation (core/hints.py, English or Portuguese following the interface language), turned off with `set hints off` in compare.py and in the format tests.
- Undocumented details get a `VERIFY` mark in the code and a case in compat/do/.
- Before opening a PR: `pytest` passing and `python tools/compare.py` with no regressions.
- Work on a branch and open a PR for review. Do not push directly to main.

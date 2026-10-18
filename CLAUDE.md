# Instructions for development sessions

## Project rules

- Project: OpenDTA, a do-file interpreter with a Stata-style interface. The target version is in the goal of ROADMAP.md; elsewhere, don't name the version.
- Follow docs/clean-room-policy.md: no StataCorp code, files or graphical elements, and no decompiling. Use only manuals, public help, the .dta specification and observation of behavior.
- Language: everything in the repository is in English: documentation, code comments, docstrings, code names and commit messages. Each documentation page (README, ROADMAP, docs/) has a Portuguese copy with the `.pt-BR.md` suffix and an **English | Português** selector at the top; when a page changes, update both. The interface is bilingual through src/opendta/i18n.py.
- Error messages and output imitate Stata (in English); after an error message, OpenDTA adds its own explanation (core/hints.py, English or Portuguese following the interface language), turned off with `set hints off` in compare.py and in the format tests.
- Undocumented details get a `VERIFY` mark in the code and a case in compat/do/.
- Before opening a PR: `pytest` passing and `python tools/compare.py` with no regressions.
- Work on a branch and open a PR for review. Do not push directly to main.

## Working guidelines

These bias toward caution over speed. For trivial tasks, use judgment.

### 1. Think before coding

- State assumptions explicitly. If uncertain, ask.
- If a request has more than one reading, present them; don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If the manual does not settle how Stata behaves, say so: don't guess silently. Mark it `VERIFY` and add a compat/do/ case instead.

### 2. Simplicity first

- Minimum code that solves the problem. No features, options or configurability beyond what was asked.
- No abstractions for single-use code, no error handling for impossible scenarios.
- If 200 lines could be 50, rewrite it.

### 3. Surgical changes

- Touch only what the request needs. Don't "improve" adjacent code, comments or formatting, and don't refactor what isn't broken.
- Match the existing style of the module.
- Remove imports, variables and functions that your own change made unused; mention unrelated dead code instead of deleting it.
- Every changed line should trace back to the request.

### 4. Goal-driven execution

Turn each task into something verifiable before writing code:

- "Fix a difference from Stata" → reproduce it with a do-file in compat/do/ (or a pytest case), then make the output match.
- "Add a command or option" → write the tests from the manual's examples, then make them pass.
- "Refactor X" → `pytest` and `tools/compare.py` give the same result before and after.

For multi-step tasks, state a brief plan:

```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
```

A known difference goes in compat/known_differences.txt only with the reason it cannot match; it is not a way to make a failing case pass.

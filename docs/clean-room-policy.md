# Development policy ("clean room")

**English** | [Português](clean-room-policy.pt-BR.md)

OpenDTA reimplements Stata's behavior from public sources. It uses no StataCorp code.

## May be used as a source

- The Stata manuals (PDFs), the help files (`help`) and the public FAQs
- The public specification of the `.dta` format
- Papers, books and teaching materials about Stata and about the statistical methods
- Observation of behavior: running a do-file in Stata 14 and comparing its output with OpenDTA's
- Open libraries with a compatible license (for example, ReadStat for `.dta`)

## May not be used

- `.ado`, `.mo`, `.mlib`, `.sthlp`, `.dlg` or `.scheme` files distributed with Stata, neither to copy excerpts nor as a basis for adaptation
- Decompiling, disassembling or inspecting Stata's executables and libraries
- Stata's icons, logos and other graphical elements

## In practice

- Each command is written from the manual. When the manual does not define a detail (for example, when `%g` switches from fixed to exponential notation), the code gets a `VERIFY` mark and a case in `compat/do/`. The actual behavior is fixed by the log produced in Stata.
- SSC packages and users' do-files may be used as test cases, respecting their licenses.
- The interface reproduces the layout and workflow of Stata 14, but the icons and visual identity are OpenDTA's own.

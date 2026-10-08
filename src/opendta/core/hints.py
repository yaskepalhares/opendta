"""Explicações de erro: o que deu errado e como resolver.

Melhoria do OpenDTA em relação ao Stata: depois da mensagem de erro
(mantida como no Stata, assim como o código r(#)), vêm uma ou duas linhas
explicando o problema no contexto do comando e sugerindo uma correção:

    . replace nota = 9 in 2
    Obs. nos. out of range
      → "in 2" refers to observation 2, but the dataset has 0 observations.
        Hint: create observations with "set obs 2" or load data with "use".
    r(198);

As explicações ficam em inglês, como o resto da interface, e podem ser
desligadas com `set hints off` (o compare.py faz isso para comparar com os
logs do Stata). Só há explicação quando ela acrescenta algo; para erros
sem regra conhecida, nada é mostrado.
"""

from __future__ import annotations

import difflib
import os
import re
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..session import Session
    from .errors import StataError


def _close(word: str, candidates, n: int = 3) -> list[str]:
    return difflib.get_close_matches(word, list(candidates), n=n, cutoff=0.6)


def _did_you_mean(word: str, candidates) -> str:
    near = _close(word, candidates)
    if not near:
        return ""
    return "did you mean " + ", ".join(f'"{x}"' for x in near) + "?"


def _join(items: list[str], limit: int = 8) -> str:
    shown = items[:limit]
    more = f" (and {len(items) - limit} more)" if len(items) > limit else ""
    return ", ".join(shown) + more


def _in_qualifier(text: str) -> str | None:
    m = re.search(r"\bin\s+([-\w/]+)", text)
    return m.group(1) if m else None


def explain(err: "StataError", s: "Session") -> list[str]:
    """Linhas de explicação para o erro (vazio se não houver o que dizer).
    A primeira é a explicação; a segunda, se houver, começa com 'Hint:'."""
    if getattr(err, "hint", None) is not None:
        return [ln for ln in err.hint.split("\n") if ln]
    msg = err.message or ""
    if not msg:
        return []
    command, text = getattr(err, "context", None) or ("", "")
    ds = s.data
    names = ds.names
    rc = err.rc

    def two(what: str, hint: str = "") -> list[str]:
        return [what] + ([f"Hint: {hint}"] if hint else [])

    # -- observações e faixas ------------------------------------------------------------
    if msg == "Obs. nos. out of range":
        q = _in_qualifier(text)
        n = ds.nobs
        what = (f'"in {q}" asks for observations that do not exist; the dataset has '
                f"{n:,} observation{'s' if n != 1 else ''}." if q else
                f"The observation numbers are outside 1 to {n:,}.")
        if n == 0:
            return two(what, 'there is no data in memory: create observations with "set obs #" '
                       'or load a dataset with "use".')
        return two(what, f'valid numbers go from 1 to {n:,} ("in f/l" means all of them); '
                   'add observations with "set obs #".')

    # -- variáveis ---------------------------------------------------------------------------
    m = re.fullmatch(r"variable (\S+) not found", msg)
    if m:
        name = m.group(1)
        if not names:
            return two(f'There is no variable "{name}": the dataset in memory has no variables.',
                       'load data with "use" or create variables with "generate".')
        return two(f'No variable is named "{name}" (or starts with "{name}").',
                   _did_you_mean(name, names) or 'type "describe" to list the variables.')
    m = re.fullmatch(r"(\S+) ambiguous abbreviation", msg)
    if m:
        ab = m.group(1)
        matches = [n for n in names if n.startswith(ab.rstrip("*~"))]
        return two(f'"{ab}" is the beginning of more than one variable: {_join(matches)}.',
                   "type more letters so that only one variable matches.")
    m = re.fullmatch(r"variable (\S+) already defined", msg)
    if m:
        name = m.group(1)
        return two(f'A variable named "{name}" already exists.',
                   f'to change its values use "replace {name} = ..."; to start over, '
                   f'"drop {name}" first.')
    m = re.fullmatch(r"(\S+) not found", msg)
    if m and rc == 111:
        name = m.group(1)
        hint = _did_you_mean(name, names + list(s.scalars))
        return two(f'"{name}" is not a variable, scalar or known name.',
                   hint or "text values must be in double quotes, e.g. \"abc\".")
    if msg == "type mismatch":
        what = "The command mixed text (string) and numbers."
        if command in ("replace", "generate", "gen", "egen"):
            return two(what, "a numeric variable only takes numbers and a string variable only "
                       'text in quotes; convert with real("12") or string(12).')
        return two(what, 'text must be in double quotes and compared with text; convert '
                   'with real() or string().')

    # -- comandos e sintaxe ------------------------------------------------------------------
    m = re.fullmatch(r"command (\S+) is unrecognized", msg)
    if m:
        from ..commands.registry import REGISTRY
        word = m.group(1)
        return two(f'"{word}" is not a command known to OpenDTA.',
                   _did_you_mean(word, REGISTRY) or "check the spelling; some commands "
                   "are not implemented in OpenDTA yet (see the ROADMAP).")
    m = re.fullmatch(r"option (\S+) not allowed", msg)
    if m:
        return two(f'"{m.group(1)}" is not an option of {command or "this command"}.',
                   "check the spelling of the options after the comma.")
    m = re.fullmatch(r"unknown function (\S+)\(\)", msg)
    if m:
        from ..lang.functions import FUNCTIONS
        name = m.group(1)
        return two(f'"{name}()" is not a function known to OpenDTA.',
                   _did_you_mean(name, FUNCTIONS) or "check the spelling of the function.")
    if msg == "varlist required":
        return two(f"{command or 'This command'} needs one or more variable names.",
                   f'for example "{command or "command"} x y"; "_all" means every variable.')
    if msg == "varlist not allowed":
        return two(f"{command or 'This command'} does not take variable names here.")
    if msg in ("invalid syntax", "unmatched quote") or msg.startswith("too many '('"):
        problems = []
        if text.count('"') % 2:
            problems.append("a double quote is not closed")
        if text.count("(") != text.count(")"):
            problems.append("parentheses are unbalanced")
        if text.count("[") != text.count("]"):
            problems.append("brackets are unbalanced")
        if text.count(",") > 1 and re.search(r",[^,]*,", text.split('"')[0] if '"' in text else text):
            problems.append("there is more than one comma (options come after a single comma)")
        if problems:
            return two("OpenDTA could not read this command: " + "; ".join(problems) + ".")
        if command:
            return two(f"OpenDTA could not read this {command} command.",
                       f'check the order of the parts: {command} [varlist] [= exp] '
                       '[if] [in] [using file] [, options].')
        return []
    if re.fullmatch(r"= ?exp(ression)? required", msg):
        return two(f"{command or 'This command'} needs an expression after an equals sign.",
                   f'for example "{command or "generate"} y = x * 2".')

    # -- dados em memória ----------------------------------------------------------------------
    if msg == "no; data in memory would be lost":
        return two("The data in memory have changes that were not saved.",
                   'save them first ("save filename") or add the clear option to discard '
                   'them, e.g. "use file, clear".')
    if msg in ("no observations", "no variables defined"):
        return two("There is no data for this command to work with.",
                   'load a dataset with "use" or "import", or create data with "input" or '
                   '"set obs".')
    if msg.startswith("observations must be at least"):
        return two("set obs can only add observations, never remove them.",
                   'to remove observations use "drop in" or "keep in".')

    # -- arquivos ------------------------------------------------------------------------------
    m = re.fullmatch(r"file (.+) not found", msg)
    if m:
        target = Path(m.group(1)).expanduser()
        folder = target.parent if str(target.parent) not in ("", ".") else Path(os.getcwd())
        where = folder if target.parent.is_absolute() else Path(os.getcwd()) / target.parent
        try:
            siblings = [p.name for p in where.iterdir()] if where.is_dir() else []
        except OSError:
            siblings = []
        near = _close(target.name, siblings)
        what = f'OpenDTA looked for "{target.name}" in {where.resolve() if where.exists() else where}.'
        if near:
            return two(what, "files with similar names there: " + ", ".join(near) + ".")
        return two(what, 'check the name and the folder; "pwd" shows the working directory and '
                   '"cd" changes it.')
    m = re.fullmatch(r"file (.+) already exists", msg)
    if m:
        return two(f'A file named "{m.group(1)}" already exists and was left untouched.',
                   'add the replace option to overwrite it, e.g. "save file, replace".')
    if msg.startswith("file") and "could not be opened" in msg:
        return two("The operating system did not allow the file to be written or read.",
                   "check that the folder exists, that you can write to it and that the file "
                   "is not open in another program.")
    if msg.startswith("file not Stata format"):
        return two("The file is not a .dta dataset, or it is damaged.",
                   'for text files use "import delimited"; for spreadsheets, "import excel".')
    if msg == "assertion is false":
        return two("The condition given to assert is not true for every observation.",
                   'use "list if !(condition)" to see the observations that break it.')
    return []

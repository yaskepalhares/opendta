"""Error explanations: what went wrong and how to fix it.

OpenDTA improvement over Stata: after the error message
(kept as in Stata, as is the r(#) code), one or two lines
explain the problem in the context of the command and suggest a fix:

    . replace nota = 9 in 2
    Obs. nos. out of range
      → "in 2" refers to observation 2, but the dataset has 0 observations.
        Hint: create observations with "set obs 2" or load data with "use".
    r(198);

The explanations come out in the interface language (English by default; Portuguese
with language pt, see i18n.py) and can be turned off with `set hints off` (compare.py does
this to compare with the Stata logs). There is an explanation only when it adds something; for
errors without a known rule, nothing is shown.
"""

from __future__ import annotations

import difflib
import os
import re
from pathlib import Path
from typing import TYPE_CHECKING

from ..i18n import L

if TYPE_CHECKING:
    from ..session import Session
    from .errors import StataError


def _close(word: str, candidates, n: int = 3) -> list[str]:
    return difflib.get_close_matches(word, list(candidates), n=n, cutoff=0.6)


def _did_you_mean(word: str, candidates) -> str:
    near = _close(word, candidates)
    if not near:
        return ""
    return L("did you mean ", "você quis dizer ") + ", ".join(f'"{x}"' for x in near) + "?"


def _join(items: list[str], limit: int = 8) -> str:
    shown = items[:limit]
    more = (L(f" (and {len(items) - limit} more)", f" (e mais {len(items) - limit})")
            if len(items) > limit else "")
    return ", ".join(shown) + more


# commands that do not follow the standard syntax [varlist] [=exp] [if] [in] ...
_FREE_SYNTAX = {"matrix", "display", "local", "global", "file", "program", "forvalues",
                "foreach", "while", "if", "scalar", "set", "label", "notes", "char",
                "syntax", "gettoken", "tokenize", "args", "import", "export", "macro",
                "return", "ereturn", "sreturn", "adopath", "sysdir", "infix", "infile"}


def _in_qualifier(text: str) -> str | None:
    m = re.search(r"\bin\s+([-\w/]+)", text)
    return m.group(1) if m else None


def explain(err: "StataError", s: "Session") -> list[str]:
    """Explanation lines for the error (empty if there is nothing to say).
    The first is the explanation; the second, if any, starts with 'Hint:'."""
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
        return [what] + ([L("Hint: ", "Dica: ") + hint] if hint else [])

    # -- observations and ranges ------------------------------------------------------------
    if msg == "Obs. nos. out of range":
        q = _in_qualifier(text)
        n = ds.nobs
        plural = n != 1
        what = (L(f'"in {q}" asks for observations that do not exist; the dataset has '
                  f"{n:,} observation{'s' if plural else ''}.",
                  f'"in {q}" pede observações que não existem; a base tem '
                  f"{n:,} observaç{'ões' if plural else 'ão'}.") if q else
                L(f"The observation numbers are outside 1 to {n:,}.",
                  f"Os números de observação estão fora de 1 a {n:,}."))
        if n == 0:
            return two(what, L('there is no data in memory: create observations with "set obs #" '
                               'or load a dataset with "use".',
                               'não há dados na memória: crie observações com "set obs #" '
                               'ou abra uma base com "use".'))
        return two(what, L(f'valid numbers go from 1 to {n:,} ("in f/l" means all of them); '
                           'add observations with "set obs #".',
                           f'os números válidos vão de 1 a {n:,} ("in f/l" quer dizer todas); '
                           'acrescente observações com "set obs #".'))

    # -- variables ---------------------------------------------------------------------------
    m = re.fullmatch(r"variable (\S+) not found", msg)
    if m:
        name = m.group(1)
        if not names:
            return two(L(f'There is no variable "{name}": the dataset in memory has no variables.',
                         f'Não existe a variável "{name}": a base na memória não tem variáveis.'),
                       L('load data with "use" or create variables with "generate".',
                         'abra uma base com "use" ou crie variáveis com "generate".'))
        return two(L(f'No variable is named "{name}" (or starts with "{name}").',
                     f'Nenhuma variável se chama "{name}" (nem começa com "{name}").'),
                   _did_you_mean(name, names) or L('type "describe" to list the variables.',
                                                   'digite "describe" para listar as variáveis.'))
    m = re.fullmatch(r"(\S+) ambiguous abbreviation", msg)
    if m:
        ab = m.group(1)
        matches = [n for n in names if n.startswith(ab.rstrip("*~"))]
        return two(L(f'"{ab}" is the beginning of more than one variable: {_join(matches)}.',
                     f'"{ab}" é o começo de mais de uma variável: {_join(matches)}.'),
                   L("type more letters so that only one variable matches.",
                     "digite mais letras, até sobrar só uma variável."))
    m = re.fullmatch(r"variable (\S+) already defined", msg)
    if m:
        name = m.group(1)
        return two(L(f'A variable named "{name}" already exists.',
                     f'Já existe uma variável chamada "{name}".'),
                   L(f'to change its values use "replace {name} = ..."; to start over, '
                     f'"drop {name}" first.',
                     f'para mudar os valores use "replace {name} = ..."; para recomeçar, '
                     f'"drop {name}" antes.'))
    m = re.fullmatch(r"(\S+) not found", msg)
    if m and rc == 111:
        name = m.group(1)
        hint = _did_you_mean(name, names + list(s.scalars))
        return two(L(f'"{name}" is not a variable, scalar or known name.',
                     f'"{name}" não é variável, escalar nem nome conhecido.'),
                   hint or L("text values must be in double quotes, e.g. \"abc\".",
                             "textos vão entre aspas duplas, por exemplo \"abc\"."))
    if msg == "type mismatch":
        what = L("The command mixed text (string) and numbers.",
                 "O comando misturou texto (string) e números.")
        if command in ("replace", "generate", "gen", "egen"):
            return two(what, L("a numeric variable only takes numbers and a string variable only "
                               'text in quotes; convert with real("12") or string(12).',
                               "variável numérica só recebe números e variável string só texto "
                               'entre aspas; converta com real("12") ou string(12).'))
        return two(what, L('text must be in double quotes and compared with text; convert '
                           'with real() or string().',
                           'texto vai entre aspas duplas e se compara com texto; converta '
                           'com real() ou string().'))

    # -- comandos e sintaxe ------------------------------------------------------------------
    m = re.fullmatch(r"command (\S+) is unrecognized", msg)
    if m:
        from ..commands.registry import REGISTRY
        word = m.group(1)
        return two(L(f'"{word}" is not a command known to OpenDTA.',
                     f'"{word}" não é um comando conhecido pelo OpenDTA.'),
                   _did_you_mean(word, REGISTRY) or L(
                       "check the spelling; some commands are not implemented in OpenDTA "
                       "yet (see the ROADMAP).",
                       "confira a grafia; alguns comandos ainda não existem no OpenDTA "
                       "(veja o ROADMAP)."))
    m = re.fullmatch(r"option (\S+) not allowed", msg)
    if m:
        return two(L(f'"{m.group(1)}" is not an option of {command or "this command"}.',
                     f'"{m.group(1)}" não é opção de {command or "este comando"}.'),
                   L("check the spelling of the options after the comma.",
                     "confira a grafia das opções depois da vírgula."))
    m = re.fullmatch(r"unknown function (\S+)\(\)", msg)
    if m:
        from ..lang.functions import FUNCTIONS
        name = m.group(1)
        return two(L(f'"{name}()" is not a function known to OpenDTA.',
                     f'"{name}()" não é uma função conhecida pelo OpenDTA.'),
                   _did_you_mean(name, FUNCTIONS) or L("check the spelling of the function.",
                                                       "confira a grafia da função."))
    if msg == "varlist required":
        return two(L(f"{command or 'This command'} needs one or more variable names.",
                     f"{command or 'Este comando'} precisa de um ou mais nomes de variáveis."),
                   L(f'for example "{command or "command"} x y"; "_all" means every variable.',
                     f'por exemplo "{command or "comando"} x y"; "_all" quer dizer todas as variáveis.'))
    if msg == "varlist not allowed":
        return two(L(f"{command or 'This command'} does not take variable names here.",
                     f"{command or 'Este comando'} não aceita nomes de variáveis aqui."))
    if msg in ("invalid syntax", "unmatched quote") or msg.startswith("too many '('"):
        problems = []
        if text.count('"') % 2:
            problems.append(L("a double quote is not closed", "há aspas duplas sem fechar"))
        if text.count("(") != text.count(")"):
            problems.append(L("parentheses are unbalanced", "os parênteses não fecham"))
        if text.count("[") != text.count("]"):
            problems.append(L("brackets are unbalanced", "os colchetes não fecham"))
        if text.count(",") > 1 and re.search(r",[^,]*,", text.split('"')[0] if '"' in text else text):
            problems.append(L("there is more than one comma (options come after a single comma)",
                              "há mais de uma vírgula (as opções vêm depois de uma vírgula só)"))
        if problems:
            return two(L("OpenDTA could not read this command: ",
                         "O OpenDTA não conseguiu ler este comando: ") + "; ".join(problems) + ".")
        if command in _FREE_SYNTAX:
            return two(L(f"OpenDTA could not read this {command} command.",
                         f"O OpenDTA não conseguiu ler este comando {command}."),
                       L("check quotes, parentheses and the order of its parts.",
                         "confira aspas, parênteses e a ordem das partes."))
        if command:
            return two(L(f"OpenDTA could not read this {command} command.",
                         f"O OpenDTA não conseguiu ler este comando {command}."),
                       L(f'check the order of the parts: {command} [varlist] [= exp] '
                         '[if] [in] [using file] [, options].',
                         f'confira a ordem das partes: {command} [varlist] [= exp] '
                         '[if] [in] [using arquivo] [, opções].'))
        return []
    if re.fullmatch(r"= ?exp(ression)? required", msg):
        return two(L(f"{command or 'This command'} needs an expression after an equals sign.",
                     f"{command or 'Este comando'} precisa de uma expressão depois do sinal de igual."),
                   L(f'for example "{command or "generate"} y = x * 2".',
                     f'por exemplo "{command or "generate"} y = x * 2".'))

    # -- data in memory ----------------------------------------------------------------------
    if msg == "no; data in memory would be lost":
        return two(L("The data in memory have changes that were not saved.",
                     "Os dados na memória têm alterações que não foram salvas."),
                   L('save them first ("save filename") or add the clear option to discard '
                     'them, e.g. "use file, clear".',
                     'salve antes ("save arquivo") ou acrescente a opção clear para descartá-las, '
                     'por exemplo "use arquivo, clear".'))
    if msg in ("no observations", "no variables defined"):
        return two(L("There is no data for this command to work with.",
                     "Não há dados para este comando trabalhar."),
                   L('load a dataset with "use" or "import", or create data with "input" or '
                     '"set obs".',
                     'abra uma base com "use" ou "import", ou crie dados com "input" ou '
                     '"set obs".'))
    if msg.startswith("observations must be at least"):
        return two(L("set obs can only add observations, never remove them.",
                     "set obs só acrescenta observações, nunca remove."),
                   L('to remove observations use "drop in" or "keep in".',
                     'para remover observações use "drop in" ou "keep in".'))

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
        place = where.resolve() if where.exists() else where
        what = L(f'OpenDTA looked for "{target.name}" in {place}.',
                 f'O OpenDTA procurou "{target.name}" em {place}.')
        if near:
            return two(what, L("files with similar names there: ",
                               "arquivos com nomes parecidos lá: ") + ", ".join(near) + ".")
        return two(what, L('check the name and the folder; "pwd" shows the working directory and '
                           '"cd" changes it.',
                           'confira o nome e a pasta; "pwd" mostra a pasta de trabalho e '
                           '"cd" muda.'))
    m = re.fullmatch(r"file (.+) already exists", msg)
    if m:
        return two(L(f'A file named "{m.group(1)}" already exists and was left untouched.',
                     f'Já existe um arquivo chamado "{m.group(1)}", e ele não foi alterado.'),
                   L('add the replace option to overwrite it, e.g. "save file, replace".',
                     'acrescente a opção replace para substituí-lo, por exemplo "save arquivo, replace".'))
    if msg.startswith("file") and "could not be opened" in msg:
        return two(L("The operating system did not allow the file to be written or read.",
                     "O sistema operacional não deixou gravar ou ler o arquivo."),
                   L("check that the folder exists, that you can write to it and that the file "
                     "is not open in another program.",
                     "confira se a pasta existe, se você pode gravar nela e se o arquivo "
                     "não está aberto em outro programa."))
    if msg.startswith("file not Stata format"):
        return two(L("The file is not a .dta dataset, or it is damaged.",
                     "O arquivo não é uma base .dta, ou está danificado."),
                   L('for text files use "import delimited"; for spreadsheets, "import excel".',
                     'para arquivos de texto use "import delimited"; para planilhas, "import excel".'))
    if msg == "assertion is false":
        return two(L("The condition given to assert is not true for every observation.",
                     "A condição do assert não vale para todas as observações."),
                   L('use "list if !(condition)" to see the observations that break it.',
                     'use "list if !(condição)" para ver as observações que a violam.'))
    return []

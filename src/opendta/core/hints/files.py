"""Explanations of errors with files: not found (with similar names, or the same name in another
format), already existing, not writable, not a dataset."""

from __future__ import annotations

import os
import re
from pathlib import Path

from ...i18n import L
from .base import Ctx, Explanation, close, rule

_IMPORTERS = {".csv": "import delimited", ".txt": "import delimited", ".tsv": "import delimited",
              ".xlsx": "import excel", ".xls": "import excel", ".sav": "import spss",
              ".sas7bdat": "import sas"}


@rule(r"file (.+) not found")
def _not_found(ctx: Ctx, m) -> Explanation:
    target = Path(m.group(1).strip('"')).expanduser()
    where = target.parent if target.parent.is_absolute() else Path(os.getcwd()) / target.parent
    try:
        siblings = [p.name for p in where.iterdir()] if where.is_dir() else []
    except OSError:
        siblings = []
    place = where.resolve() if where.exists() else where
    what = L(f'OpenDTA looked for "{target.name}" in {place}, and it is not there.',
             f'O OpenDTA procurou "{target.name}" em {place}, e ele não está lá.')
    why = L('Without a full path, files are looked for in the working directory (shown by "pwd" '
            "and at the bottom of the window). When the name has no extension, .dta is added.",
            'Sem o caminho completo, os arquivos são procurados na pasta de trabalho (mostrada pelo '
            '"pwd" e no rodapé da janela). Quando o nome não tem extensão, o .dta é acrescentado.')
    other = [n for n in siblings if Path(n).stem == target.stem and Path(n).suffix.lower() in _IMPORTERS]
    if other and ctx.command in ("use", "merge", "append", ""):
        imp = _IMPORTERS[Path(other[0]).suffix.lower()]
        return Explanation(
            what + L(f' There is "{other[0]}", which is not a Stata dataset.',
                     f' Existe "{other[0]}", que não é uma base do Stata.'),
            L("use, merge and append read .dta files; other formats are read with import and can "
              "then be saved as .dta.",
              "use, merge e append leem arquivos .dta; outros formatos são lidos com import e "
              "depois podem ser salvos como .dta."),
            f'{imp} "{other[0]}", clear', "import")
    near = close(target.name, siblings)
    fix = (L("files with similar names there: ", "arquivos com nomes parecidos lá: ") + ", ".join(near) + "."
           if near else
           L('check the name and the folder; "cd" changes the working directory (or File > Change '
             "working directory), and \"dir\" lists the files.",
             'confira o nome e a pasta; o "cd" muda a pasta de trabalho (ou Arquivo > Mudar pasta '
             'de trabalho), e o "dir" lista os arquivos.'))
    return Explanation(what, why, fix, "files")


@rule(r"file (.+) already exists")
def _exists(ctx: Ctx, m) -> Explanation:
    name = m.group(1)
    return Explanation(
        L(f'A file named "{name}" already exists and was left untouched.',
          f'Já existe um arquivo chamado "{name}", e ele não foi alterado.'),
        L("Stata never overwrites a file unless you ask, so that a typo cannot erase your work.",
          "O Stata nunca sobrescreve um arquivo sem você pedir, para que um erro de digitação não "
          "apague seu trabalho."),
        L(f'add the replace option to overwrite it: "{ctx.text.strip().split(",")[0]}, replace".',
          f'acrescente a opção replace para substituí-lo: "{ctx.text.strip().split(",")[0]}, replace".'))


@rule(r"file .* could not be opened.*")
def _not_opened(ctx: Ctx, _m) -> Explanation:
    return Explanation(
        L("The operating system did not allow the file to be written or read.",
          "O sistema operacional não deixou gravar ou ler o arquivo."),
        L("This happens when the folder does not exist, when you cannot write to it, or when the "
          "file is open in another program (Excel locks the files it has open).",
          "Isso acontece quando a pasta não existe, quando você não pode gravar nela, ou quando o "
          "arquivo está aberto em outro programa (o Excel trava os arquivos que tem abertos)."),
        L("check the folder, close the file in other programs, or save it under another name.",
          "confira a pasta, feche o arquivo nos outros programas, ou salve com outro nome."),
        "files")


@rule(r"file not Stata format.*")
def _not_dta(ctx: Ctx, _m) -> Explanation:
    return Explanation(
        L("The file is not a Stata dataset (.dta), or it is damaged.",
          "O arquivo não é uma base do Stata (.dta), ou está danificado."),
        L("use reads only .dta files. Spreadsheets and text files have other formats and are read "
          "with import.",
          "O use lê só arquivos .dta. Planilhas e arquivos de texto têm outros formatos e são "
          "lidos com import."),
        L('for text files (.csv, .txt) use "import delimited"; for spreadsheets, "import excel".',
          'para arquivos de texto (.csv, .txt) use "import delimited"; para planilhas, "import excel".'),
        "import")


@rule(r"invalid file specification")
def _bad_spec(ctx: Ctx, _m) -> Explanation | None:
    if not re.search(r"\s", ctx.args.split(",")[0].strip()):
        return None
    return Explanation(
        L("The file name has spaces, so it must be in double quotes.",
          "O nome do arquivo tem espaços, então precisa ficar entre aspas duplas."),
        "",
        L(f'for example: {ctx.command} "my file.dta"', f'por exemplo: {ctx.command} "meu arquivo.dta"'),
        "files")

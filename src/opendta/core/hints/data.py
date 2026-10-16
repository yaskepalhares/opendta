"""Explanations of errors about the data in memory: observations, variables and their types,
sorting, merging and assert."""

from __future__ import annotations

import re

import numpy as np

from ...i18n import L
from .base import Ctx, Explanation, did_you_mean, fmt_value, join, obs, rule

# ---------------------------------------------------------------------------
# observations
# ---------------------------------------------------------------------------


@rule(r"Obs\. nos\. out of range")
def _out_of_range(ctx: Ctx, _m) -> Explanation:
    q = re.search(r"\bin\s+([-\w/]+)", ctx.text)
    n = ctx.ds.nobs
    what = (L(f'"in {q.group(1)}" asks for observations that do not exist; the dataset has {obs(n)}.',
              f'"in {q.group(1)}" pede observações que não existem; a base tem {obs(n)}.') if q else
            L(f"The observation numbers are outside 1 to {n:,}.",
              f"Os números de observação estão fora de 1 a {n:,}."))
    why = L('Observations are numbered from 1 to the number of observations (_N). "in" takes one '
            'number or a range: "in 5", "in 1/10"; "f" is the first and "l" the last, and '
            '"-5/l" means the last five.',
            'As observações são numeradas de 1 até o número de observações (_N). O "in" recebe um '
            'número ou um intervalo: "in 5", "in 1/10"; "f" é a primeira, "l" a última, e '
            '"-5/l" quer dizer as cinco últimas.')
    if n == 0:
        return Explanation(what, why, L('there is no data in memory: create observations with '
                                        '"set obs #" or load a dataset with "use".',
                                        'não há dados na memória: crie observações com "set obs #" '
                                        'ou abra uma base com "use".'))
    return Explanation(what, why, L(f'valid numbers go from 1 to {n:,} ("in f/l" means all of them); '
                                    'add observations with "set obs #".',
                                    f'os números válidos vão de 1 a {n:,} ("in f/l" quer dizer '
                                    'todas); acrescente observações com "set obs #".'))


@rule(r"no; data in memory would be lost")
def _would_be_lost(ctx: Ctx, _m) -> Explanation:
    return Explanation(
        L("The data in memory have changes that were not saved, and this command would replace "
          "them.",
          "Os dados na memória têm alterações não salvas, e este comando os substituiria."),
        L("Stata keeps one dataset in memory at a time. Before opening another one, it asks you to "
          "save the changes or to say that they can be discarded.",
          "O Stata guarda uma base de cada vez na memória. Antes de abrir outra, pede que você "
          "salve as alterações ou diga que elas podem ser descartadas."),
        L('save them first ("save filename, replace") or add the clear option to discard them, '
          'e.g. "use file, clear".',
          'salve antes ("save arquivo, replace") ou acrescente a opção clear para descartá-las, '
          'por exemplo "use arquivo, clear".'),
        "use")


@rule(r"no observations|no variables defined")
def _no_data(ctx: Ctx, _m) -> Explanation | None:
    if ctx.ds.nobs and ctx.ds.nvars:
        return None                  # the data exist: estimation.py explains the empty sample
    return Explanation(
        L("There is no data in memory for this command to work with.",
          "Não há dados na memória para este comando trabalhar."),
        L("Commands work on the dataset in memory; it is empty now (after clear, or before any "
          "use).",
          "Os comandos trabalham sobre a base que está na memória; ela está vazia agora (depois "
          "de um clear, ou antes de qualquer use)."),
        L('load a dataset with "use" or "import", or create data with "input" or "set obs".',
          'abra uma base com "use" ou "import", ou crie dados com "input" ou "set obs".'),
        "use")


@rule(r"(?:obs|observations) must be at least .*")
def _set_obs(ctx: Ctx, _m) -> Explanation:
    return Explanation(
        L(f"set obs can only add observations: the dataset already has {obs(ctx.ds.nobs)}.",
          f"O set obs só acrescenta observações: a base já tem {obs(ctx.ds.nobs)}."),
        "",
        L('to remove observations use "drop in" or "keep in" (for example "keep in 1/100").',
          'para remover observações use "drop in" ou "keep in" (por exemplo "keep in 1/100").'),
        "datamgmt")


# ---------------------------------------------------------------------------
# variables
# ---------------------------------------------------------------------------

_WHY_NAMES = (
    "Variable names must be written exactly, with upper and lower case as in the dataset; Stata "
    "also accepts the beginning of a name when only one variable starts that way.",
    "Os nomes das variáveis precisam ser escritos exatamente, com maiúsculas e minúsculas como "
    "na base; o Stata também aceita o começo de um nome quando só uma variável começa assim.")

_CMP = r"(==|!=|~=|<=|>=|<|>)"


def _compared_with(ctx: Ctx, word: str) -> str | None:
    """The variable that `word` is compared with in the command (x == word or word == x)."""
    for pat in (rf"\b(\w+)\s*{_CMP}\s*{re.escape(word)}\b", rf"\b{re.escape(word)}\s*{_CMP}\s*(\w+)\b"):
        m = re.search(pat, ctx.text)
        if m:
            other = m.group(1) if m.group(1) != word and not m.group(1).startswith(("=", "!", "<", ">", "~")) \
                else m.group(2)
            found = ctx.variables(other)
            if found:
                return found[0]
    return None


def _value_meant(ctx: Ctx, word: str) -> Explanation | None:
    """`if sex == male` with male a value, not a variable: quotes or the value label."""
    var = _compared_with(ctx, word)
    if var is None:
        return None
    v = ctx.ds.get(var)
    why = L(f"Without quotes, Stata reads {word} as the name of a variable. Text values go in "
            "double quotes; a numeric variable with value labels is compared with its numeric "
            'code, or with "label":labelname.',
            f"Sem aspas, o Stata lê {word} como o nome de uma variável. Valores de texto vão entre "
            "aspas duplas; uma variável numérica com rótulos de valor se compara com o código "
            'numérico, ou com "rótulo":nomedorotulo.')
    if v.is_string:
        return Explanation(
            L(f'"{word}" is not a variable; {var} is a string (text) variable, so its values go in '
              "quotes.",
              f'"{word}" não é uma variável; {var} é uma variável de texto (string), então os '
              "valores vão entre aspas."),
            why, re.sub(rf"\b{re.escape(word)}\b", f'"{word}"', ctx.text.strip(), count=1), "operators")
    labels = ctx.ds.value_labels.get(v.value_label, {}) if v.value_label else {}
    codes = [k for k, lab in labels.items() if lab.lower() == word.lower()]
    if codes:
        code = fmt_value(codes[0])
        return Explanation(
            L(f'"{word}" is not a variable: it is the value label of code {code} of {var} '
              f"(label {v.value_label}).",
              f'"{word}" não é uma variável: é o rótulo do código {code} de {var} (rótulo '
              f"{v.value_label})."),
            why,
            L(f'compare with the code: "... {var} == {code}", or with the label: '
              f'"... {var} == "{labels[codes[0]]}":{v.value_label}".',
              f'compare com o código: "... {var} == {code}", ou com o rótulo: '
              f'"... {var} == "{labels[codes[0]]}":{v.value_label}".'),
            "label")
    return None


@rule(r"variable (\S+) not found")
def _variable_not_found(ctx: Ctx, m) -> Explanation:
    name = m.group(1)
    if not ctx.names:
        return Explanation(
            L(f'There is no variable "{name}": the dataset in memory has no variables.',
              f'Não existe a variável "{name}": a base na memória não tem variáveis.'),
            L("Commands work on the dataset in memory, which is empty now.",
              "Os comandos trabalham sobre a base que está na memória, que está vazia agora."),
            L('load data with "use" or create variables with "generate".',
              'abra uma base com "use" ou crie variáveis com "generate".'), "use")
    meant = _value_meant(ctx, name)
    if meant is not None:
        return meant
    return Explanation(
        L(f'No variable is named "{name}" (or starts with "{name}").',
          f'Nenhuma variável se chama "{name}" (nem começa com "{name}").'),
        L(*_WHY_NAMES),
        did_you_mean(name, ctx.names) or did_you_mean(name.lower(), [n.lower() for n in ctx.names])
        or L('type "describe" to list the variables.', 'digite "describe" para listar as variáveis.'),
        "describe")


@rule(r"(\S+) not found", rc=111)
def _name_not_found(ctx: Ctx, m) -> Explanation | None:
    name = m.group(1)
    if " " in name or name.startswith(("help", "file")):
        return None
    meant = _value_meant(ctx, name)
    if meant is not None:
        return meant
    return Explanation(
        L(f'"{name}" is not a variable, scalar or other known name.',
          f'"{name}" não é variável, escalar nem outro nome conhecido.'),
        L("A word without quotes inside an expression is read as a name (of a variable, scalar "
          "or matrix); text values must be in double quotes.",
          "Uma palavra sem aspas dentro de uma expressão é lida como nome (de variável, escalar "
          "ou matriz); valores de texto precisam ficar entre aspas duplas."),
        did_you_mean(name, ctx.names + list(ctx.s.scalars))
        or L(f'if it is a text value, write "{name}" with the quotes.',
             f'se for um valor de texto, escreva "{name}" com as aspas.'))


@rule(r"(\S+) ambiguous abbreviation")
def _ambiguous(ctx: Ctx, m) -> Explanation:
    ab = m.group(1)
    matches = [n for n in ctx.names if n.startswith(ab.rstrip("*~"))]
    return Explanation(
        L(f'"{ab}" is the beginning of more than one variable: {join(matches)}.',
          f'"{ab}" é o começo de mais de uma variável: {join(matches)}.'),
        L("Stata accepts a shortened variable name only when a single variable starts that way.",
          "O Stata aceita o nome abreviado só quando uma única variável começa assim."),
        L("type more letters, or the full name, so that only one variable matches.",
          "digite mais letras, ou o nome inteiro, até sobrar só uma variável."))


@rule(r"variable (\S+) already defined")
def _already_defined(ctx: Ctx, m) -> Explanation:
    name = m.group(1)
    if name == "_merge":
        return Explanation(
            L("A previous merge left the variable _merge in the data, and merge needs to create it "
              "again.",
              "Um merge anterior deixou a variável _merge na base, e o merge precisa criá-la de "
              "novo."),
            L("_merge tells where each observation came from (1 only in memory, 2 only in the "
              "using file, 3 in both). Check it with tabulate and then drop it.",
              "_merge diz de onde veio cada observação (1 só na memória, 2 só no arquivo using, 3 "
              "nos dois). Confira com tabulate e depois apague."),
            L('run "tabulate _merge", then "drop _merge", and merge again (or use the '
              "generate(newname) option of merge).",
              'rode "tabulate _merge", depois "drop _merge", e faça o merge de novo (ou use a '
              "opção generate(novonome) do merge)."),
            "merge")
    return Explanation(
        L(f'A variable named "{name}" already exists.', f'Já existe uma variável chamada "{name}".'),
        L("generate only creates new variables; replace changes the values of a variable that "
          "already exists.",
          "O generate só cria variáveis novas; o replace muda os valores de uma variável que já "
          "existe."),
        L(f'to change its values use "replace {name} = ..."; to start over, "drop {name}" first; '
          "or choose another name.",
          f'para mudar os valores use "replace {name} = ..."; para recomeçar, "drop {name}" antes; '
          "ou escolha outro nome."),
        "generate")


# ---------------------------------------------------------------------------
# strings and numbers
# ---------------------------------------------------------------------------

_WHY_TYPES = (
    "Stata has two kinds of variables: numeric (numbers) and string (text). They cannot be "
    "mixed in a calculation or comparison, and text values are always written in double quotes.",
    "O Stata tem dois tipos de variáveis: numéricas (números) e string (texto). Elas não se "
    "misturam num cálculo ou comparação, e valores de texto sempre vão entre aspas duplas.")


def _type_mismatch_detail(ctx: Ctx) -> Explanation | None:
    text = ctx.text.strip()
    # numeric variable compared with text in quotes
    for m in re.finditer(rf'\b(\w+)\s*{_CMP}\s*"([^"]*)"', text):
        var = (ctx.variables(m.group(1)) or [None])[0]
        if var and not ctx.is_string(var):
            fixed = text.replace(m.group(0), f"{m.group(1)} {m.group(2)} {m.group(3)}"
                                 if re.fullmatch(r"-?\d+(\.\d*)?", m.group(3)) else m.group(0), 1)
            what = L(f'{var} is numeric, but "{m.group(3)}" is text in quotes.',
                     f'{var} é numérica, mas "{m.group(3)}" é texto entre aspas.')
            # the text is the value label of a code: compare with the code
            lab = ctx.ds.get(var).value_label
            codes = [k for k, v in ctx.ds.value_labels.get(lab, {}).items() if v == m.group(3)] if lab else []
            if codes:
                fixed = text.replace(m.group(0), f"{m.group(1)} {m.group(2)} {codes[0]}", 1)
                what = L(f'{var} is numeric: "{m.group(3)}" is the value label of code {codes[0]}.',
                         f'{var} é numérica: "{m.group(3)}" é o rótulo de valor do código {codes[0]}.')
                fix = L(f'{fixed}  (or {m.group(1)} {m.group(2)} "{m.group(3)}":{lab})',
                        f'{fixed}  (ou {m.group(1)} {m.group(2)} "{m.group(3)}":{lab})')
                return Explanation(what, L(*_WHY_TYPES), fix, "operators")
            fix = (fixed if fixed != text else
                   L(f'compare {var} with a number; if it has value labels, use "label":labelname.',
                     f'compare {var} com um número; se ela tiver rótulos de valor, use '
                     '"rótulo":nomedorotulo.'))
            return Explanation(what, L(*_WHY_TYPES), fix, "operators")
    # string variable compared with a number
    for m in re.finditer(rf"\b(\w+)\s*{_CMP}\s*(-?\d+(?:\.\d*)?)(?!\w)", text):
        var = (ctx.variables(m.group(1)) or [None])[0]
        if var and ctx.is_string(var):
            fixed = text.replace(m.group(0), f'{m.group(1)} {m.group(2)} "{m.group(3)}"', 1)
            return Explanation(
                L(f"{var} is a string (text) variable, but {m.group(3)} is a number.",
                  f"{var} é uma variável de texto (string), mas {m.group(3)} é um número."),
                L(*_WHY_TYPES),
                L(f"compare with text: {fixed}  (or, if {var} holds numbers written as text, "
                  f'convert it first with "destring {var}, replace")',
                  f"compare com texto: {fixed}  (ou, se {var} guarda números escritos como texto, "
                  f'converta antes com "destring {var}, replace")'),
                "operators")
    # generate/replace: the target and the expression
    if ctx.command in ("generate", "replace") and "=" in ctx.args:
        target = ctx.args.split("=")[0].split()[-1]
        exp_vars = ctx.variables(ctx.args.split("=", 1)[1].split(" if ")[0])
        strs = [v for v in exp_vars if ctx.is_string(v)]
        nums = [v for v in exp_vars if not ctx.is_string(v)]
        if strs and nums:
            return Explanation(
                L(f"The expression mixes numeric variables ({join(nums)}) with string variables "
                  f"({join(strs)}).",
                  f"A expressão mistura variáveis numéricas ({join(nums)}) com variáveis de texto "
                  f"({join(strs)})."),
                L(*_WHY_TYPES),
                L(f'convert the text to numbers with real({strs[0]}) (or "destring {strs[0]}, '
                  f'replace"), or the numbers to text with string({nums[0]}).',
                  f'converta o texto em números com real({strs[0]}) (ou "destring {strs[0]}, '
                  f'replace"), ou os números em texto com string({nums[0]}).'),
                "destring")
        if ctx.command == "replace" and ctx.ds.has(target):
            if ctx.is_string(target):
                return Explanation(
                    L(f"{target} is a string (text) variable, but the new value is a number.",
                      f"{target} é uma variável de texto (string), mas o novo valor é um número."),
                    L(*_WHY_TYPES),
                    L(f'put the value in quotes (replace {target} = "...") or convert {target} to '
                      f'numbers first with "destring {target}, replace".',
                      f'coloque o valor entre aspas (replace {target} = "...") ou converta {target} '
                      f'em números antes com "destring {target}, replace".'),
                    "destring")
            return Explanation(
                L(f"{target} is numeric, but the new value is text.",
                  f"{target} é numérica, mas o novo valor é texto."),
                L(*_WHY_TYPES),
                L(f'give it a number, or create a string variable instead: generate {target}_s = "...".',
                  f'dê um número, ou crie uma variável de texto: generate {target}_s = "...".'),
                "generate")
    return None


@rule(r"type mismatch")
def _type_mismatch(ctx: Ctx, _m) -> Explanation:
    detail = _type_mismatch_detail(ctx)
    if detail is not None:
        return detail
    return Explanation(
        L("The command mixed text (string) and numbers.", "O comando misturou texto (string) e números."),
        L(*_WHY_TYPES),
        L('convert with real("12") (text to number) or string(12) (number to text); "describe" '
          "shows which variables are strings (str#).",
          'converta com real("12") (texto para número) ou string(12) (número para texto); o '
          '"describe" mostra quais variáveis são texto (str#).'),
        "destring")


@rule(r"string variables not allowed in varlist;\s*(\S+) is a string variable")
def _string_in_varlist(ctx: Ctx, m) -> Explanation:
    var = m.group(1)
    x = ctx.ds.get(var).data if ctx.ds.has(var) else np.array([], dtype=object)
    vals = [v for v in x[:200] if isinstance(v, str) and v != ""]
    digits = bool(vals) and all(re.fullmatch(r"\s*-?[\d.,]+\s*", v) for v in vals)
    fix = (L(f'its values look like numbers stored as text: convert with "destring {var}, replace".',
             f'os valores parecem números guardados como texto: converta com "destring {var}, '
             'replace".') if digits else
           L(f'for categories, create a numeric version with "encode {var}, generate({var}_n)" and '
             f"use i.{var}_n in the model.",
             f'para categorias, crie uma versão numérica com "encode {var}, generate({var}_n)" e use '
             f"i.{var}_n no modelo."))
    return Explanation(
        L(f"{var} is a string (text) variable, and {ctx.command or 'this command'} needs numbers.",
          f"{var} é uma variável de texto (string), e o {ctx.command or 'comando'} precisa de "
          "números."),
        L("Statistics and models are computed on numeric variables. Text must be converted first: "
          "destring when the text is made of digits, encode when it holds categories.",
          "Estatísticas e modelos são calculados sobre variáveis numéricas. O texto precisa ser "
          "convertido antes: destring quando o texto é feito de números, encode quando guarda "
          "categorias."),
        fix, "recode")


@rule(r"not possible with (numeric|string) variable")
def _not_possible(ctx: Ctx, m) -> Explanation:
    kind = m.group(1)
    num = kind == "numeric"
    return Explanation(
        L(f"{ctx.command or 'This command'} does not work with a {kind} variable.",
          f"O {ctx.command or 'comando'} não funciona com uma variável "
          f"{'numérica' if num else 'de texto'}."),
        L(*_WHY_TYPES),
        L('"describe" shows the type of each variable; "tostring" and "destring" convert '
          "between the two.",
          'o "describe" mostra o tipo de cada variável; "tostring" e "destring" convertem de um '
          "para o outro."),
        "recode")


# ---------------------------------------------------------------------------
# sorting, by and merge
# ---------------------------------------------------------------------------

@rule(r"not sorted")
def _not_sorted(ctx: Ctx, _m) -> Explanation:
    m = re.match(r"\s*by\s+([^:,]+)", ctx.text)
    byvars = m.group(1).strip() if m else "group"
    rest = ctx.text.split(":", 1)[1].strip() if ":" in ctx.text else "..."
    return Explanation(
        L(f"by {byvars}: needs the data sorted by {byvars}, and they are not.",
          f"O by {byvars}: precisa dos dados ordenados por {byvars}, e eles não estão."),
        L("by repeats the command for each group of consecutive observations with the same value, "
          "so the observations of a group must be next to each other. bysort sorts and runs in "
          "one step.",
          "O by repete o comando para cada grupo de observações seguidas com o mesmo valor; por "
          "isso as observações de um grupo precisam estar juntas. O bysort ordena e roda de uma "
          "vez."),
        f"bysort {byvars}: {rest}", "by")


@rule(r"variables? (.+?) (?:does|do) not uniquely identify (?:the )?observations(?: in the (master|using) data)?")
def _not_unique(ctx: Ctx, m) -> Explanation:
    keys = m.group(1).split()
    where = m.group(2) or "master"
    n_dup = ""
    if where == "master" and all(ctx.ds.has(k) for k in keys) and ctx.ds.nobs:
        cols = [ctx.ds.get(k).data for k in keys]
        _, counts = np.unique(np.array(list(zip(*[c.tolist() for c in cols])), dtype=object).astype(str),
                              axis=0, return_counts=True)
        n_dup = L(f" ({int((counts > 1).sum()):,} key values appear more than once)",
                  f" ({int((counts > 1).sum()):,} valores da chave aparecem mais de uma vez)")
    data = L("in memory" if where == "master" else "in the using file",
             "na memória" if where == "master" else "no arquivo using")
    k = " ".join(keys)
    return Explanation(
        L(f"In the data {data}, {k} repeats across observations{n_dup}.",
          f"Nos dados {data}, {k} se repete entre observações{n_dup}."),
        L("1:1 merge needs each value of the key once in each dataset. When one side has several "
          "observations per key (people in a household, years of a country), use m:1 (repeated in "
          "memory) or 1:m (repeated in the using file).",
          "O merge 1:1 precisa de cada valor da chave uma única vez em cada base. Quando um lado "
          "tem várias observações por chave (pessoas num domicílio, anos de um país), use m:1 "
          "(repetida na memória) ou 1:m (repetida no arquivo using)."),
        L(f'see the repeated values with "duplicates report {k}" and "duplicates list {k}"; '
          "then choose m:1 or 1:m, or fix the key.",
          f'veja os valores repetidos com "duplicates report {k}" e "duplicates list {k}"; '
          "depois escolha m:1 ou 1:m, ou corrija a chave."),
        "merge")


# ---------------------------------------------------------------------------
# assert
# ---------------------------------------------------------------------------

@rule(r"assertion is false")
def _assert(ctx: Ctx, _m) -> Explanation:
    cond = ctx.args.split(",")[0]
    p = ctx.parsed()
    what = L("The condition given to assert is not true for every observation.",
             "A condição do assert não vale para todas as observações.")
    if p is not None and p.exp is None:
        try:
            from ...commands._util import eval_vector
            from ...lang.vexpr import truth
            body = re.split(r"\b(if|in)\b", cond)[0].strip()
            v = np.broadcast_to(np.asarray(truth(eval_vector(ctx.s, body)), dtype=bool), (ctx.ds.nobs,))
            bad = np.flatnonzero(ctx.sample() & ~v) + 1
            if len(bad):
                what = L(f"The condition is false in {obs(len(bad))} of {ctx.ds.nobs:,} (for example "
                         f"observation{'s' if len(bad) > 1 else ''} {join([str(i) for i in bad[:5]], 5)}).",
                         f"A condição é falsa em {obs(len(bad))} de {ctx.ds.nobs:,} (por exemplo "
                         f"{'as observações' if len(bad) > 1 else 'a observação'} "
                         f"{join([str(i) for i in bad[:5]], 5)}).")
                body_ = body
                return Explanation(
                    what,
                    L("assert stops a do-file when the data are not as expected. Remember that a "
                      "missing value (.) counts as larger than any number, so x > 0 is true when x "
                      "is missing.",
                      "O assert interrompe o do-file quando os dados não são como o esperado. "
                      "Lembre que o valor faltante (.) conta como maior que qualquer número: x > 0 "
                      "é verdadeiro quando x está faltando."),
                    L(f'list the observations that break it: "list if !({body_})".',
                      f'liste as observações que a violam: "list if !({body_})".'),
                    "assert")
        except Exception:  # noqa: BLE001 - fall back to the general text
            pass
    return Explanation(
        what, "",
        L('use "list if !(condition)" to see the observations that break it.',
          'use "list if !(condição)" para ver as observações que a violam.'),
        "assert")


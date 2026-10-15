"""Explanations of errors in the way a command is written: "=" for "==", quotes and parentheses
left open, unknown commands, functions and options, invalid names. The catch-all for
"invalid syntax" is the last rule of the package (this module is imported last)."""

from __future__ import annotations

import re
import unicodedata
from functools import lru_cache

from ...i18n import L
from .base import Ctx, Explanation, did_you_mean, join, rule

# ---------------------------------------------------------------------------
# the syntax line of a command, from its help page
# ---------------------------------------------------------------------------


@lru_cache(maxsize=None)
def _syntax_lines(path: str) -> tuple[str, ...]:
    from ..smcl import render_text
    try:
        src = open(path, encoding="utf-8").read()
    except OSError:
        return ()
    m = re.search(r"\{title:Syntax\}(.*?)(\{title:|\Z)", src, re.S)
    if not m:
        return ()
    text = render_text(m.group(1), 200)
    return tuple(" ".join(ln.split()) for ln in text.splitlines() if ln.strip())


def syntax_of(ctx: Ctx, command: str | None = None) -> str:
    """The syntax line of the command in its help page ("" if there is none)."""
    from ...commands.helpcmd import find_help
    cmd = command or ctx.command
    if not cmd:
        return ""
    path = find_help(ctx.s, cmd)
    lines = _syntax_lines(str(path)) if path else ()
    for ln in lines:
        if ln.split()[0] == cmd:
            return ln
    return lines[0] if lines else ""


def _order() -> str:
    return L("the command, then the variables, = expression, if, in, and the options after a single "
             "comma",
             "o comando, depois as variáveis, = expressão, if, in, e as opções depois de uma única "
             "vírgula")

# ---------------------------------------------------------------------------
# "=" where "==" is meant (and the opposite)
# ---------------------------------------------------------------------------

_SINGLE_EQ = re.compile(r"(?<![=<>!~])=(?!=)")


def _if_clause(text: str) -> tuple[int, int] | None:
    """Start and end of the condition after `if` (outside quotes, up to the options' comma)."""
    m = re.search(r"\bif\b", text)
    if not m:
        return None
    start, depth, quote = m.end(), 0, False
    for i in range(start, len(text)):
        ch = text[i]
        if ch == '"':
            quote = not quote
        elif not quote and ch in "([":
            depth += 1
        elif not quote and ch in ")]":
            depth -= 1
        elif not quote and depth == 0 and (ch == "," or text.startswith(" in ", i)):
            return start, i
    return start, len(text)


def _unquoted(text: str) -> str:
    return re.sub(r'"[^"]*"', lambda m: "_" * len(m.group(0)), text)


@rule(r"invalid syntax")
def _equals_in_if(ctx: Ctx, _m) -> Explanation | None:
    span = _if_clause(_unquoted(ctx.text))
    if span is None:
        return None
    a, b = span
    cond = ctx.text[a:b]
    if not _SINGLE_EQ.search(_unquoted(cond)):
        return None
    fixed = ctx.text[:a] + _SINGLE_EQ.sub("==", cond) + ctx.text[b:]
    return Explanation(
        L(f'The condition after "if" ({cond.strip()}) uses a single "=".',
          f'A condição depois do "if" ({cond.strip()}) usa um "=" só.'),
        L('In Stata, a single "=" gives a value (as in "generate y = 1"); to ask whether two '
          'things are equal, use a double "==". The other comparisons are != (different), <, >, '
          '<= and >=.',
          'No Stata, um "=" só atribui um valor (como em "generate y = 1"); para perguntar se '
          'duas coisas são iguais, use dois: "==". As outras comparações são != (diferente), <, '
          '>, <= e >=.'),
        fixed.strip(), "operators")


@rule(r"= ?exp(ression)? required")
def _exp_required(ctx: Ctx, _m) -> Explanation:
    cmd = ctx.command or "generate"
    if "==" in ctx.args.split(" if ")[0]:
        fixed = ctx.text.replace("==", "=", 1)
        return Explanation(
            L(f'"==" compares two values; to give a value to a variable, {cmd} uses a single "=".',
              f'"==" compara dois valores; para dar um valor a uma variável, o {cmd} usa um "=" só.'),
            L('"=" assigns and "==" compares: "generate y = 1" creates y, while "if y == 1" asks '
              'whether y is 1.',
              '"=" atribui e "==" compara: "generate y = 1" cria y, enquanto "if y == 1" pergunta '
              'se y vale 1.'),
            fixed.strip(), "operators")
    return Explanation(
        L(f"{cmd} needs an equals sign followed by an expression (the value to store).",
          f"O {cmd} precisa de um sinal de igual seguido de uma expressão (o valor a guardar)."),
        L("An expression is any calculation with numbers, variables and functions, such as "
          "x * 2, ln(income) or (age > 60).",
          "Uma expressão é qualquer cálculo com números, variáveis e funções, como x * 2, "
          "ln(renda) ou (idade > 60)."),
        L(f'for example "{cmd} y = x * 2".', f'por exemplo "{cmd} y = x * 2".'))


# ---------------------------------------------------------------------------
# quotes, parentheses and commas
# ---------------------------------------------------------------------------

def _balance_problems(text: str) -> list[str]:
    problems = []
    if text.count('"') % 2:
        problems.append(L("a double quote is not closed", "há aspas duplas sem fechar"))
    bare = _unquoted(text) if text.count('"') % 2 == 0 else text
    if bare.count("(") != bare.count(")"):
        problems.append(L("parentheses are unbalanced", "os parênteses não fecham"))
    if bare.count("[") != bare.count("]"):
        problems.append(L("brackets are unbalanced", "os colchetes não fecham"))
    if bare.count(",") > 1 and re.search(r",[^,()]*,", re.sub(r"\([^()]*\)", "", bare)):
        problems.append(L("there is more than one comma (options come after a single comma)",
                          "há mais de uma vírgula (as opções vêm depois de uma vírgula só)"))
    return problems


@rule(r"invalid syntax|unmatched quote|too many '\('.*|parentheses unbalanced")
def _balance(ctx: Ctx, _m) -> Explanation | None:
    problems = _balance_problems(ctx.text)
    if not problems:
        return None
    return Explanation(
        L("OpenDTA could not read this command: ", "O OpenDTA não conseguiu ler este comando: ")
        + "; ".join(problems) + ".",
        L("Every opening quote, parenthesis or bracket needs its closing one on the same line; "
          "otherwise the rest of the line is read as part of it.",
          "Toda aspa, parêntese ou colchete aberto precisa ser fechado na mesma linha; senão, o "
          "resto da linha é lido como parte dele."),
        L("count the quotes and parentheses; text values look like \"abc\".",
          "conte as aspas e os parênteses; textos ficam assim: \"abc\"."))


# ---------------------------------------------------------------------------
# common slips that make a standard command unreadable
# ---------------------------------------------------------------------------

_GEN_LIKE = {"generate", "egen", "replace"}
_TYPES = re.compile(r"^(byte|int|long|float|double|str\d+|strL)$")


@rule(r"invalid syntax")
def _space_in_new_name(ctx: Ctx, _m) -> Explanation | None:
    if ctx.command not in _GEN_LIKE or "=" not in ctx.args:
        return None
    words = [w for w in ctx.args.split("=")[0].split() if not _TYPES.match(w)]
    if len(words) < 2:
        return None
    name = "_".join(words)
    return Explanation(
        L(f'"{" ".join(words)}" has a space: a variable name must be a single word.',
          f'"{" ".join(words)}" tem espaço: o nome de uma variável precisa ser uma palavra só.'),
        L("Names are made of letters, digits and _ (underscore), start with a letter or _, and "
          "have up to 32 characters. Use the variable label for a longer description.",
          "Nomes têm letras, números e _ (sublinhado), começam com letra ou _ e têm até 32 "
          "caracteres. Para uma descrição mais longa, use o rótulo da variável."),
        L(f'join the words with _: "{ctx.command} {name} = ..."; then, if you want, '
          f'label variable {name} "{" ".join(words)}".',
          f'junte as palavras com _: "{ctx.command} {name} = ..."; depois, se quiser, '
          f'label variable {name} "{" ".join(words)}".'),
        "generate")


@rule(r"invalid syntax")
def _empty_if(ctx: Ctx, _m) -> Explanation | None:
    if not re.search(r"\bif\s*(,|$)", ctx.text.strip()):
        return None
    return Explanation(
        L('"if" is not followed by a condition.', 'Falta a condição depois do "if".'),
        L('"if" keeps only the observations where a condition is true, so it needs one, such '
          'as "if age >= 18" or "if sex == 1".',
          'O "if" mantém só as observações em que uma condição é verdadeira; por isso precisa de '
          'uma, como "if idade >= 18" ou "if sexo == 1".'),
        L(f'write the condition after "if": "{ctx.text.strip()} age >= 18".',
          f'escreva a condição depois do "if": "{ctx.text.strip()} idade >= 18".'),
        "operators")


@rule(r"invalid syntax|too few variables specified")
def _rename(ctx: Ctx, _m) -> Explanation | None:
    if ctx.command != "rename" or len(ctx.args.split()) >= 2:
        return None
    old = ctx.args.split()[0] if ctx.args.split() else "old"
    return Explanation(
        L("rename needs two names: the current one and the new one.",
          "O rename precisa de dois nomes: o atual e o novo."),
        "",
        L(f'"rename {old} newname".', f'"rename {old} novonome".'))


@rule(r"invalid syntax")
def _unquoted_path(ctx: Ctx, _m) -> Explanation | None:
    m = re.search(r"\b(using|use|save|saveold|do|run|cd)\s+([^\",]+\s[^\",]+)", ctx.text)
    if not m or '"' in ctx.text:
        return None
    path = m.group(2).strip()
    if not re.search(r"[/\\.]|\s", path):
        return None
    word = m.group(1)
    fixed = ctx.text.replace(path, f'"{path}"', 1)
    return Explanation(
        L(f'The file name "{path}" has spaces and no quotes, so Stata reads it as several words.',
          f'O nome de arquivo "{path}" tem espaços e não está entre aspas; o Stata lê como várias '
          'palavras.'),
        L(f"Put file and folder names that have spaces in double quotes after {word}.",
          f"Coloque nomes de arquivos e pastas com espaço entre aspas duplas depois de {word}."),
        fixed.strip(), "files")


# ---------------------------------------------------------------------------
# unknown commands and functions
# ---------------------------------------------------------------------------

# names of other statistical software and the Stata command that does the same
_OTHER_COMMANDS = {
    "summary": "summarize", "describe_data": "describe", "print": "list or display",
    "head": "list in 1/5", "tail": "list in -5/l", "str": "describe", "glimpse": "describe",
    "nrow": "count", "length": "count", "mutate": "generate", "transform": "generate",
    "compute": "generate", "filter": "keep if ...", "subset": "keep if ...", "select": "keep",
    "arrange": "sort", "order_by": "sort", "lm": "regress", "ols": "regress", "anova_test": "anova",
    "t.test": "ttest", "chisq": "tabulate x y, chi2", "crosstab": "tabulate x y",
    "crosstabs": "tabulate x y", "freq": "tabulate", "frequencies": "tabulate",
    "frequency": "tabulate", "descriptives": "summarize", "cor": "correlate", "corr_test": "pwcorr",
    "median": "summarize, detail", "quantile": "summarize, detail", "setwd": "cd", "getwd": "pwd",
    "view": "browse", "unique": "levelsof", "boxplot": "graph box", "plot": "twoway scatter",
    "barplot": "graph bar", "hist": "histogram", "read": "use", "load": "use", "write": "save",
    "read_csv": "import delimited", "read_excel": "import excel", "write_csv": "export delimited",
    "rm": "drop (variables) or erase (files)", "del": "drop", "delete": "drop",
    "library": "(nothing: every command is always available)", "require": "(nothing)",
}


@rule(r"command (\S+) is unrecognized")
def _unknown_command(ctx: Ctx, m) -> Explanation:
    from ...commands.registry import REGISTRY
    word = m.group(1)
    rest = ctx.text.strip()[len(word):].lstrip()
    if rest.startswith("=") and not rest.startswith("=="):
        verb = "replace" if ctx.ds.has(word) else "generate"
        return Explanation(
            L(f'"{word}" is not a command: to {"change" if verb == "replace" else "create"} a '
              f"variable, the line must start with {verb}.",
              f'"{word}" não é um comando: para {"mudar" if verb == "replace" else "criar"} uma '
              f"variável, a linha precisa começar com {verb}."),
            L("Every Stata line starts with a command. generate creates a new variable; replace "
              "changes the values of one that exists.",
              "Toda linha do Stata começa com um comando. generate cria uma variável nova; "
              "replace muda os valores de uma que já existe."),
            f"{verb} {ctx.text.strip()}", "generate")
    if word.lower() in _OTHER_COMMANDS:
        alt = _OTHER_COMMANDS[word.lower()]
        return Explanation(
            L(f'"{word}" comes from another program; Stata\'s equivalent is {alt}.',
              f'"{word}" vem de outro programa; o equivalente no Stata é {alt}.'),
            L("Each statistical package has its own command names.",
              "Cada programa estatístico tem seus próprios nomes de comando."),
            L(f'try "{alt.split(" or ")[0]}" (type "help {alt.split()[0]}" for its options).',
              f'tente "{alt.split(" or ")[0]}" (digite "help {alt.split()[0]}" para ver as opções).'),
            alt.split()[0] if alt.split()[0].isalpha() else "")
    near = did_you_mean(word, REGISTRY)
    return Explanation(
        L(f'"{word}" is not a command known to OpenDTA.',
          f'"{word}" não é um comando conhecido pelo OpenDTA.'),
        L("The first word of each line is the command. Commands can be shortened (sum for "
          "summarize, gen for generate), but the letters must be in the right order.",
          "A primeira palavra de cada linha é o comando. Os comandos podem ser abreviados (sum "
          "para summarize, gen para generate), mas as letras precisam estar na ordem certa."),
        near or L('check the spelling; type "help" to see the commands, and note that some '
                  "Stata commands are not in OpenDTA yet (see the ROADMAP).",
                  'confira a grafia; digite "help" para ver os comandos. Alguns comandos do Stata '
                  "ainda não existem no OpenDTA (veja o ROADMAP)."),
        "opendta")


# functions of egen that people try in generate
_EGEN_ONLY = {"mean", "median", "sd", "count", "total", "rowmean", "rowtotal", "rowmax", "rowmin",
              "rowsd", "rowmiss", "rownonmiss", "rank", "group", "tag", "mode", "pctile", "iqr",
              "std", "cut", "concat", "fill", "seq", "diff", "kurt", "skew", "mad", "mdev",
              "anycount", "anymatch", "anyvalue", "rowfirst", "rowlast", "rowmedian", "rowpctile"}

_OTHER_FUNCTIONS = {
    "nchar": "strlen()", "len": "strlen()", "toupper": "strupper()", "upper": "strupper()",
    "tolower": "strlower()", "lower": "strlower()", "paste": "the + operator (\"a\" + \"b\")",
    "paste0": "the + operator", "concat": "the + operator", "is_na": "missing()",
    "isna": "missing()", "isnull": "missing()", "is_null": "missing()", "nvl": "cond()",
    "ifelse": "cond()", "if_else": "cond()", "iif": "cond()", "as_numeric": "real()",
    "as_character": "string()", "str": "string()", "num": "real()", "number": "real()",
    "grepl": "regexm() or strpos()", "contains": "strpos()", "today": "date(c(current_date), \"DMY\")",
    "log2": "ln(x)/ln(2)", "power": "the ^ operator", "pow": "the ^ operator",
}


@rule(r"unknown function (\S+)\(\)")
def _unknown_function(ctx: Ctx, m) -> Explanation:
    from ...lang.functions import FUNCTIONS
    name = m.group(1)
    if name in _EGEN_ONLY and ctx.command in ("generate", "replace"):
        fixed = re.sub(r"^\s*(gen\w*|replace)\b", "egen", ctx.text.strip())
        return Explanation(
            L(f'"{name}()" is an egen function: generate cannot use it.',
              f'"{name}()" é uma função do egen: o generate não sabe usá-la.'),
            L("generate computes each observation from the values of that same observation. "
              "Functions that look at many observations at once (mean, total, rank...) or "
              "across variables (rowmean, rowtotal...) belong to egen.",
              "O generate calcula cada observação a partir dos valores dela mesma. Funções que "
              "olham várias observações de uma vez (mean, total, rank...) ou várias variáveis "
              "(rowmean, rowtotal...) são do egen."),
            fixed if ctx.command == "generate" else L(
                f"create a new variable with egen (egen tmp = {name}(...)) and then replace.",
                f"crie uma variável nova com egen (egen tmp = {name}(...)) e depois use replace."),
            "egen")
    if name.lower() in _OTHER_FUNCTIONS or name.lower().replace(".", "_") in _OTHER_FUNCTIONS:
        alt = _OTHER_FUNCTIONS.get(name.lower(), _OTHER_FUNCTIONS.get(name.lower().replace(".", "_")))
        return Explanation(
            L(f'"{name}()" comes from another program; in Stata use {alt}.',
              f'"{name}()" vem de outro programa; no Stata use {alt}.'),
            L("Function names differ between statistical packages.",
              "Os nomes das funções mudam de um programa estatístico para outro."),
            re.sub(rf"\b{re.escape(name)}\(", alt.split()[0][:-1], ctx.text.strip(), count=1)
            if alt.endswith("()") and " " not in alt else "", "functions")
    return Explanation(
        L(f'"{name}()" is not a function known to OpenDTA.',
          f'"{name}()" não é uma função conhecida pelo OpenDTA.'),
        L("Function names are lowercase and are followed directly by parentheses, as in "
          "ln(x) or substr(name, 1, 3).",
          "Os nomes de funções são em minúsculas e vêm colados aos parênteses, como em ln(x) ou "
          "substr(nome, 1, 3)."),
        did_you_mean(name, FUNCTIONS) or L("check the spelling of the function.",
                                            "confira a grafia da função."),
        "functions")


# ---------------------------------------------------------------------------
# options, variable lists and names
# ---------------------------------------------------------------------------

@rule(r"option (\S+) not allowed")
def _bad_option(ctx: Ctx, m) -> Explanation:
    word = m.group(1)
    known = ctx.err.detail.get("options") or []
    near = did_you_mean(word, known) if known else ""
    cmd = ctx.command or L("this command", "este comando")
    fix = near or (L(f"the options of {cmd} are: ", f"as opções do {cmd} são: ") + join(sorted(known), 12)
                   + "." if known else L("check the spelling of the options after the comma.",
                                         "confira a grafia das opções depois da vírgula."))
    return Explanation(
        L(f'"{word}" is not an option of {cmd}.', f'"{word}" não é uma opção do {cmd}.'),
        L("Options come after the comma, and each command has its own list; most can be "
          "shortened (d for detail), but not misspelled.",
          "As opções vêm depois da vírgula, e cada comando tem a sua lista; a maioria pode ser "
          "abreviada (d para detail), mas não escrita com erro."),
        fix)


@rule(r"varlist required")
def _varlist_required(ctx: Ctx, _m) -> Explanation:
    cmd = ctx.command or L("This command", "Este comando")
    return Explanation(
        L(f"{cmd} needs one or more variable names.", f"O {cmd} precisa de um ou mais nomes de variáveis."),
        L('A varlist is a list of variables separated by spaces; "_all" means every variable, '
          '"x*" every variable that starts with x and "x1-x5" the variables from x1 to x5 in the '
          "dataset order.",
          'Uma varlist é uma lista de variáveis separadas por espaço; "_all" quer dizer todas, '
          '"x*" todas que começam com x e "x1-x5" as variáveis de x1 até x5 na ordem da base.'),
        L(f'for example "{ctx.command or "command"} x y".', f'por exemplo "{ctx.command or "comando"} x y".'))


@rule(r"varlist not allowed")
def _varlist_not_allowed(ctx: Ctx, _m) -> Explanation:
    cmd = ctx.command or L("This command", "Este comando")
    syn = syntax_of(ctx)
    return Explanation(
        L(f"{cmd} does not take variable names here.", f"O {cmd} não aceita nomes de variáveis aqui."),
        "",
        (L("its syntax is: ", "a sintaxe dele é: ") + syn) if syn else "")


_RESERVED = {"_all", "_b", "byte", "_coef", "_cons", "double", "float", "if", "in", "int", "long",
             "_n", "_N", "_pi", "_pred", "_rc", "_se", "_skip", "strL", "using", "with"}


def _name_problem(name: str) -> tuple[str, str]:
    """(what is wrong with the name, a valid name close to it)."""
    plain = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()
    fixed = re.sub(r"[^A-Za-z0-9_]", "_", plain)[:32] or "v"
    if fixed[0].isdigit():
        fixed = "_" + fixed
    if name in _RESERVED or re.fullmatch(r"str\d+", name):
        return L(f'"{name}" is a reserved word in Stata', f'"{name}" é uma palavra reservada do Stata'), name + "_"
    if name[:1].isdigit():
        return L("it starts with a digit", "começa com um número"), fixed
    if len(name) > 32:
        return L("it has more than 32 characters", "tem mais de 32 caracteres"), fixed
    bad = sorted({c for c in name if not (c.isascii() and (c.isalnum() or c == "_"))})
    if bad:
        return (L("it has characters that names cannot have: ", "tem caracteres que nomes não podem ter: ")
                + " ".join(bad)), fixed
    return L("it is not a valid name", "não é um nome válido"), fixed


@rule(r"(\S+) invalid name")
def _invalid_name(ctx: Ctx, m) -> Explanation:
    name = m.group(1)
    problem, fixed = _name_problem(name)
    return Explanation(
        L(f'"{name}" cannot be a name: {problem}.', f'"{name}" não pode ser um nome: {problem}.'),
        L("Names of variables, scalars and matrices use only letters without accents, digits and "
          "_; they start with a letter or _ and have up to 32 characters. Upper and lower case "
          "are different names.",
          "Nomes de variáveis, escalares e matrizes usam só letras sem acento, números e _; "
          "começam com letra ou _ e têm até 32 caracteres. Maiúsculas e minúsculas são nomes "
          "diferentes."),
        L(f'use "{fixed}" instead.', f'use "{fixed}" no lugar.'))


@rule(r"help for (\S+) not found")
def _no_help(ctx: Ctx, m) -> Explanation:
    from ...commands.helpcmd import available_topics
    topic = m.group(1)
    return Explanation(
        L(f'"{topic}" has no help page in OpenDTA.', f'"{topic}" não tem página de ajuda no OpenDTA.'),
        "",
        did_you_mean(topic, available_topics()) or L('type "help" to see the topics available.',
                                                     'digite "help" para ver os tópicos disponíveis.'),
        "opendta")


# ---------------------------------------------------------------------------
# the catch-all: invalid syntax of a standard command
# ---------------------------------------------------------------------------

@rule(r"invalid syntax")
def _invalid_syntax(ctx: Ctx, _m) -> Explanation | None:
    if not ctx.command:
        return None
    syn = syntax_of(ctx)
    return Explanation(
        L(f"OpenDTA could not read this {ctx.command} command.",
          f"O OpenDTA não conseguiu ler este comando {ctx.command}."),
        L(f"Commands are written in a fixed order: {_order()}. A word in the wrong place, a "
          "missing comma before the options or an extra symbol stops the reading.",
          f"Os comandos seguem uma ordem fixa: {_order()}. Uma palavra fora do lugar, a falta da "
          "vírgula antes das opções ou um símbolo a mais interrompem a leitura."),
        (L("compare with the syntax: ", "compare com a sintaxe: ") + syn) if syn else
        L("check the order of the parts and the spelling of the options.",
          "confira a ordem das partes e a grafia das opções."))

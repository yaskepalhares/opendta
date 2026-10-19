"""Explanations of errors in loops, blocks, programs and Mata."""

from __future__ import annotations

import re

from ...i18n import L
from .base import Ctx, Explanation, rule

_LOOP_FORMS = (
    'foreach v of varlist x y z {   (variables)  |  foreach w in a b c {   (any words)  |  '
    'forvalues i = 1/10 {   (numbers). The "{" ends the first line, the commands go on the '
    'next lines and "}" closes the block on a line of its own.',
    'foreach v of varlist x y z {   (variáveis)  |  foreach w in a b c {   (quaisquer palavras)  |  '
    'forvalues i = 1/10 {   (números). A "{" termina a primeira linha, os comandos vêm nas '
    'linhas seguintes e a "}" fecha o bloco numa linha só para ela.')


@rule(r"invalid syntax|command (?:foreach|forvalues) is unrecognized|.*\{ required|program error.*")
def _loops(ctx: Ctx, _m) -> Explanation | None:
    word = ctx.text.strip().split(None, 1)[0] if ctx.text.strip() else ""
    if ctx.command not in ("foreach", "forvalues") and word not in ("foreach", "forvalues"):
        return None
    loop = ctx.command or word
    return Explanation(
        L(f"This {loop} line is not written in one of the forms Stata knows.",
          f"Esta linha de {loop} não está escrita em uma das formas que o Stata conhece."),
        L("The loop forms are: ", "As formas do laço são: ") + L(*_LOOP_FORMS),
        L("inside the block, refer to the loop variable as `v' (backtick before, apostrophe after), "
          "for example: summarize `v'.",
          "dentro do bloco, use a variável do laço como `v' (crase antes, apóstrofo depois), por "
          "exemplo: summarize `v'."),
        "macro")


@rule(r"unexpected end of file")
def _unclosed_block(ctx: Ctx, _m) -> Explanation:
    opened = ctx.text.count("{") - ctx.text.count("}")
    return Explanation(
        L(f"A block opened with {{ was not closed{f' ({opened} still open)' if opened > 1 else ''}.",
          f"Um bloco aberto com {{ não foi fechado{f' ({opened} ainda abertos)' if opened > 1 else ''}."),
        L("foreach, forvalues, while, if/else and program open a block with { at the end of the "
          "line; Stata runs it only when it finds the } that closes it, on a line of its own.",
          "foreach, forvalues, while, if/else e program abrem um bloco com { no fim da linha; o "
          "Stata só roda quando encontra a } que o fecha, numa linha só para ela."),
        L("add the closing }. In the Command window, type the whole block (Shift+Enter starts a "
          "new line) or write it in the Do-file Editor.",
          "acrescente a } de fechamento. Na janela Command, digite o bloco inteiro (Shift+Enter "
          "começa uma nova linha) ou escreva no Do-file Editor."),
        "do")


@rule(r"unexpected \}")
def _extra_brace(ctx: Ctx, _m) -> Explanation:
    return Explanation(
        L("This } closes a block that was never opened.", "Esta } fecha um bloco que nunca foi aberto."),
        L("Each { at the end of a foreach, forvalues, if or program line pairs with one }.",
          "Cada { no fim de uma linha de foreach, forvalues, if ou program forma par com uma }."),
        L("remove the extra }, or check that the { of the block is the last character of its line.",
          "apague a } a mais, ou confira se a { do bloco é o último caractere da linha dela."),
        "do")


@rule(r"else without if")
def _else(ctx: Ctx, _m) -> Explanation:
    return Explanation(
        L("This else has no if block right before it.", "Este else não tem um bloco if logo antes."),
        L("In do-files, if and else are written as blocks, if `x' > 0 { ... } else { ... }, and "
          "the else comes right after the } of the if.",
          "Em do-files, if e else são escritos como blocos, if `x' > 0 { ... } else { ... }, e o "
          "else vem logo depois da } do if."),
        L('to choose values per observation, use "if" after a command (replace y = 1 if x > 0) or '
          "the cond() function, not if/else.",
          'para escolher valores por observação, use o "if" depois de um comando (replace y = 1 if '
          "x > 0) ou a função cond(), não if/else."),
        "do")


@rule(r"(\S+) already defined", rc=110)
def _program_defined(ctx: Ctx, m) -> Explanation | None:
    if ctx.command != "program":
        return None
    name = m.group(1)
    return Explanation(
        L(f"A program named {name} already exists in memory (probably from an earlier run of this "
          "do-file).",
          f"Já existe um programa chamado {name} na memória (provavelmente de uma execução "
          "anterior deste do-file)."),
        L("Programs stay in memory until they are dropped or OpenDTA is closed.",
          "Os programas ficam na memória até serem apagados ou o OpenDTA ser fechado."),
        L(f'put "capture program drop {name}" right before "program define {name}".',
          f'coloque "capture program drop {name}" logo antes de "program define {name}".'),
        "program")


@rule(r"label (\S+) already defined")
def _label_defined(ctx: Ctx, m) -> Explanation:
    name = m.group(1)
    return Explanation(
        L(f"The value label {name} already exists.", f"O rótulo de valores {name} já existe."),
        "",
        L(f'add "modify" to add codes to it or "replace" to start it over: '
          f'label define {name} 1 "..." , modify',
          f'acrescente "modify" para incluir códigos ou "replace" para refazê-lo: '
          f'label define {name} 1 "..." , modify'),
        "label")


# ---------------------------------------------------------------------------
# Mata (the message is the Mata traceback: "      <istmt>:  3200  conformability error")
# ---------------------------------------------------------------------------

@rule(rc=3200)
def _mata_conformability(ctx: Ctx, _m) -> Explanation:
    return Explanation(
        L("In Mata, the matrices of this operation do not have compatible sizes.",
          "No Mata, as matrizes desta operação não têm tamanhos compatíveis."),
        L("+ and - need the same rows and columns; * needs the columns of the first equal to the "
          "rows of the second; the operators with a colon (:*, :+, :/) work element by element, "
          "also between a matrix and a row or column vector.",
          "+ e - precisam das mesmas linhas e colunas; * precisa que as colunas da primeira sejam "
          "iguais às linhas da segunda; os operadores com dois-pontos (:*, :+, :/) trabalham "
          "elemento a elemento, também entre uma matriz e um vetor linha ou coluna."),
        L("check the sizes with rows(X) and cols(X); X' transposes.",
          "confira os tamanhos com rows(X) e cols(X); X' transpõe."),
        "mata")


@rule(rc=3250)
def _mata_type(ctx: Ctx, _m) -> Explanation:
    return Explanation(
        L("In Mata, this operation mixed real numbers and strings.",
          "No Mata, esta operação misturou números reais e textos."),
        L("Mata matrices are either real or string, and functions expect one kind.",
          "As matrizes do Mata são de números ou de texto, e cada função espera um tipo."),
        L("convert with strtoreal() and strofreal(); eltype(X) shows the type.",
          "converta com strtoreal() e strofreal(); eltype(X) mostra o tipo."),
        "mata")


@rule(rc=3499)
def _mata_not_found(ctx: Ctx, _m) -> Explanation:
    m = re.search(r"(\w+) not found", ctx.msg)
    name = m.group(1) if m else "x"
    return Explanation(
        L(f"In Mata, {name} does not exist yet.", f"No Mata, {name} ainda não existe."),
        L("Mata variables exist only after a value is assigned to them; Stata variables are not "
          "Mata variables.",
          "Variáveis do Mata só existem depois que recebem um valor; variáveis do Stata não são "
          "variáveis do Mata."),
        L(f'assign it first ({name} = ...); to bring a Stata variable, {name} = st_data(., "{name}").',
          f'atribua antes ({name} = ...); para trazer uma variável do Stata, '
          f'{name} = st_data(., "{name}").'),
        "mata")


@rule(r"unrecognized command", rc=199)
def _plugin_call(ctx: Ctx, _m) -> Explanation | None:
    if ctx.command != "plugin" and not ctx.text.strip().startswith("plugin"):
        return None
    return Explanation(
        L("OpenDTA does not run compiled Stata plugins, so there is no plugin to call.",
          "O OpenDTA não executa plugins compilados do Stata, então não há plugin para chamar."),
        L("A plugin is a file of machine code written for Stata (program name, plugin loads it and "
          "plugin call runs it); OpenDTA cannot load it.",
          "Um plugin é um arquivo de código de máquina escrito para o Stata (program nome, plugin "
          "o carrega e plugin call o executa); o OpenDTA não consegue carregá-lo."),
        L("do the same computation in a do-file, an ado-file or Mata.",
          "faça o mesmo cálculo em um do-file, um ado-file ou em Mata."),
        "progutil")

"""Explanations of errors in statistics and models: empty or small samples (with the missing
values of each variable), outcomes that do not vary, panel and time settings, estimates that
do not exist yet and matrices that do not conform."""

from __future__ import annotations

import re

import numpy as np

from ...i18n import L
from ..missing import SYSMISS
from .base import Ctx, Explanation, and_join, join, obs, rule, values_summary

_WHY_SAMPLE = (
    "Statistics and models use only the observations that have a value in every variable of "
    "the command (and that pass its if and in). A variable missing (.) in many observations "
    "shrinks that sample; a model also needs more observations than coefficients.",
    "Estatísticas e modelos usam só as observações que têm valor em todas as variáveis do "
    "comando (e que passam pelo if e pelo in). Uma variável faltante (.) em muitas observações "
    "encolhe essa amostra; um modelo também precisa de mais observações do que coeficientes.")


def _missing(ctx: Ctx, name: str, rows: np.ndarray) -> np.ndarray:
    x = ctx.ds.get(name).data[rows]
    if ctx.is_string(name):
        return np.array([v == "" for v in x], dtype=bool)
    return np.asarray(x, dtype=float) >= SYSMISS


def _model_vars(ctx: Ctx) -> list[str]:
    p = ctx.parsed()
    if p is None:
        return []
    return ctx.variables(re.sub(r"\b[icbo]\d*\.|\b[LDFS]\d*\.", "", p.varlist))


@rule(r"no observations|insufficient observations")
def _empty_sample(ctx: Ctx, _m) -> Explanation | None:
    if not ctx.ds.nobs or not ctx.ds.nvars:
        return None
    p = ctx.parsed()
    rows = ctx.sample()
    n_sel = int(rows.sum())
    vars_ = _model_vars(ctx)
    cond = " ".join(x for x in ((f"if {p.if_}" if p and p.if_ else ""), (f"in {p.in_}" if p and p.in_ else "")) if x)
    if cond and n_sel == 0:
        return Explanation(
            L(f'No observation satisfies "{cond}", so nothing is left to analyze.',
              f'Nenhuma observação satisfaz "{cond}", então não sobra nada para analisar.'),
            L(*_WHY_SAMPLE),
            L(f'check the condition with "count {cond}" and the values with "tabulate"; text '
              'values must match exactly, including upper case and spaces.',
              f'confira a condição com "count {cond}" e os valores com "tabulate"; valores de '
              "texto precisam bater exatamente, inclusive maiúsculas e espaços."),
            "missing")
    if not vars_:
        return None
    idx = np.flatnonzero(rows)
    miss = {v: _missing(ctx, v, idx) for v in vars_}
    complete = int((~np.any(np.column_stack(list(miss.values())), axis=1)).sum()) if idx.size else 0
    worst = sorted(((int(m.sum()), v) for v, m in miss.items() if m.any()), reverse=True)
    strings = [v for v in vars_ if ctx.is_string(v)]
    scope = L(f" selected by {cond}", f" selecionadas por {cond}") if cond else ""
    allv = and_join(vars_)
    if complete == 0:
        what = L(f"No observation{scope} has a value in all of {allv}",
                 f"Nenhuma observação{scope.replace('selecionadas', 'selecionada')} tem valor em "
                 f"todas as variáveis {allv}")
    else:
        what = L(f"Only {complete:,} of the {n_sel:,} observations{scope} "
                 f"{'has' if complete == 1 else 'have'} a value in all of {allv}",
                 f"Só {complete:,} das {n_sel:,} observações{scope} "
                 f"{'tem' if complete == 1 else 'têm'} valor em todas as variáveis {allv}")
    if worst:
        parts = [L(f"{v} is missing in {k:,}", f"{v} falta em {k:,}") for k, v in worst[:4]]
        what += ": " + "; ".join(parts) + "."
    else:
        what += "."
    fix = L(f'see where the data are missing with "misstable summarize {" ".join(vars_)}"; leave out '
            "or fix the variables with many missing values.",
            f'veja onde faltam dados com "misstable summarize {" ".join(vars_)}"; retire ou corrija '
            "as variáveis com muitos valores faltantes.")
    if strings:
        fix = L(f"{join(strings)} {'is a string' if len(strings) == 1 else 'are string'} (text) "
                f'variable{"s" if len(strings) > 1 else ""}: convert with "destring" or "encode". ',
                f"{join(strings)} {'é variável' if len(strings) == 1 else 'são variáveis'} de texto: "
                'converta com "destring" ou "encode". ') + fix
    return Explanation(what, L(*_WHY_SAMPLE), fix, "missing")


def _outcome(ctx: Ctx) -> tuple[str, np.ndarray] | None:
    vars_ = _model_vars(ctx)
    if not vars_ or ctx.is_string(vars_[0]):
        return None
    rows = ctx.sample()
    idx = np.flatnonzero(rows)
    ok = ~np.any(np.column_stack([_missing(ctx, v, idx) for v in vars_]), axis=1)
    return vars_[0], np.asarray(ctx.ds.get(vars_[0]).data[idx][ok], dtype=float)


@rule(r"outcome does not vary.*")
def _outcome_constant(ctx: Ctx, _m) -> Explanation:
    out = _outcome(ctx)
    dep = out[0] if out else "y"
    what = L(f"The outcome {dep} has a single kind of value in the sample.",
             f"O desfecho {dep} tem um só tipo de valor na amostra.")
    if out is not None and out[1].size:
        y = out[1]
        kind = (L("is always 0 (negative)", "vale sempre 0 (negativo)") if np.all(y == 0) else
                L("is never 0 (always positive)", "nunca vale 0 (sempre positivo)"))
        what = L(f"In the sample used ({obs(y.size)}), {dep} {kind}; its values: {values_summary(y)}.",
                 f"Na amostra usada ({obs(y.size)}), {dep} {kind}; os valores: {values_summary(y)}.")
    return Explanation(
        what,
        L(f"{ctx.command or 'logit'} models a yes/no outcome: 0 means no, any other value means "
          "yes. To estimate anything, the sample needs both.",
          f"O {ctx.command or 'logit'} modela um desfecho sim/não: 0 quer dizer não, qualquer outro "
          "valor quer dizer sim. Para estimar, a amostra precisa ter os dois."),
        L(f'check with "tabulate {dep}"; for a count or a measurement, create a 0/1 variable, e.g. '
          f'"generate high = ({dep} > 10) if !missing({dep})".',
          f'confira com "tabulate {dep}"; para uma contagem ou medida, crie uma variável 0/1, por '
          f'exemplo "generate alto = ({dep} > 10) if !missing({dep})".'),
        "logit")


@rule(r"(\S+) is not a 0/1 variable")
def _not_binary(ctx: Ctx, m) -> Explanation:
    var = m.group(1)
    vals = ""
    if ctx.ds.has(var) and not ctx.is_string(var):
        x = np.asarray(ctx.ds.get(var).data, dtype=float)
        vals = L(f" Its values: {values_summary(x[x < SYSMISS])}.", f" Os valores: {values_summary(x[x < SYSMISS])}.")
    return Explanation(
        L(f"{var} should take only the values 0 and 1.{vals}",
          f"{var} deveria ter só os valores 0 e 1.{vals}"),
        L("This command needs a yes/no variable coded 0 (no) and 1 (yes).",
          "Este comando precisa de uma variável sim/não codificada como 0 (não) e 1 (sim)."),
        L(f'recode it, e.g. "generate {var}01 = ({var} == 2) if !missing({var})" when 2 means yes.',
          f'recodifique, por exemplo "generate {var}01 = ({var} == 2) if !missing({var})" quando 2 '
          "quer dizer sim."))


_ID_NAMES = re.compile(r"(^|_)(id|pid|hhid|cod\w*|code|country|pais|uf|state|estado|munic\w*|firm\w*|"
                       r"empresa|person|pessoa|escola|school)($|_)", re.I)
_TIME_NAMES = re.compile(r"^(year|ano|time|tempo|t|wave|onda|period|periodo|date|data|month|mes|"
                         r"quarter|trimestre|semana|week)\d*$", re.I)


def _guess(ctx: Ctx, rx: re.Pattern, default: str) -> str:
    return next((n for n in ctx.names if rx.search(n)), default)


@rule(r"panel variable not set.*")
def _xtset(ctx: Ctx, _m) -> Explanation:
    pid = _guess(ctx, _ID_NAMES, L("idvar", "varid"))
    tvar = _guess(ctx, _TIME_NAMES, L("timevar", "vartempo"))
    return Explanation(
        L(f"{ctx.command or 'This command'} works on panel data, and the data have not been "
          "declared as a panel.",
          f"O {ctx.command or 'comando'} trabalha com dados em painel, e a base ainda não foi "
          "declarada como painel."),
        L("Panel (xt) commands need to know which variable identifies each unit (person, firm, "
          "city) and, for most of them, which one identifies time.",
          "Os comandos de painel (xt) precisam saber qual variável identifica cada unidade "
          "(pessoa, empresa, cidade) e, na maioria, qual identifica o tempo."),
        L(f'declare it once with "xtset {pid} {tvar}" and run the command again.',
          f'declare uma vez com "xtset {pid} {tvar}" e rode o comando de novo.'),
        "xtset")


@rule(r"time variable not set.*")
def _tsset(ctx: Ctx, _m) -> Explanation:
    tvar = _guess(ctx, _TIME_NAMES, L("timevar", "vartempo"))
    return Explanation(
        L(f"{ctx.command or 'This command'} needs a time variable, and none has been declared.",
          f"O {ctx.command or 'comando'} precisa de uma variável de tempo, e nenhuma foi declarada."),
        L("Time-series commands and the operators L. (lag), F. (lead) and D. (difference) need to "
          "know the order of the observations in time.",
          "Os comandos de séries temporais e os operadores L. (defasagem), F. (avanço) e D. "
          "(diferença) precisam saber a ordem das observações no tempo."),
        L(f'declare it with "tsset {tvar}" (panel data: "xtset idvar {tvar}").',
          f'declare com "tsset {tvar}" (dados em painel: "xtset varid {tvar}").'),
        "xtset")


@rule(r"data not st")
def _not_st(ctx: Ctx, _m) -> Explanation:
    return Explanation(
        L("The st commands need the data declared as survival-time data first.",
          "Os comandos st precisam que os dados sejam declarados como de sobrevivência antes."),
        L("stset tells Stata which variable is the follow-up time and which marks the event "
          "(death, failure...); every st command reads that declaration.",
          "O stset diz ao Stata qual variável é o tempo de seguimento e qual marca o evento "
          "(óbito, falha...); todos os comandos st leem essa declaração."),
        L('run "stset timevar, failure(failvar)" (add id() if subjects have several records).',
          'rode "stset tempo, failure(evento)" (com id() se cada sujeito tiver vários registros).'),
        "stset")


@rule(r"last estimates not found")
def _no_estimates(ctx: Ctx, _m) -> Explanation:
    return Explanation(
        L(f"{ctx.command or 'This command'} works on the results of a model, and no model has been "
          "fitted yet (or the results were cleared).",
          f"O {ctx.command or 'comando'} trabalha sobre os resultados de um modelo, e nenhum modelo "
          "foi estimado ainda (ou os resultados foram apagados)."),
        L("predict, test, margins, estat and lincom use the last estimation command; commands such "
          "as clear or use with clear erase it.",
          "predict, test, margins, estat e lincom usam o último comando de estimação; comandos "
          "como clear ou use com clear o apagam."),
        L(f'fit the model first (for example "regress y x") and then run {ctx.command or "the command"}; '
          'saved models come back with "estimates restore name".',
          f'estime o modelo antes (por exemplo "regress y x") e depois rode o '
          f'{ctx.command or "comando"}; modelos guardados voltam com "estimates restore nome".'),
        "estimates")


@rule(r"conformability error", rc=503)
def _conformability(ctx: Ctx, _m) -> Explanation:
    sizes = []
    for name in dict.fromkeys(re.findall(r"\b[A-Za-z_]\w*\b", ctx.args)):
        mat = getattr(ctx.s, "matrices", {}).get(name)
        if mat is not None:
            r, c = mat.data.shape
            sizes.append(f"{name} {r}×{c}")
    what = L("The matrices do not have compatible sizes for this operation",
             "As matrizes não têm tamanhos compatíveis para esta operação")
    what += (": " + ", ".join(sizes) + ".") if sizes else "."
    return Explanation(
        what,
        L("A + B and A - B need the same number of rows and columns; A * B needs the columns of "
          "A to equal the rows of B; A' turns rows into columns.",
          "A + B e A - B precisam do mesmo número de linhas e colunas; A * B precisa que as "
          "colunas de A sejam iguais às linhas de B; A' transforma linhas em colunas."),
        L('check the sizes with "matrix list"; often a transpose (A\') fixes the product.',
          'confira os tamanhos com "matrix list"; muitas vezes uma transposta (A\') resolve o '
          "produto."),
        "matrix")


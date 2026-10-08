# Roadmap

Meta: cobrir os comandos do Stata 14, incluindo gráficos, com a interface o mais próxima possível da original.

As fases estão em ordem de dependência. Cada uma termina com algo que dá para abrir, testar e comparar com o Stata 14 (ver [docs/como-testar.md](docs/como-testar.md)). Os prazos não são estimados de antemão. Cada fase é dividida em pull requests revisáveis.

## Fase 0: fundação ✅ (este PR)

- [x] Lexer de do-files: comentários, `///`, `#delimit`
- [x] Macros locais e globais, expansão aninhada, `=exp`, incremento, funções estendidas básicas
- [x] Expressões: precedência, regras de missing, strings, primeira leva de funções
- [x] `display` com formatos e diretivas
- [x] Controle de fluxo: `forvalues`, `foreach`, `while`, `if`/`else`, `continue`
- [x] Prefixos `quietly`, `noisily` e `capture`; `do`, `run` e `include`
- [x] Códigos de erro `r(#)` e eco de comandos e blocos no padrão do Stata
- [x] Interface: Results, Command, Review, Variables, Properties, menus e barra de ferramentas
- [x] Console, modo batch e ferramenta de comparação com logs do Stata
- [x] Testes automáticos e integração contínua

## Fase 1: dados em memória

- [ ] Conjunto de dados com tipos `byte`, `int`, `long`, `float`, `double`, `str#` e `strL`, com a precisão de `float` reproduzida
- [ ] Leitura e gravação de `.dta`: formatos 114/115 (Stata 10–12), 117 (13) e 118 (14), com rótulos de valor, notas e características
- [ ] varlists: abreviações, curingas (`*`, `?`, `~`), intervalos (`a-c`), `_all`
- [ ] Qualificadores `if` e `in`, pesos (`fweight`, `aweight`, `pweight`, `iweight`) e o prefixo `by`/`bysort`
- [ ] `_n`, `_N`, `x[_n-1]` e expressões vetorizadas
- [ ] `use`, `save`, `clear`, `set obs`, `describe`, `list`, `generate`, `replace`, `drop`, `keep`, `rename`, `order`, `sort`, `gsort`, `label` (variable/define/values/list), `format`, `count`, `compress`, `recast`, `notes`, `browse`/`edit`
- [ ] `insheet`, `import delimited`, `import excel`, `export delimited`, `infile`, `input`
- [ ] Janelas Variables e Properties ligadas aos dados; Data Editor e Data Browser

## Fase 2: linguagem de programação completa

- [ ] `program define`/`end`, `syntax`, `args`, `gettoken`, `marksample`, `markout`
- [ ] `return`, `ereturn`, `sreturn`; `tempvar`, `tempname`, `tempfile`; `preserve`/`restore`
- [ ] Todas as funções estendidas de macro (`: variable label`, `: type`, `: dir`, `: list`…)
- [ ] `matrix` (definição, operadores, funções, `matrix list`)
- [ ] Busca de `.ado` no adopath (permite rodar pacotes do SSC); `which`, `findfile`
- [ ] `log using`/`log close`, `cmdlog`; Viewer com renderizador de SMCL e `help`
- [ ] `assert`, `set trace`, `creturn list` completo, `timer`

## Fase 3: estatística descritiva e manipulação

- [ ] `summarize` (`detail`), `tabulate` (uma e duas vias, `chi2`, `exact`, `row`/`col`), `tabstat`, `table`, `ci`, `ttest`, `prtest`, `correlate`, `pwcorr`, `centile`, `pctile`, `xtile`
- [ ] `egen` (todas as funções), `collapse`, `contract`, `reshape`, `merge`, `append`, `joinby`, `cross`, `expand`, `fillin`
- [ ] `duplicates`, `recode`, `encode`/`decode`, `destring`/`tostring`, `split`, `mvencode`/`mvdecode`, `isid`, `levelsof`, `distinct`
- [ ] Gerador de números aleatórios: `runiform()`, `rnormal()` etc. A sequência exata do Stata não é documentada, então os resultados serão estatisticamente equivalentes, mas não idênticos.

## Fase 4: Mata

- [ ] Interpretador de Mata: tipos, operadores matriciais, funções, controle de fluxo
- [ ] Integração com os dados: `st_data`, `st_view`, `st_store`, `st_local`, `st_numscalar`, `st_matrix`
- [ ] Biblioteca de funções Mata usada pelos comandos de estimação

## Fase 5: estimação

- [ ] Variáveis fatoriais (`i.`, `c.`, `#`, `##`) e operadores de séries temporais (`L.`, `F.`, `D.`, `S.`); `tsset`, `xtset`
- [ ] `regress`, `predict`, `test`, `testparm`, `lincom`, `nlcom`, `margins`, `estimates store/table/restore`
- [ ] `logit`, `logistic`, `probit`, `ologit`, `oprobit`, `mlogit`, `poisson`, `nbreg`, `glm`, `tobit`, `ivregress`, `xtreg`, `xtlogit`, `areg`
- [ ] Erros-padrão robustos, por cluster e bootstrap; motor `ml` (máxima verossimilhança)

## Fase 6: amostras complexas

- [ ] `svyset` (estratos, UPAs, pesos, FPC, múltiplos estágios), `svydescribe`
- [ ] Prefixo `svy:` com `mean`, `proportion`, `ratio`, `total`, `tabulate`, `regress`, `logit`, `logistic`, `poisson`
- [ ] Linearização de Taylor, BRR, jackknife, bootstrap e SDR; `subpop()`, `over()`, `estat effects`

## Fase 7: gráficos

- [ ] Motor de gráficos com sistema de *schemes* (`s2color`, `s1mono`…) e janela Graph
- [ ] `twoway` (`scatter`, `line`, `connected`, `lfit`, `qfit`, `lowess`, `area`, `bar`, `rcap`, `function`…), `graph bar`, `graph box`, `graph pie`, `graph dot`, `histogram`, `kdensity`, `graph matrix`
- [ ] Opções de título, eixos, legendas e marcadores; `by()`
- [ ] `graph combine`, `graph export` (PNG, PDF, SVG), `graph save`/`use` (formato próprio)
- [ ] Graph Editor (subconjunto)

## Fase 8: restante do Stata 14 e ferramentas da interface

- [ ] Sobrevivência (`stset`, `sts`, `stcox`, `streg`), séries temporais (`arima`, `var`), multinível (`mixed`, `melogit`), `sem`, `power`, `mi` (imputação múltipla) e demais comandos, priorizados pelo uso
- [ ] Do-file Editor com realce de sintaxe e execução de seleção; Variables Manager; diálogos dos menus; Preferences; impressão

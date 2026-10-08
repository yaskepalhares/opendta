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

## Fase 1: dados em memória e arquivos (em andamento)

**1a ✅**
- [x] Conjunto de dados com tipos `byte`, `int`, `long`, `float`, `double`, `str#` e `strL`, com a precisão de `float` reproduzida
- [x] varlists: abreviações, curingas (`*`, `?`, `~`), intervalos (`a-c`), `_all`
- [x] Qualificadores `if` e `in` e o prefixo `by`/`bysort` (pesos são reconhecidos na sintaxe; o uso chega com os comandos estatísticos)
- [x] `_n`, `_N`, `x[_n-1]`, `sum()` e expressões vetorizadas; `replace` sequencial com subscritos
- [x] `clear`, `set obs`, `describe`, `list`, `generate`, `replace`, `drop`, `keep`, `rename`, `order`, `sort`, `gsort`, `label` (variable/define/values/list/dir/drop/data), `format`, `count`, `compress`, `input`
- [x] Janelas Variables e Properties ligadas aos dados

**1b ✅**
- [x] Leitura e gravação próprias de `.dta` nos formatos 117 (13), 118 (14) e 119, com rótulos de valor, notas e características; verificação cruzada com a ReadStat nos testes
- [x] Leitura de `.dta` antigos (até 115, Stata 12 e anteriores) pela ReadStat (`pip install pyreadstat`)
- [x] `use`, `save`, `saveold`, `recast`, `notes`, `char`, `type`, `erase`, `browse`/`edit`
- [x] `insheet`, `outsheet`, `import delimited`, `export delimited`
- [x] File > Open/Save, arrastar e soltar `.dta`, Data Browser somente leitura

**1c ✅**
- [x] Formatos 113–115 (Stata 8–12): leitura própria; gravação 114/115 com `saveold, version(11|12)`
- [x] `import excel`/`export excel` (.xlsx; .xls com o pacote xlrd)
- [x] `infile` (formato livre e com dicionário), `infix`, `file open/write/read/close`
- [x] Data Editor com edição de células (gera `replace`/`set obs`/`generate`) e Variables Manager

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

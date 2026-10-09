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

## Fase 2: linguagem de programação completa ✅

- [x] `program define`/`end`, `syntax`, `args`, `gettoken`, `marksample`, `markout`
- [x] `return`, `ereturn`, `sreturn`; `tempvar`, `tempname`, `tempfile`; `preserve`/`restore`
- [x] Funções estendidas de macro (`: variable label`, `: type`, `: format`, `: label`, `: char`, `: dir`, `: list`, `: subinstr`, `: permname`, `: r(scalars)`…)
- [x] `matrix` (definição, operadores, funções, submatrizes, `matrix list`, nomes de linhas e colunas)
- [x] Busca de `.ado` no adopath (permite rodar pacotes do SSC); `sysdir`, `adopath`, `which`, `findfile`, `discard`
- [x] `log using`/`log close` (texto e SMCL), `cmdlog`; Viewer com renderizador de SMCL, `help`, `view`, `search`
- [x] `assert`, `set trace`, `creturn list`, `timer`
- [ ] Páginas de help para todos os comandos (hoje: as principais, escritas do zero)

## Fase 3: estatística descritiva e manipulação ✅

**3a ✅**
- [x] `summarize` (`detail`, `meanonly`, pesos), `tabstat`, `ci`, `ttest`, `prtest`, `correlate`, `pwcorr`, `centile`, `pctile`, `_pctile`, `xtile`
- [x] `tabulate` (uma e duas vias, `row`/`col`/`cell`/`expected`, `chi2`, `lrchi2`, `V`, `gamma`, `taub`, `exact` por enumeração), `tab1`, `tab2`, `table`

**3b ✅**
- [x] `egen`: funções de grupo (`count`, `mean`, `sd`, `total`, `min`, `max`, `median`, `mode`, `pctile`, `iqr`, `skew`, `kurt`, `mad`, `mdev`, `std`, `rank`), `group`, `tag`, `seq`, `fill`, `cut`, `concat`, `ends`, `diff` e funções de linha (`row*`, `any*`)
- [x] `collapse` (todas as estatísticas, pesos, `cw`), `contract`, `expand`, `fillin`
- [x] `merge` (1:1, m:1, 1:m, m:m, `_n`, `update`/`replace`, `keep()`, `assert()`, `keepusing()`), `append`, `joinby`, `cross`
- [x] `reshape long`/`wide` (com `@`, `string` e repetição da última especificação)

**3c ✅**
- [x] `duplicates` (report/examples/list/tag/drop), `isid`, `levelsof`
- [x] `recode`, `encode`/`decode`, `destring`/`tostring`, `split`, `mvencode`/`mvdecode`
- O `distinct` não é comando oficial: é um pacote do SSC e roda pelo adopath (fase 2).

**3d ✅**
- [x] Gerador Mersenne Twister de 64 bits (o mt64 do Stata 14), `set seed` (número ou estado), `c(seed)`, `set rng`
- [x] `runiform()`, `runiformint()`, `rnormal()`, `rbinomial()`, `rpoisson()`, `rchi2()`, `rt()`, `rbeta()`, `rgamma()`, `rexponential()`, `rlogistic()`, `rweibull()`, `rnbinomial()`, `rhypergeometric()`. A conversão da semente e os algoritmos de cada distribuição não são documentados: os resultados são estatisticamente equivalentes, mas não idênticos aos do Stata (casos em `compat/do/0307_aleatorios.do`).

**Conferência com o Stata 14**: os 28 casos de `compat/do/` têm referência gerada no Stata/SE 14.0; 25 batem linha a linha. As diferenças conhecidas estão em `compat/diferencas_conhecidas.txt`: o progresso do algoritmo de rede do `tabulate, exact` (tabelas maiores que 2×2) e, nos números aleatórios, `rnormal()`, `runiformint()`, `rpoisson()` e `rt()` (`runiform()`, `rbinomial()`, `rexponential()` e `rchi2(2)` já reproduzem a sequência exata do Stata).

## Fase 4: Mata (em andamento)

**4a–4c ✅ (primeira versão, aguardando conferência com o Stata)**
- [x] Lexer e analisador com a precedência do Mata; tipos real, string, complex, pointer e struct; escalares, vetores e matrizes
- [x] Operadores (aritméticos, com dois-pontos, transposta, Kronecker, junções `,` e `\`, sequências `..` e `::`, relacionais, lógicos, `?:`), subscritos `[i,j]` e `[|...|]`, missing propagados como no Stata
- [x] Controle de fluxo (`if`, `for`, `while`, `do`, `break`, `continue`, `return`), funções do usuário com tipos, argumentos opcionais e passagem por referência, `struct`, ponteiros
- [x] Blocos `mata`/`mata:` ... `end` com o eco do Stata, `mata: instrução`, `mata clear/describe/drop`; erros no formato do Mata (`<istmt>:  3499  x not found`)
- [x] Biblioteca: ~150 funções (matemáticas, matriciais, estatísticas, missing, strings, `printf`/`sprintf`)
- [x] Integração com os dados: `st_data`, `st_sdata`, `st_view` (cópia), `st_store`, `st_addvar`, `st_local`, `st_global`, `st_numscalar`, `st_matrix`, `stata()`

**4d**
- [x] Exibição, mensagens e erros conferidos com o Stata 14 (casos `compat/do/0401`–`0405`): números em `%12.0g`, colunas de largura única, cabeçalho do bloco, erros de compilação e de execução, linhas puladas no `mata:`
- [x] `st_view` que grava nos dados ao atribuir a elementos; arrays associativos (`asarray`); `optimize()` com avaliadores d0/d1/d2 e técnicas nr/bfgs
- [ ] Conferir `optimize()`, `st_view` e `asarray` (casos `0406` e `0407`)
- [ ] Classes; `moptimize()` (junto com o `ml` da fase 5)

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
- [ ] Do-file Editor com realce de sintaxe e execução de seleção; diálogos dos menus; impressão

## Fase 9 (última): idioma da interface

- [ ] Seletor de idioma em Preferences: inglês (padrão) ou português, valendo para menus, janelas, diálogos, explicações de erro e mensagens do OpenDTA
- [ ] Com o português ligado, as mensagens no formato do Stata podem continuar em inglês (para quem compara com logs e manuais) ou ser traduzidas, à escolha do usuário

## Melhorias em relação ao Stata (ao longo das fases)

- [x] Explicação dos erros: depois da mensagem e antes do r(#), uma linha diz o que deu errado no contexto do comando e outra sugere como resolver (nomes parecidos, arquivos parecidos na pasta, número de observações…). Em inglês; `set hints off` desliga
- [x] Expoente sobrescrito na notação científica exibida (`1.0000×10¹⁰`, `2.5×10⁻⁵`), mantendo a largura do formato; `set superscript off` (ou Help → Superscript Exponents) volta ao `e+10` do Stata. `string()`, macros, `file write`, exportação e `sprintf()`/`strofreal()` do Mata ficam sempre na forma do Stata

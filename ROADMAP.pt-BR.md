# Roadmap

[English](ROADMAP.md) | **Português**

Meta: cobrir os comandos do Stata 14, incluindo gráficos, com a interface o mais próxima possível da original.

**A versão 1.0** sai quando todos os comandos do Stata 14 estiverem cobertos. Os comandos que ainda faltam estão em [Faltam para a 1.0](#faltam-para-a-10); o que vem depois está em [Depois da 1.0](#depois-da-10).

As fases estão em ordem de dependência. Cada uma termina com algo que dá para abrir, testar e comparar com o Stata (ver [docs/testing.pt-BR.md](docs/testing.pt-BR.md)). Os prazos não são estimados de antemão. Cada fase é dividida em pull requests revisáveis.

[![Comandos implementados por área](https://yaskepalhares.github.io/opendta/roadmap-coverage.pt-BR.svg)](https://yaskepalhares.github.io/opendta/?lang=pt)

Um bloco por área de [Faltam para a 1.0](#faltam-para-a-10): o tamanho é o número de comandos e a cor, a parcela já implementada (opções de comandos existentes não entram na conta). Ele é redesenhado a cada merge no `main`. **[Abra a versão interativa](https://yaskepalhares.github.io/opendta/?lang=pt)** para ver os comandos de cada área e procurar um comando.

## Fase 0: fundação ✅

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

## Fase 1: dados em memória e arquivos ✅

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
- [x] Páginas de help para todos os comandos, escritas do zero

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
- [x] `sample` (percentual ou `count`, `if`/`in`, `by()` e o prefixo by)
- [x] `recode`, `encode`/`decode`, `destring`/`tostring`, `split`, `mvencode`/`mvdecode`
- O `distinct` não é comando oficial: é um pacote do SSC e roda pelo adopath (fase 2).

**3d ✅**
- [x] Gerador Mersenne Twister de 64 bits (mt64), `set seed` (número ou estado), `c(seed)`, `set rng`
- [x] `runiform()`, `runiformint()`, `rnormal()`, `rbinomial()`, `rpoisson()`, `rchi2()`, `rt()`, `rbeta()`, `rgamma()`, `rexponential()`, `rlogistic()`, `rweibull()`, `rnbinomial()`, `rhypergeometric()`. A conversão da semente e os algoritmos de cada distribuição não são documentados: os resultados são estatisticamente equivalentes, mas não idênticos aos do Stata (casos em `compat/do/0307_aleatorios.do`).

**Conferência com o Stata**: os 28 casos de `compat/do/` têm referência gerada no Stata; 25 batem linha a linha. As diferenças conhecidas estão em `compat/known_differences.txt`: o progresso do algoritmo de rede do `tabulate, exact` (tabelas maiores que 2×2) e, nos números aleatórios, `rnormal()`, `runiformint()`, `rpoisson()` e `rt()` (`runiform()`, `rbinomial()`, `rexponential()` e `rchi2(2)` já reproduzem a sequência exata do Stata).

## Fase 4: Mata ✅

**4a–4c ✅**
- [x] Lexer e analisador com a precedência do Mata; tipos real, string, complex, pointer e struct; escalares, vetores e matrizes
- [x] Operadores (aritméticos, com dois-pontos, transposta, Kronecker, junções `,` e `\`, sequências `..` e `::`, relacionais, lógicos, `?:`), subscritos `[i,j]` e `[|...|]`, missing propagados como no Stata
- [x] Controle de fluxo (`if`, `for`, `while`, `do`, `break`, `continue`, `return`), funções do usuário com tipos, argumentos opcionais e passagem por referência, `struct`, ponteiros
- [x] Blocos `mata`/`mata:` ... `end` com o eco do Stata, `mata: instrução`, `mata clear/describe/drop`; erros no formato do Mata (`<istmt>:  3499  x not found`)
- [x] Biblioteca: ~150 funções (matemáticas, matriciais, estatísticas, missing, strings, `printf`/`sprintf`)
- [x] Integração com os dados: `st_data`, `st_sdata`, `st_view` (cópia), `st_store`, `st_addvar`, `st_local`, `st_global`, `st_numscalar`, `st_matrix`, `stata()`

**4d ✅**
- [x] Exibição, mensagens e erros conferidos com o Stata (casos `compat/do/0401`–`0405`): números em `%12.0g`, colunas de largura única, cabeçalho do bloco, erros de compilação e de execução, linhas puladas no `mata:`
- [x] `st_view` que grava nos dados ao atribuir a elementos; arrays associativos (`asarray`); `optimize()` com avaliadores d0/d1/d2 e técnicas nr/bfgs
- [x] Conferir `optimize()`, `st_view` e `asarray` (casos `0406` e `0407`; no `optimize()` d0, só o ruído das derivadas numéricas difere)
- `moptimize()` passou para a fase 5, junto com o `ml`; classes do Mata, para a fase 8

## Fase 5: estimação ✅

- [x] Variáveis fatoriais (`i.`, `c.`, `#`, `##`) e operadores de séries temporais (`L.`, `F.`, `D.`, `S.`); `tsset`, `xtset`
- [x] `regress`, `predict`, `test`, `testparm`, `lincom`, `nlcom`, `margins`, `estimates store/table/restore`
- [x] `logit`, `logistic`, `probit`, `ologit`, `oprobit`, `mlogit`, `poisson`, `nbreg`, `glm`, `tobit`, `ivregress`, `xtreg`, `areg`
- [x] Erros-padrão robustos, por cluster e bootstrap (`vce(bootstrap)` e prefixo `bootstrap`); motor `ml` (máxima verossimilhança)
- [x] `xtlogit` e `xtpoisson` (fe e re), `qreg`, `lrtest`, `hausman`, `xtsum`, `xtdescribe`
- [x] `moptimize()` do Mata (avaliadores lf, d0–d2, lf0–lf2 e gf0–gf2; variância oim, opg e robusta; o caso `0408` conferido com o log do Stata 14)
- Conferência com os logs do Stata: casos `0501`–`0513` conferidos

## Fase 6: amostras complexas ✅

- [x] `svyset` (estratos, UPAs, pesos, FPC), `svydescribe`
- [x] Prefixo `svy:` com `mean`, `proportion`, `ratio`, `total`, `tabulate`, `regress`, `logit`, `logistic`, `poisson`; `subpop()` e `over()`
- [x] Linearização de Taylor; correção de Rao–Scott no `svy: tabulate` de duas vias
- [x] Vários estágios, BRR, jackknife, bootstrap e SDR; `estat effects`; `svy:` também com `ologit`, `oprobit`, `mlogit`, `nbreg`, `glm` e `tobit`
- Conferência com os logs do Stata: casos `0601`, `0602` e `0603` conferidos

## Fase 7: gráficos ✅

- [x] Motor de gráficos com esquemas próprios (`s2color`, `s1mono`… são aceitos e mapeados para eles) e janela Graph
- [x] `twoway` (`scatter`, `line`, `connected`, `lfit`, `qfit`, `lfitci`, `qfitci`, `area`, `bar`, `spike`, `dropline`, `rcap`, `rarea`, `rbar`, `rline`, `function`, `histogram`, `kdensity`, `scatteri`), `tsline`, `graph bar`/`hbar`, `graph box`/`hbox`, `graph pie`, `histogram`, `kdensity`, `graph matrix`
- [x] Opções de título, eixos, legendas e marcadores; `by()`
- [x] `graph combine`, `graph export` (PNG, PDF, SVG, EPS, PS, TIFF, JPEG), `graph save`/`use` (formato próprio), `graph dir`/`drop`/`rename`/`copy`/`display`/`describe`
- [x] `lowess`, `graph dot`, `marginsplot`, gráficos de diagnóstico (`rvfplot`, `rvpplot`, `qnorm`, `pnorm`)
- [x] Graph Editor simples (títulos, legenda, esquema, tamanho, grade, escalas logarítmicas), com desfazer
- Conferência com os logs do Stata: caso `0701` (só as mensagens de texto; o desenho é próprio)

## Fase 8: restante dos comandos do núcleo e ferramentas da interface ✅

- [x] Do-file Editor com abas, realce de sintaxe, localizar/substituir e execução de seleção (`doedit`); menus Data, Graphics e Statistics; impressão
- [x] Testes não paramétricos e de normalidade (`ranksum`, `signrank`, `kwallis`, `spearman`, `swilk`, `sktest`, `sdtest`) e comandos imediatos (`cii`, `ttesti`, `prtesti`, `tabi`)
- [x] `oneway`, `anova`, `tabulate, summarize()`, `estat` (`ic`, `vce`, `summarize`, `hettest`, `ovtest`, `classification`, `gof`, `bootstrap`)
- [x] `codebook`, `misstable`, `putexcel`
- [x] Classes do Mata (`class`, `extends`, métodos virtuais, `super`, `this`, membros estáticos, construtores; o caso `0409` conferido com o log do Stata 14)
- [x] Análise de sobrevivência: `stset` (com `id()`, `time0()`, `origin()`, `enter()`, `exit()`, `scale()`), `st`, `streset`, `stdescribe`, `stsum`, `sts list`/`generate`/`test`/`graph` (Kaplan-Meier, Nelson-Aalen, log-rank, Wilcoxon, Tarone-Ware, Peto-Peto), `stcox` (empates de Breslow e Efron, estratos, robusto), `estat phtest`, `streg` (exponencial, Weibull, Gompertz, lognormal, loglogística) e as opções de `predict`; os casos `0901`–`0902` conferidos com o log do Stata 14
- [x] Modelos multinível: `mixed` (ML e REML, níveis aninhados, inclinações aleatórias, covariâncias identity/independent/exchangeable/unstructured; `estat icc`, `group`, `recovariance`; `predict`) e `melogit`, `meprobit`, `mepoisson` (quadratura de Gauss-Hermite adaptativa); os casos `1001`–`1002` conferidos com o log do Stata 14
- [x] Séries temporais: `arima` (ML exata com filtro de Kalman, variância OPG), `var`, `varsoc`, `vargranger`, `dfuller` (valores críticos de Fuller, p-valores de MacKinnon), `corrgram`, `newey`; o caso `1101` conferido com o log do Stata 14
- [x] Imputação múltipla: `mi set` (mlong, flong), `register`, `describe`, `impute` (`regress`, `pmm`, `logit`, `ologit`, `mlogit`, `poisson`, `chained`), `estimate` (regras de Rubin, gl de Barnard-Rubin), `passive`, `xeq`, `extract`, `misstable`; caso `1201` (os valores imputados diferem dos do Stata: diferença conhecida)
- [x] `sem` (modelos de caminhos, AFC e variáveis latentes por ML, solução padronizada, `estat gof`, `predict`), `power` (`onemean`, `twomeans`, `pairedmeans`, `oneproportion`, `twoproportions`, `onecorrelation`, `twocorrelations`) e `irf` (`create`, `table`, `graph`, `describe`, `rename`, `drop` depois do `var`); os casos `1301`–`1303` conferidos com o log do Stata 14
- [x] Comandos completados: `anova` (somas de quadrados sequenciais, termos aninhados, termos de erro, medidas repetidas), `xtgee` e o estimador `pa` de `xtreg`/`xtlogit`/`xtpoisson`, `xtreg, be` e `mle`, `tabulate, summarize()` de duas entradas, `misstable tree`/`nested`, `ci`/`cii variances`, `prtest`/`prtesti` de duas amostras, `reshape i/j/xij/error`, a sintaxe antiga do `merge`, `sort ... in`, `list, sepby()`, `display _request()`, `export excel` para `.xls`, `set rng kiss32`; `replace` sequencial com operadores de séries temporais (`replace y = 2*L.y`); o caso `0805` conferido com o log do Stata 14
- [x] Diálogos dos menus (38 comandos de Statistics e Graphics), com prévia do comando e OK, Submit e Copy
- [x] Páginas de ajuda para todos os comandos
- [x] Pacote para macOS (`tools/build_macos.sh`: `OpenDTA.app` e `.dmg`; "Abrir com" para `.do` e `.dta`)
- Conferência com os logs do Stata: casos `0801`–`0804` conferidos

## Fase 9: idioma da interface ✅

- [x] Janela de configurações com barra lateral e busca: tema (claro, escuro ou do sistema), ícone do app, fonte com prévia, idioma, explicações de erro, expoentes sobrescritos, numerics e atualizações, aplicados na hora
- [x] Seletor de idioma em Settings: inglês (padrão) ou português, valendo para menus, janelas, diálogos, explicações de erro e mensagens do OpenDTA
- [x] Com o português ligado, as mensagens no formato do Stata podem continuar em inglês (para quem compara com logs e manuais) ou ser traduzidas, à escolha do usuário (`set language`, `set messages`)

## Fase 10: epidemiologia ✅

- [x] Tabelas para epidemiologistas: `cs`/`csi` (diferença e razão de riscos, frações atribuíveis ou prevenidas), `cc`/`cci` (razão de chances com limites exatos, de Cornfield, de Woolf ou baseados no teste; teste exato de Fisher), `ir`/`iri` (diferença e razão de taxas, limites exatos, mid-p), `mcc`/`mcci` (teste de McNemar e comparações pareadas), com estratificação `by()`: estimativas brutas e de Mantel-Haenszel e testes de homogeneidade
- [x] `tabodds` (chances por nível de exposição, testes de homogeneidade e de tendência, razões de chances ajustadas por Mantel-Haenszel, gráfico) e `mhodds` (razão de chances de Mantel-Haenszel, teste de escore para tendência)
- [x] `symmetry` (Bowker, Stuart-Maxwell, Bickenboller, significância exata), `dstdize` (padronização direta) e `istdize` (padronização indireta, SMR com limites exatos)
- [x] Concordância entre avaliadores: `kap` (kappa de Cohen, kappa ponderado; kappa de Fleiss para vários avaliadores) e `kappa` (formato de frequências)
- [x] Análise ROC: `roctab` (área com erros de DeLong ou de Hanley-McNeil, relatório detalhado, curva), `roccomp` (igualdade de áreas, correlacionadas ou independentes), `rocfit` (modelo binormal), `lroc` e `lsens` depois de `logit`, `logistic` e `probit`
- Conferido com o log do Stata 14: casos `1401` e `1402` (as diferenças que restam, e o motivo, estão em `compat/known_differences.txt`)

## Fase 11: manipulação de dados ✅

- [x] Variáveis e observações: `clonevar`, `insobs`, `range`, `separate`, `ipolate`, `expandcl`, `stack`, `xpose`, `drawnorm`, `corr2data`
- [x] Descrever e comparar dados: `ds`, `lookfor`, `inspect`, `labelbook`, `compare`, `cf`, `datasignature`, `snapshot`
- [x] Arquivos e sistema: `copy`, `mkdir`, `rmdir`, `dir`/`ls`, `shell` (e `!`), `filefilter`, `zipfile`, `unzipfile`, `checksum`, `outfile`, `unicode`
- [x] `statsby`, `putmata`; `sysuse` (dados de exemplo próprios do OpenDTA) e `webuse` (baixa sob demanda)
- [x] `icd9` e `icd10` com as listas de códigos de domínio público ICD-9-CM e ICD-10-CM do CMS americano (check, clean, generate, lookup, search); `icd9p` pelo formato
- [x] `odbc` pelo pacote opcional `pyodbc`; calendários de dias úteis (`bcal`, o formato `%tb`, `bofd()` e `dofb()`)
- Conferido com o log do Stata 14: caso `1501` (bate)

## Fase 12: epidemiologia, conclusão ✅

- [x] `cs`, `cc` e `ir` com `by()`: estimativas combinadas diretamente (`pool`) com o seu teste de homogeneidade, estimativas padronizadas interna e externamente (`istandard`, `estandard`) ou com pesos de uma variável (`standard()`), diferenças de riscos e de taxas padronizadas (`rd`, `ird`), registros de frequências binomiais (`binomial()`) e os testes de Breslow-Day e de Tarone do `cc` (`bd`, `tarone`)
- [x] `mhodds` com `by()`: a razão de chances dentro de cada grupo, a estimativa controlada pelos grupos e um teste de homogeneidade
- [x] `roctab` com o erro padrão de Bamber (`bamber`), limites binomiais exatos (`binomial`) e a curva de Lorenz com os índices de Gini e de Pietra (`lorenz`); `roccomp` com áreas binormais (`binormal`)
- [x] `rocreg`: estimação ROC não paramétrica com controles ajustados por covariáveis (estratificados ou por regressão linear, valores percentis empíricos ou normais) e os seus resumos (AUC, ROC(t), ROC inversa, AUC parcial), e o modelo paramétrico GLM probit com covariáveis da ROC; erros padrão por bootstrap com intervalos normal, percentil e com correção de viés
- Conferido com o log do Stata 14: caso `1601` (as diferenças que restam, e o motivo, estão em `compat/known_differences.txt`)

## Fase 13: manipulação de dados, opções ✅

- [x] `file` com leitura e gravação binária (`%1b`, `%2b`, `%4b`, suas formas com e sem sinal, `%4z`, `%8z`, `%#s`), `file seek`, `file set byteorder`, `file query`, `read write` e a função `byteorder()`
- [x] `infile` com `byvariable()`; `putexcel` para `.xls`
- [x] `import sas` (com `bcat()`), `import spss` e `import sasxport`, e um `export sasxport` básico, pelo pacote opcional `pyreadstat`
- [x] `describe` com `numbers`, `varlist`, `replace` e `using`; `list` com `header()`, `nodotz`, `mean`, `sum`, `N`, `labvar()` e `nocompress`
- Conferido com o log do Stata 14: caso `1701` (bate)

## Fase 14: programação ✅

- [x] `unab`, `tsunab`, `fvunab`, `unabcmd`, `fvexpand`, `fvrevar`, `tsrevar` e `_rmcoll`, com operadores de variáveis fatoriais e de séries temporais (inclusive as formas específicas como `i(2,3)bn.rep78`)
- [x] `mark`, `markin` e `markout` (com `sysmissok`); `marksample` e `syntax varlist(fv ts)` com variáveis fatoriais e operadores de séries temporais
- [x] `postfile`, `post`, `postclose` e `postutil`
- [x] `matlist` (especificações de colunas e linhas, bordas, formatos, títulos) e `tabdisp` (de uma a três variáveis, `by()`, totais, `concise`, `dotz`)
- [x] `program, sortpreserve`; `program nome, plugin` (explica que os plugins do Stata não podem ser carregados)
- [x] `pause`, `sleep`, `window` (no modo lote as caixas de diálogo agem como canceladas) e `viewsource` (mostra o código-fonte do próprio OpenDTA)
- [x] `signestimationsample`, `checkestimationsample`, `makecns` com matriz de restrições e `matcproc`
- Conferido com o log do Stata 14: caso `1801` (as diferenças que restam, e o motivo, estão em `compat/known_differences.txt`)

## Fase 15: análise de sobrevivência ✅

- [x] `strate`, `stptime`, `stci` (percentis, médias restrita e estendida), `stir`, `stmh`, `stmc` e `ltable`
- [x] `stsplit`, `stjoin`, `stfill`, `stgen`, `stvary`, `stbase` e `snapspan`
- [x] Dados de contagem no tempo: `ctset`, `ct`, `cttost`, `sttoct`; amostragem caso-controle com `sttocc`
- [x] `stcrreg` (riscos competitivos de Fine-Gray) e `stcurve` depois de `stcox`, `streg` e `stcrreg`
- [x] `streg, distribution(ggamma)` e `frailty(gamma|invgaussian)`, com `shared()`; os empates `exactm` e `exactp` do `stcox`
- Conferido com o log do Stata: caso `1901`

## Fase 16: modelos multinível ✅

- [x] `meglm` (famílias gaussian, bernoulli, binomial, gamma, nbinomial, ordinal e poisson com suas ligações), `mecloglog`, `menbreg` (dispersão mean e constant), `meologit`, `meoprobit`, `meqrlogit` e `meqrpoisson`; `predict` e `estat icc`/`recovariance`/`group` depois deles
- [x] `mestreg` (exponencial, Weibull, lognormal, loglogística e gama com um efeito aleatório por grupo)
- [x] Efeitos aleatórios cruzados no `mixed` (`|| _all: R.fator`)
- Conferido com o log do Stata 14: caso `1902` (as diferenças que restam, e o motivo, estão em `compat/known_differences.txt`)

## Fase 17: estatística básica ✅

- [x] Prefixos de reamostragem `jackknife` (observações ou conglomerados, `mse`, `dof()`), `permute` e `simulate`
- [x] `ml` com os avaliadores `lf`, `d0`, `d1` e `d2` (`ml model`, `init`, `search`, `maximize`, `display`, `clear`, `score`; variâncias `oim`, `opg`, `robust` e `cluster()`; derivadas numéricas), `mleval`, `mlsum`, `mlvecsum`, `mlmatsum` e `mlexp`
- [x] Construção de modelos: `mkspline` (linear, `pctile`, `cubic`, `marginal`), `xi`, `nestreg`, `stepwise` e `fp`/`mfp`
- [x] Depois da estimação: `testnl`, `predictnl`, `linktest`, `contrast`, `pwcompare`, `pwmean` e `suest`; `margins` depois de `tobit`, `ologit` e `oprobit`
- [x] Transformações: `boxcox` (modelo do lado esquerdo), `ladder`, `gladder` e `lnskew0`
- [x] Testes e correlações: `bitest`/`bitesti`, `ksmirnov`, `median`, `nptrend`, `runtest`, `tetrachoric`, `ktau` e `pcorr`
- [x] Concordância e tamanhos de efeito: `icc`, `loneway`, `esize`/`esizei` e `ameans`
- [x] Exibições e gráficos: `cumul`, `stem`, `lv`, `smooth`, `quantile`, `symplot`, `serrbar`, `dotplot`, `sunflower`, `lpoly`, `cusum`
- Conferido com o log do Stata 14: o caso `1903` confere, exceto pelas diferenças de `compat/known_differences.txt` (passos iniciais de `ml`/`ologit`, a busca do `lnskew0`, o pseudosigma de profundidade 1 do `lv`, os sorteios de `permute`/`simulate` e a tabela de caixas do `sunflower`)

## Fase 18: séries temporais ✅

- [x] Ferramentas de dados: `tsappend`, `tsfill`, `tsreport` e o prefixo `rolling`
- [x] Filtros e suavizadores: `tsfilter` (`hp`, `bk`, `cf`, `bw`) e `tssmooth` (`ma`, `exponential`, `dexponential`, `hwinters`, `shwinters`, `nl`)
- [x] Análise espectral e testes de ruído branco: `pergram`, `cumsp`, `xcorr`, `wntestb`, `wntestq` e `psdensity` depois de `arima` e `arfima`
- [x] Raízes unitárias e regressão com erros AR(1): `dfgls`, `pperron` e `prais` (Prais-Winsten, Cochrane-Orcutt, `rhotype()`, `ssesearch`, `twostep`)
- [x] Pós-estimação de VAR: `varlmar`, `varnorm`, `varstable`, `varwle`, `varbasic` e `fcast` (depois de `var`, `svar` e `vec`)
- [x] Cointegração e VAR estrutural: `vecrank`, `vec` e `svar` (matrizes de curto prazo `aeq()`/`beq()`); `irf` depois de `svar` e `vec`, com erros-padrão bootstrap (`bs`, `bsp`), e `irf ctable`, `cgraph`, `ograph` e `add`
- [x] Variância condicional e memória longa: `arch` (ARCH, GARCH e os demais termos de variância, ARCH-M, erros ARMA, distribuições t e GED, `het()`), `arfima` (ML exata) e `mgarch` (`dvech`, `ccc` e `dcc`)
- [x] `arima` com `sarima()`, `mar()`, `mma()`, `condition` e `diffuse`
- [x] Modelos de espaço de estados pelo filtro de Kalman com estados iniciais difusos: `ucm`, `sspace` e `dfactor`
- Conferido com o log do Stata 14: o caso `1904` confere, exceto pelas diferenças em `compat/known_differences.txt` (os valores críticos do `dfgls` por defasagem, os logs de iteração de `prais, ssesearch` e `svar`, os erros padrão do `fcast` a partir do quarto passo, as tabelas bootstrap do `irf`, e as estimativas do `ucm` com `cycle()` e do `mgarch dcc`, em que o Stata chega a outro ponto)

## Fase 19: dados em painel ✅

- [x] Modelos de efeitos aleatórios por quadratura de Gauss-Hermite (adaptativa `mvaghermite` e `ghermite`, `intpoints()`): `xtprobit`, `xtcloglog`, `xtologit`, `xtoprobit`, `xttobit`, `xtintreg` e `xtnbreg` (efeitos aleatórios, efeitos fixos condicionais e `pa`), com o teste LR do efeito aleatório e `vce(robust)`/`vce(cluster)`; `xtlogit, re` passou ao mesmo motor; `quadchk`
- [x] Estimadores lineares de painel: `xtgls` (`panels()`, `corr()`, `igls`), `xtpcse`, `xtregar` (`fe`, `re`, `rhotype()`, `rhof()`, `twostep`), `xtrc`, `xthtaylor` (`endog()`, `constant()`, `varying()`, `amacurdy`) e `xtivreg` (efeitos aleatórios G2SLS e EC2SLS, `be`, `fe`, `fd`)
- [x] `xtfrontier`: ineficiência invariante no tempo (`ti`) e com decaimento variável (`tvd`), fronteiras de produção e de custo (`cost`)
- [x] `xtunitroot`: `ht`, `breitung`, `hadri` e `fisher`
- [x] Painel dinâmico por GMM: `xtabond`, `xtdpdsys` e `xtdpd` (um e dois passos, erros robustos com a correção de Windmeijer, `pre()`, `endogenous()`, `inst()`, `maxldep()`), com `estat abond` e `estat sargan`
- [x] Descrição: `xtline` (com `overlay`), `xttab` e `xtdata` (`fe`, `be`)
- Conferido com o log do Stata 14: o caso `1905` bate, exceto pelas diferenças de `compat/known_differences.txt` (os caminhos de iteração dos maximizadores numéricos e do `igls`, o `xtgls` com painéis correlacionados e o `igls` com `psar1`, os componentes de variância do `xtregar, re`, o teste F de `u_i` do `xtivreg, fe`, o R² entre do `xtivreg, fd`, a verossimilhança plana do `xtfrontier, tvd` e a estatística de Breitung); `xtunitroot llc` e `ips`, que o Stata 14 executa, ainda não foram implementados

## Fase 20: estatística multivariada ✅

- [x] Componentes principais e análise fatorial: `pca`, `factor` (`pf`, `pcf`, `ipf`, `ml`), `rotate` (`varimax`, `quartimax`, `equamax`, `parsimax`, `promax`, `oblimin`, `quartimin`), `alpha`, `screeplot`, `loadingplot` e `scoreplot`, com `estat` (`kmo`, `smc`, `anti`, `loadings`, `correlation`, `residuals`, `summarize`, `structure`, `common`) e `predict`
- [x] Testes e modelos multivariados: `hotelling`, `mvreg`, `manova` (Wilks, Pillai, Lawley-Hotelling e Roy), `canon` e `mvtest` (`means`, `covariances`, `correlations`, `normality`)
- [x] Análise de agrupamentos: `cluster` (`kmeans`, `kmedians`, as sete ligações hierárquicas, `generate`, `dendrogram`, `stop` com `calinski` e `duda`, `list`, `dir`, `drop`, `rename`, `notes`) e `matrix dissimilarity`
- [x] Análise discriminante: `discrim` (`lda`, `qda`, `knn`, `logistic`) e `candisc`, com `estat` (`classtable`, `errorrate`, `list`, `grdistances`, `covariance`) e `predict`
- [x] Análise de correspondência e escalonamento: `ca`, `mca`, `mds` e `mdsmat` (clássico e moderno), `procrustes` ortogonal e `biplot`, com `estat` e `predict` onde se aplicam
- Conferido com o log do Stata 14: o caso `1906` confere, exceto pelas diferenças em `compat/known_differences.txt` (os caminhos de iteração de `factor, ml`, `discrim logistic` e `mds, method(modern)`, e os números do `factor, ml` que decorrem de onde o `ml` do Stata para)

## Fase 21: modelos, parte 1 (binários e contagem, censurados e seleção, variáveis instrumentais) ✅

- [x] Resposta binária e fracionária: `cloglog`, `scobit`, `hetprobit` (`het()`), `fracreg` (`logit`, `probit`, `cloglog`), `binreg` (`or`, `rr`, `hr`, `rd`) e `glogit`
- [x] Dados de contagem: `zip` e `zinb` (`inflate()`, `probit`, `vuong`, `irr`), `tpoisson`, `tnbreg` e `gnbreg` (`lnalpha()`)
- [x] Modelos censurados, truncados e de seleção: `intreg`, `truncreg` (`ll()`, `ul()`), `heckman` (máxima verossimilhança e `twostep`), `heckprobit`, `heckoprobit` e `biprobit` (lista de variáveis compartilhada e forma `(y1 = x) (y2 = z)`), com a função de distribuição normal bivariada em `core/bvnorm.py` (identidade de Plackett, derivadas exatas)
- [x] Variáveis instrumentais: `ivregress liml` (classe k) e `gmm` (dois passos, `wmatrix()`), `ivprobit` e `ivtobit` (máxima verossimilhança, teste de Wald de exogeneidade), `ivpoisson gmm` (momentos aditivos e `multiplicative`, `onestep`) e `cfunction`
- [x] Um motor de máxima verossimilhança de várias equações (`estimation/mleq.py`) com as verossimilhanças separadas dos comandos (`estimation/likelihoods.py`), reaproveitado pelos modelos de escolha da próxima parte
- Conferido com o log do Stata 14: o caso `1907` confere, exceto pelas diferenças em `compat/known_differences.txt` (os caminhos de iteração e o ponto de parada de cada otimizador; blocos cujos dados levam um parâmetro à fronteira, de modo que as estimativas do Stata dependem do caminho dele: `scobit`, `zinb`, `tnbreg`, `gnbreg`; `binreg, hr`, em que o IRLS do Stata não converge com os dados do caso; `glogit`, rodado com contagens agrupadas diferentes das do Stata; a variância robusta de `ivregress liml` e de `ivpoisson cfunction`). Achados e corrigidos no caminho: `binreg` é IRLS com a variância da informação esperada, `ivpoisson cfunction` usa o resíduo multiplicativo, `ivregress gmm` com matriz de pesos robusta mostra a variância sanduíche, e os modelos só com constante de `tpoisson` e `gnbreg` seguem os pontos de referência do Stata
- Ainda não conferidos com o Stata: `glogit` (dados do caso), `fracreg cloglog`, `heckoprobit` e as opções que o caso `1907` não exercita

## Fase 22: modelos, parte 2 (escolha, exatos, sistemas e efeitos de tratamento) ✅

- [x] Sistemas e restrições: `constraint`, `cnsreg`, `sureg` (com `isure`), `reg3` (`2sls`, `3sls`, `ols`, `sure`), `eivreg`, `vwls` e `rreg` (Huber e depois biquadrado), com a álgebra em `core/syseq.py`
- [x] Modelos de escolha: `clogit` (qualquer número de alternativas escolhidas por grupo), `rologit` (`exactp`, `breslow`, `efron` e `none` para empates), `asclogit`, `slogit`, `mprobit`, `asmprobit` e `asroprobit` (simulador GHK com sorteios de Hammersley), com as verossimilhanças em `core/choice.py`
- [x] Métodos exatos: `exlogistic` e `expoisson` (máxima verossimilhança condicional, estimativas medianas não viesadas, intervalos de confiança e valores de p exatos, a partir da distribuição exata das estatísticas suficientes em `core/exact.py`)
- [x] Modelos não lineares: `nl`, `nlsur` (NLS, FGNLS e IFGNLS), `gmm` (um e dois passos, matrizes de pesos robusta e por conglomerado) com expressões substituíveis, e `frontier` (`hnormal`, `exponential`, `tnormal`)
- [x] Efeitos de tratamento: `teffects ra`, `ipw`, `ipwra`, `aipw`, `nnmatch` e `psmatch` (variância de Abadie-Imbens), `eteffects` (função de controle), `stteffects ra`, `ipw` e `ipwra`, `etregress` (máxima verossimilhança e `twostep`) e `etpoisson` (Gauss-Hermite), com os estimadores em `core/treatfx.py`
- Conferido com o log do Stata 14: caso `1908` (as diferenças e o motivo de cada uma estão em `compat/known_differences.txt`)
- Ainda não conferidos com o Stata: o layout das tabelas de todos os comandos desta parte (marcados com `VERIFY` no código), a variância de `eivreg`, os pesos de `stteffects` e as opções que o caso não exercita

## Fase 23: outras suítes (sem, gsem, irt, estimação bayesiana, power, mi impute e svy) ✅

- [x] `sem` com `method(mlmv)` (verossimilhança normal por padrão de ausência, valores iniciais por EM em `core/mvnmiss.py`) e `method(adf)`, `group()` e `ginvariant()`, rótulos de caminho (`@nome`) e os subcomandos de `estat` `eqgof`, `eqtest`, `framework`, `ggof`, `ginvariant`, `mindices`, `residuals`, `scoretests`, `stable`, `stdize:`, `summarize` e `teffects`
- [x] `gsem` com as famílias `gaussian`, `bernoulli`, `binomial`, `poisson`, `nbinomial`, `gamma` e `ordinal`, variáveis latentes por observação e por conglomerado (`M[id]`), quadratura de Gauss-Hermite adaptativa e `predict`
- [x] `irt` (`1pl`, `2pl`, `3pl`, `grm`, `nrm`, `pcm`, `gpcm`, `rsm` e `hybrid`) sobre o `gsem`, com `estat report` e `irtgraph` (`icc`, `iif`, `tcc`, `tif`)
- [x] `bayesmh` (Metropolis-Hastings de passeio aleatório adaptativo em `core/mcmc.py`, oito verossimilhanças e dez prioris), o prefixo `bayes` para `regress`, `logit`, `probit` e `poisson`, `bayesstats` (`summary`, `ess`, `ic`), `bayestest` (`interval`, `model`) e `bayesgraph`
- [x] `power` com `onevariance`, `twovariances`, `pairedproportions` e `oneway`, os testes `lrchi2` e `fisher` do `twoproportions` e a opção `graph`
- [x] `mi impute mvn` (aumento de dados), `monotone`, `truncreg`, `intreg` e `nbreg`
- [x] `svy:` com `cloglog`, `intreg`, `truncreg`, `heckman`, `heckprobit`, `zip`, `zinb`, `gnbreg`, `hetprobit`, `scobit`, `biprobit`, `fracreg`, `slogit`, `mprobit`, `clogit` e `ivregress`
- Conferido com o log do Stata 14: caso `1909` (as diferenças e o motivo de cada uma estão em `compat/known_differences.txt`). Achados e corrigidos no caminho: a variância do `sem, method(adf)` é a Hessiana numérica da discrepância, não a aproximação de Gauss-Newton; o `ivregress` ignorava os pesos; os intervalos de confiança de `sigma`, `rho`, `alpha` e `lambda` no `svy:` usam os graus de liberdade do desenho; o `irt 3pl` tem um pseudoacerto comum, a menos que se dê `sepguessing`. Os sorteios do `bayesmh` e do `mi impute` vêm do gerador do próprio OpenDTA e não reproduzem os do Stata

## Fase 24: programação, parte 2 ✅

- [x] `makecns` com lista de números de restrições (restrições que não podem ser lidas ficam de fora com uma nota), as opções `displaycns` e `nocnsnotes`, e `matrix dispCns`
- [x] `constraint #` sem a palavra `define`
- [x] `matlist` com nomes de equações: `showcoleq()` (`first`, `each`, `combined`, `lcombined`, `rcombined`), `roweqonly`, `coleqonly`, `keepcoleq`, `colorcoleq()` e `lines(eq)`; `matrix rownames` e `colnames` definem a equação dos nomes escritos `eq:nome`
- [x] `window` na interface gráfica: caixas de diálogo de arquivo `fopen` e `fsave`, caixas de mensagem `stopbox`, `push` para a janela Review, o menu User montado com `window menu`, e `manage` (`minimize`, `restore`, `maintitle`, `forward`, `close`)
- Conferido com o log do Stata 14: caso `1910`. Achados e corrigidos no caminho: o `makecns` pula números sem restrição, mantém as restrições redundantes e não mostra nada com `r`; o `matrix dispCns` exige resultados gravados por um comando de estimação; o `matlist` mostra as equações das colunas como faixas, uma linha de rótulo para cada equação das linhas e, nos painéis seguintes, deixa de fora a metade superior vazia de uma matriz simétrica

## Fase 25: manipulação de dados, parte 2 ✅

- [x] `table` com variável de supercoluna (`table rowvar colvar supercolvar`), a opção `scol` e `by()` com até quatro variáveis
- [x] `list` com `string()` e o formato `display` (escolhido também quando a tabela não cabe na linha)
- [x] `import sasxport` com `describe` e `member()` (leitor próprio do OpenDTA de todos os membros de um arquivo XPORT, `io/xport.py`) e os rótulos de valores do `formats.xpf`; `export sasxport` com `rename` e `vallabfile()` (`xpf`, `sascode`, `both`, `none`)
- Conferido com o log do Stata 14: caso `1911`. Achados e corrigidos no caminho: o `table` mantém todas as linhas e todas as combinações das variáveis do `by()`, com as supercolunas separadas por 4 espaços; o `list, string()` termina com `..` os valores cortados; o `import sasxport, describe` tem o leiaute e os resultados guardados do Stata

## Fase 26: epidemiologia, parte 2 ✅

- [x] `rocreg` com `probit ml`: o modelo ROC binormal por máxima verossimilhança, com covariáveis da ROC (`roccov()`), covariáveis de controle (`ctrlcov()`, linear) e a AUC pelo método delta
- [x] `rocreg` com `cluster()`: erros-padrão robustos a clusters depois de `ml` e reamostragem de clusters inteiros no bootstrap
- Conferido com o log do Stata 14: caso `1912` (as diferenças e o motivo de cada uma estão em `compat/known_differences.txt`). Achados e corrigidos no caminho: o `probit ml` ajusta, para cada classificador, o modelo de ajuste por covariáveis dos controles (médias e desvios-padrão de casos e controles) e tira dele os parâmetros da ROC e a AUC pelo método delta; vários classificadores entram num só modelo, com o teste de AUCs iguais

## Fase 27: estatística básica, parte 2 ✅

- [x] `ml` com os métodos `lf0`, `lf1`, `lf2`, `gf0`, `gf1` e `gf2`, `constraints()` (maximização no subespaço restrito, `e(Cns)`) e `svy` (pesos do desenho e variância linearizada)
- [x] `boxcox` com `model(rhsonly)`, `model(lambda)` e `model(theta)` e `notrans()`
- [x] `icc` com `testvalue()` e dados não balanceados no modelo de um fator; `loneway` com `mean`, `median` e `exact`
- [x] `esize` e `esizei` com `unequal` e `welch`; `pwmean` com `unequal` (graus de liberdade de Satterthwaite por comparação); `ameans` com `fweight` e `aweight`
- [x] `lv` com `generate` e `tail()`; `stem` com `width()` (e uma correção para dados decimais guardados como `float`)
- Conferido com o log do Stata 14: caso `1913` (as diferenças e o motivo de cada uma estão em `compat/known_differences.txt`). Achados e corrigidos no caminho: o `ml` com restrições grava e(Cns) e os resultados de survey do Stata; o `boxcox` itera por Newton–Raphson a partir de 1 e transforma todas as variáveis da varlist; o `icc` deixa de fora os alvos com avaliações faltantes; o `loneway` guarda o r() do Stata; o `esize, unequal/welch` informa os graus de liberdade

## Fase 28: dados em painel, parte 2 ✅

- [x] `xtivreg` com `first` (regressões de primeiro estágio com o mesmo estimador) e `regress` (todos os regressores exógenos)
- [x] `xtdata` com `re` e `ratio()` (a transformação de MQG de efeitos aleatórios, com uma variável `constant`)
- [x] `xtrc` com `betas` (melhores preditores lineares dos coeficientes de cada painel, `e(beta_ps)` e `e(V_ps)`)
- [x] `xtunitroot fisher` com `drift`, `demean` e `pperron`; `xtunitroot hadri` com `kernel()` (variâncias de longo prazo de Bartlett, Parzen e espectral quadrática com um número dado de defasagens)
- [x] `xtabond`, `xtdpdsys` e `xtdpd` com `artests()`; `xtdpd` com `hascons`
- Conferido com o log do Stata 14: caso `1914` (as diferenças e o motivo de cada uma estão em `compat/known_differences.txt`). Achados e corrigidos no caminho: as regressões de primeiro estágio do `xtivreg, first` seguem o estimador (o G2SLS tem cabeçalho próprio); o `xtdata, re` mantém a variável de painel e acrescenta a constante; o `xtunitroot fisher` e o `hadri` seguem o layout e o r() do Stata; o `xtabond` descarta instrumentos duplicados e a constante do `hascons` com uma nota

## Fase 29: séries temporais, parte 2 ✅

- [x] `prais` com `vce(robust)`, `vce(hc2)`, `vce(hc3)` e `vce(cluster)` (variância semirrobusta da regressão transformada)
- [x] `arch` com `constraints()` e `from()`; `arfima` com `vce(robust)`
- [x] `svar` com `aconstraints()` e `bconstraints()`; `sspace` e `dfactor` com `constraints()` nos coeficientes
- [x] `predict` depois de `ucm`, `sspace` e `dfactor` (previsões um passo à frente, resíduos e resíduos padronizados); a subopção `frequency()` do `cycle()` do `ucm`; `psdensity` depois de `ucm`
- Conferido com o log do Stata 14: caso `1915` (as diferenças e o motivo de cada uma estão em `compat/known_differences.txt`). Achados e corrigidos no caminho: a log-verossimilhança do `prais` deixa de fora o termo do jacobiano; o `arfima` refina as estimativas numa segunda maximização; o `svar` grava todas as entradas de A e B em e(b), com variância zero para as fixas, e o e(Cns); o `ucm` confere a frequência do ciclo e o teste de Wald cobre a frequência e o amortecimento; o `sspace` e o `dfactor` padronizam os resíduos em conjunto

## Fase 30: modelos multinível, parte 2 ✅

- [x] `meglm` e os outros comandos `me` com níveis aninhados de efeitos aleatórios (`|| school: || class:`), por quadratura adaptativa em produto sobre os níveis
- [x] Pesos nos comandos `me` e no `mestreg`: `fweight`, `pweight` e `iweight` no nível da observação e as opções `fweight()` e `pweight()` de cada equação de efeitos aleatórios
- [x] `vce(robust)` e `vce(cluster clustvar)` (o padrão com pesos amostrais): variância sanduíche com os escores dos grupos do nível mais alto
- [x] `predict` depois do `mestreg` com `median`, `mean` e `hazard` em todas as distribuições
- Conferido com o log do Stata 14: caso `1916` (as diferenças e o motivo de cada uma estão em `compat/known_differences.txt`). Achados e corrigidos no caminho: com níveis aninhados os pontos de quadratura de um nível inferior acompanham os efeitos do nível acima e as derivadas seguem os pontos, o que dá as estimativas e os erros-padrão do Stata; uma variância que vai a zero não infla mais a verossimilhança; o `meqrlogit` e o `meqrpoisson` aceitam níveis aninhados; os pesos de frequência entram na tabela de grupos; a mediana do `mestreg` usa só os efeitos fixos

## Fase 31: estatística multivariada, parte 2 ✅

- [x] `clustermat` (agrupamento hierárquico de uma matriz de dissimilaridades, com `shape()`, `add`, `clear` e `labelvar()`); as medidas binárias de similaridade do `cluster` e do `matrix dissimilarity` (`matching`, `Jaccard`, `Russell`, `Hamann`, `Dice`, `antiDice`, `Sneath`, `Rogers`, `Ochiai`, `Yule`, `Anderberg`, `Kulczynski`, `Pearson`, `Gower2`); `cluster kmeans` e `kmedians` com `keepcenters`
- [x] `manovatest` (testes de termos depois do `manova`, com termos de erro depois de `/` e `ytransform()`)
- [x] `camat`, `cabiplot`, `caprojection`, `mcaplot` e `mcaprojection`
- [x] `mdslong`, `mdsconfig` e `mdsshepard`; `procoverlay`; `procrustes` com `transform(unrestricted)` e `transform(oblique)`
- Conferido com o log do Stata 14: o caso `1917` bate, exceto o histórico de iterações de `mdslong, method(modern)` (o otimizador do Stata não é documentado; ver `compat/known_differences.txt`). Encontrado e corrigido no caminho: as variáveis `_ord` e `_hgt` de uma análise hierárquica são por posição no dendrograma (o id da k-ésima folha e a altura da junção com a seguinte), o grupo com a menor observação vem primeiro, as medidas binárias de similaridade guardam as alturas como similaridades e `cluster generate` recusa um número de grupos que os empates tornam ambíguo; `clustermat` informa as observações novas e o nome da análise e põe a variável de rótulos por último; `manovatest` com vários termos testa-os juntos na codificação por nível de base, para quando o termo de erro tem menos graus de liberdade que as variáveis dependentes e tem `showorder`; `list` abrevia um nome de modo que ele corresponda só àquela variável (`_clus~id`)

## Fase 32: modelos, parte 3 ✅

- [x] `tpoisson` e `tnbreg` com `ul()` e limites dados por variáveis; `tnbreg` e `nbreg` com `dispersion(constant)` (o `nbreg` aceitava a opção e ajustava o modelo de dispersão média)
- [x] `glogit` com `or`; `glm` e `binreg` com `mu()`
- [x] `offset()` no `intreg`, `truncreg`, `heckman`, `heckprobit` e `heckoprobit` (antes aceito e ignorado), `offset1()` e `offset2()` no `biprobit`; a linha do offset na tabela de coeficientes desses modelos e dos modelos de contagem; `heckman` com `mills()` (duas etapas e máxima verossimilhança)
- [x] `sureg` e `reg3` com `constraints()`
- [x] `ivregress gmm` com `wmatrix(hac kernel #)`, `center` e `igmm`
- [x] `frontier` com `uhet()`, `vhet()` e `cm()`
- Conferido com o log do Stata 14: caso `1918` (as diferenças e o motivo de cada uma estão em `compat/known_differences.txt`). O Stata 14 não tem `ul()` no `tpoisson` e no `tnbreg` (veio depois): o OpenDTA mantém a opção como extensão. Encontrado e corrigido no caminho: as verossimilhanças binomiais negativas somam ln Gama(y + r) - ln Gama(r) termo a termo quando a dispersão vai a zero (a diferença das duas perdia dígitos e punha a verossimilhança NB1 acima da Poisson); a probabilidade da janela de uma contagem com limite superior é somada termo a termo; `tnbreg` e `nbreg` com `dispersion(constant)` começam o modelo completo do ajuste só com constante e o `tnbreg` guarda o e() do Stata; `sureg` e `reg3` com restrições tiram a primeira covariância dos resíduos do ajuste restrito equação por equação e guardam o e() do Stata (`Sigma`, `Cns`, as estatísticas de cada equação), e o `2sls` restrito usa estatísticas z; `ivregress gmm` com matriz de pesos HAC: a largura de banda do núcleo espectral quadrático é # + 1, a nota do VCE HAC, e(J), e(hac_lag), e(S) e e(W), e o histórico do `igmm`; `frontier` com `uhet()` ou `vhet()` mostra as duas variâncias como equações, e com `cm()` os parâmetros `/lnsigma2` e `/ilgtgamma`

## Fase 33: poder e tamanho de amostra, parte 2 ✅

- [x] `power twoway` (F não central, efeitos de linha, coluna e interação, `cellweights()`, `showcellsizes`), `power repeated` (efeitos entre, intra e entre-intra, simetria composta ou `covmatrix()`, correção de Greenhouse-Geisser), `power cmh`, `power mcc` (com `m()` e `compare`), `power trend`, `power cox`, `power exponential` e `power logrank`; tamanho de amostra, poder e tamanho de efeito, listas de valores em tabelas
- Estatística em `core/power_suites.py`; comando em `commands/power3.py`.
- Conferido com o log do Stata 14: caso `1919` (as diferenças e o motivo de cada uma estão em `compat/known_differences.txt`). Achados e corrigidos no caminho: o `return list` tem ordem própria em cada método, com `r(N_a)` (o tamanho antes de arredondar) e os itens do solucionador, e `r(init)` é a raiz da cauda na direção do efeito em `cmh`, `mcc` e `trend`; `power repeated` divide as variâncias efetivas pelo número de medidas repetidas (`Var_be` = variância x (1 + (M - 1) corr) / M, `Var_w` é a variância das médias) e multiplica a não centralidade por epsilon; `power twoway` com `n()` usa sujeitos inteiros por célula, mostra a variância ponderada com `cellweights()` e lista `Nij` quando as células diferem; `power cmh` toma N como múltiplo da soma dos pesos dos estratos, divide grupos iguais ao meio (fracionário quando N é ímpar), arredonda para cima o grupo experimental com `grratios()`, não mostra a mensagem de iteração com `onesided` e tem a tabela de grupo por estrato no layout do Stata; `power cox` e `power logrank` calculam N a partir dos eventos antes de arredondar e mostram E = teto de N (1 - Pr_w) Pr_E; a largura dos rótulos cresce com o maior rótulo e as tabelas mostram o poder pedido quando o tamanho é calculado
- Ainda não disponíveis nesses métodos: `graph()`, `saving()`, `show`, `st1()`, `n1()` e `n2()`

## Distribuição e documentação

- [x] Licença MIT e avisos de licença dos componentes de terceiros
- [x] Download para macOS (Apple Silicon e Intel) com a Release gerada pelo GitHub Actions a cada tag `v*`
- [x] Toda a documentação publicada no GitHub em inglês, com cópia em português (`*.pt-BR.md`) e o seletor **English | Português** no topo de cada página; comentários do código e mensagens de commit em inglês
- [x] Busca e instalação de atualizações dentro do app (Releases do GitHub, download verificado, verificação diária que pode ser desligada)
- [x] Download para Windows (64 bits, Intel/AMD): instalador por usuário e `.zip` portátil, gerados pelo mesmo workflow
- [x] Versão 1.0: README curto com captura de tela, documentação em `docs/` (primeiros passos, recursos, testes, desenvolvimento), `CONTRIBUTING.md`, modelos de issue (diferença com o Stata, bug, pedido de comando) e de pull request; o workflow de release confere se a tag coincide com a versão

## Melhorias em relação ao Stata (ao longo das fases)

- [x] Explicação dos erros: depois da mensagem e antes do r(#), uma linha diz o que deu errado no contexto do comando e outra sugere como resolver (nomes parecidos, arquivos parecidos na pasta, número de observações…). Em inglês ou português, conforme o idioma da interface; `set hints off` desliga
- [x] Explicações para iniciantes: o que deu errado com os dados à mão (faltantes por variável, variáveis de texto ou numéricas, valores do desfecho, contradições do assert, tamanhos de matrizes), por que a regra existe, o comando corrigido (= no lugar de ==, funções do egen no generate, aspas para texto e rótulos de valor, bysort, nomes sem espaço ou acento...) e um link clicável para a ajuda; novas páginas de ajuda para operadores, funções, valores faltantes e by
- [x] `set numerics precise`: QR no `regress`, `D.`/`S.` em double, passo final de Newton no `ml`/`optimize()`, Richardson nas derivadas numéricas e variáveis novas em double (o padrão `stata` reproduz os arredondamentos do Stata)
- [ ] Opções estatísticas além do Stata: `vce(cr2)` com graus de liberdade de Satterthwaite, *wild cluster bootstrap*, `logit, firth`
- [x] Expoente sobrescrito na notação científica exibida (`1.0000×10¹⁰`, `2.5×10⁻⁵`), mantendo a largura do formato; `set superscript off` (ou Help → Superscript Exponents) volta ao `e+10` do Stata. `string()`, macros, `file write`, exportação e `sprintf()`/`strofreal()` do Mata ficam sempre na forma do Stata
- [x] O Review pode listar os comandos executados pelos do-files (`do`, `run`, `include` e os botões Do e Run do Editor de do-files) sob a linha do do, com os do-files aninhados e os comandos com erro em vermelho, e as seleções do editor identificadas pelo arquivo (Settings > General > Session; desligado por padrão, como no Stata, em que o Review mostra só `do "<arquivo temporário>"`)

## Faltam para a 1.0

Comandos e opções do Stata 14 que o OpenDTA ainda não tem, por área. Todos são necessários para a 1.0; a ordem é a dos blocos de [docs/next-steps.pt-BR.md](docs/next-steps.pt-BR.md), ajustada pelos pedidos abertos com o modelo [Command request](https://github.com/yaskepalhares/opendta/issues/new?template=command-request.md).

- **Manipulação de dados**: `icd9p` com a lista de procedimentos da ICD-9-CM.
- **Modelos**: `ivprobit` e `ivtobit` com `twostep` e com mais de um regressor endógeno; `ivpoisson cfunction` com mais de um regressor endógeno; `ivregress gmm` com `wmatrix(hac kernel opt)`; `biprobit` com `partial` e `constraints()`; `asmprobit` e `asroprobit` com `altwise`, alternativas não balanceadas e outras estruturas de covariância; `slogit` com `dim()` acima de 1; `exlogistic` e `expoisson` com `terms`, `group()`, `binomial()`, `condvars()`, `saving()`, `exposure()` e `offset()`; as formas com programa de `nl`, `nlsur` e `gmm`, e `gmm` com `derivatives()`, `parameters()` e `igmm`; `frontier` com `constraints()`; `teffects` com tratamentos de vários níveis, outros modelos de desfecho e `predict`; `stteffects wra`; `eteffects` com outros modelos de desfecho e de tratamento; `etregress` e `etpoisson` com `cfunc`.
- **Dados em painel**: `xtunitroot` com `llc` e `ips` (precisam dos momentos tabelados dos artigos originais), `trend` no `ht`, `lags()` no `breitung` e as escolhas automáticas de defasagens do `kernel()`; `xtregar` com `lbi`; `xtdata` com `fa`; `xtabond`, `xtdpdsys` e `xtdpd` com `teffects`, `small`, `diffvars()`, `collapse`, `h()` e `vce(cluster)` (algumas podem ser opções do `xtabond2`, de usuário, e não do Stata 14: conferir no manual).
- **Séries temporais**: `vec` com `bconstraints()`, `aconstraints()` e `sindicators()`, e os valores críticos de Osterwald-Lenum para as tendências diferentes de `constant` e `rconstant`; `svar` com `lrconstraints()`; `fcast` com `bs`, `bsp`, `exog` e defasagens não consecutivas; `arfima` com `mpl` e `smemory`; `ucm` com `constraints()`; `sspace` com `covstate()` e `covobserved()` diferentes de `diagonal`, restrições nas variâncias e estados com defasagens acima de um; `dfactor` com mais de um fator, `ma()` e autorregressões nos erros; `mgarch vcc`, `distribution(t)`, `het()`, `constraints()` e `vce(robust)`; `predict` com os estados, os componentes e `smethod()` e `estat` depois de `ucm`, `sspace` e `dfactor`; `predict` e `estat` depois de `mgarch`.
- **Estatística multivariada**: `manovatest` com `test()`; `pca` com `vce(normal)`, `means` e variáveis suplementares; `factor` com `estat factors` e `estat rotatecompare`; `rotate` com `target()`, `matrix()` e os outros critérios (`oblimax`, `bentler`, `cf()`…); `manova` com termos de erro; `canon` com `estat`; `mvtest` com as hipóteses `compound`, `block` e `pairwise` e `normality` com `univariate`, `bivariate`, `by()` e Doornik-Hansen; `discrim` com `ties()`, `priors(matname)`, pesos e `lootable` depois de `logistic`, `estat correlations`, `estat grsummarize` e `estat manova`; `ca` e `mca` com pontos suplementares, `method(indicator)` e `method(johnson)` e os outros subcomandos de `estat`; `mds` com `method(nonmetric)`, `loss()`, `transform()` e `addconstant`; `mdslong` com `s2d()`.
- **Modelos multinível**: `intmethod(laplace)` com mais de um nível; efeitos aleatórios cruzados (`_all:`) nos comandos `me`.
- **Outras suítes**: `power` com `graph()`, `saving()`, `show`, `st1()` e `n1()`/`n2()` em `twoway`, `repeated`, `cmh`, `mcc`, `trend`, `cox`, `exponential` e `logrank`; `power exponential` com os resultados `r(E0_vec)`, `r(Ea_vec)` e `r(Pr_vec)`; `power repeated` com o multiplicador m de Muller e Barton em `covmatrix()` não esférico (o poder difere do Stata na 4ª casa); `bayes` com os comandos de estimação além de `regress`, `logit`, `probit` e `poisson`; `irt` com `group()`; `gsem` com variáveis latentes por observação e por conglomerado no mesmo modelo; `svy:` com `ivprobit`, `ivtobit` e os demais estimadores que ainda não aceitam pesos; `e(b)` de `sem`, `gsem` e `irt` com os parâmetros restritos (e `e(Cns)`); `svy:` guardando os resultados `e()` do estimador.
- **Programação**: o sistema de classes do Stata (as classes do Mata já funcionam).

Ficam de fora da 1.0 de propósito, como limitações conhecidas: `import haver`, que lê os bancos da Haver Analytics só pelo software da própria Haver (Windows e licença paga), e o carregamento de plugins compilados do Stata (`program nome, plugin` e `plugin call`), que são código de máquina escrito para a interface interna do Stata; o mesmo cálculo pode ser escrito num do-file, num ado-file ou em Mata.

## Depois da 1.0

O que vem depois que todos os comandos do Stata 14 estiverem cobertos, em ordem de prioridade. São melhorias em relação ao Stata para inquéritos populacionais de saúde e epidemiologia, mais as versões seguintes do Stata. A maioria não tem equivalente no Stata 14; por isso a referência é a literatura ou o R (`survey`), e não um log do Stata, e cada item diz como será validado. A estatística fica em `core/`, separada da sintaxe (veja [docs/next-steps.pt-BR.md](docs/next-steps.pt-BR.md)).

**P0: amostras complexas e epidemiologia**

1. Calibração de pesos: pós-estratificação, *raking*, calibração linear e logit, aparamento (*trimming*) e diagnóstico dos pesos (tamanho efetivo de amostra de Kish, variabilidade dos pesos).
2. Intervalos de confiança de prevalência para estimativas raras ou extremas (Clopper-Pearson e Wilson com tamanho efetivo de amostra) e um relatório de confiabilidade (CV, graus de liberdade, tamanho efetivo, UPAs por domínio), como aviso na saída do `svy:`.
3. Tabelas de resultados prontas num só comando: Tabela 1 ponderada (`n`, `%`, IC, valor de p, `by()`), regressões lado a lado (como o `esttab`) com razões de prevalência e IC em português ou inglês, exportadas para Word, Excel e LaTeX (reaproveitando o `putexcel`).
4. Razões de prevalência e prevalência ajustada: `binreg` (bloco 10) com `svy:`, e `margins`/`predict` com IC para a prevalência ajustada depois de `poisson`, `glm` e `binreg`.
5. As opções que já estão em [Melhorias em relação ao Stata](#melhorias-em-relação-ao-stata-ao-longo-das-fases): `vce(cr2)` com graus de liberdade de Satterthwaite, *wild cluster bootstrap* e `logit, firth`.
6. Verificações e dicas antes do `svy:`: estratos com uma só UPA, `if` no lugar de `subpop()`, `svyset` ou pesos ausentes (amplia o `core/hints.py`).

**P1: fluxo de trabalho com inquéritos**

7. Importador dos dicionários de inquéritos do IBGE (PNS, POF, PNAD) e do layout do Vigitel: gera o `infile` e os rótulos de valores.
8. Tamanho de amostra para inquéritos: efeito de desenho, não resposta, alocação ótima por estrato, conglomerados (ICC).
9. Índices de desigualdade em saúde (índice de inclinação e índice relativo de desigualdade, índice de concentração) com `svy:`.
10. Não resposta: ponderação pelo inverso da probabilidade de resposta, e `mi` compatível com `svyset`.
11. Escalas e índices com pesos: alfa de Cronbach e ômega, análise fatorial e TRI ponderadas (depois dos blocos 8 e 11), e escores de segurança alimentar e de dieta (EBIA, NOVA, adequação da dieta).
12. Padronização direta e indireta com o desenho amostral, e tendências entre inquéritos (*joinpoint*/APC), com testes de diferença entre domínios e correção para comparações múltiplas.

**P2: relatórios, gráficos e ferramentas**

13. Relatórios reprodutíveis: um do-file que escreve um documento Markdown/HTML/docx com comandos, tabelas e gráficos.
14. Gráficos: exportação e importação de JSON do Plotly (um subconjunto documentado de tipos de gráfico); testes de integridade das exportações SVG, PNG, JPEG e PDF; SVG com texto editável. O OpenDTA mantém o seu formato `.gph` próprio: ler o `.gph` do Stata continua fora de escopo (veja abaixo).
15. Python dentro dos do-files (blocos `python:`, com pandas e scikit-learn).
16. Bootstrap, jackknife e pesos de réplica em paralelo; ambiente e semente registrados em cada resultado.

**P3: mais adiante**

17. Estimação para pequenas áreas (modelos de área de Fay-Herriot).
18. Os comandos e as mudanças de cada versão seguinte do Stata, uma de cada vez, até o Stata 19.
19. Uma versão que roda no navegador (veja Longo prazo).

## Longo prazo

Para registro, a ordem depois da 1.0:

1. Terminar todos os comandos da versão-alvo (os blocos de [docs/next-steps.pt-BR.md](docs/next-steps.pt-BR.md)) e lançar a 1.0. Nada de [Depois da 1.0](#depois-da-10) vem antes.
2. As prioridades de [Depois da 1.0](#depois-da-10) e, em seguida, os comandos e as mudanças de cada versão seguinte do Stata, uma versão de cada vez, até o Stata 19.
3. Uma versão que roda no navegador, sem instalar nada: o núcleo estatístico (Python, NumPy e SciPy) roda em WebAssembly pelo Pyodide. Ela aproveita a separação entre o núcleo estatístico e a camada Stata que os comandos novos já seguem.

Fora do escopo: ler os arquivos `.gph` gravados pelo Stata, cujo formato não tem especificação pública ([política de desenvolvimento](docs/clean-room-policy.pt-BR.md)); o OpenDTA grava gráficos no seu próprio formato `.gph`.

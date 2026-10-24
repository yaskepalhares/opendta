# Próximos passos: fechando o roadmap

[English](next-steps.md) | **Português**

Um plano de trabalho para os comandos que ainda faltam em [Faltam para a 1.0](../ROADMAP.pt-BR.md#faltam-para-a-10), escrito para que cada bloco possa ser feito numa sessão nova sem redescobrir o código. A versão 1.0 sai quando todos os itens de Faltam para a 1.0 estiverem cobertos; depois vêm as melhorias de [Depois da 1.0](../ROADMAP.pt-BR.md#depois-da-10) e as versões mais recentes do Stata.

## Como executar um bloco

Comece uma sessão nova (uma conversa longa custa mais a cada passo) com uma instrução como:

> Implemente o bloco **Programação** conforme docs/next-steps.md. Abra um pull request e faça o merge quando o CI passar.

Todo bloco segue os mesmos passos:

1. Branch a partir do `main` (`block-<nome>`).
2. Para cada comando: a estatística em `src/opendta/core/` (veja a próxima seção), o comando em `src/opendta/commands/` (módulo novo registrado em `commands/__init__.py`), a página de ajuda e os testes.
3. Um caso de comparação `compat/do/NNNN_<bloco>.do` com os comandos e `return list`/`ereturn list`; ele fica sem referência até a rodada no Stata. Rode sob `capture noisily` todo comando cujo comportamento no Stata o manual não documenta (e faça o mesmo em qualquer caso antigo que você alterar), para que uma única rodada no Stata mostre o caso inteiro em vez de parar no primeiro erro. O Stata de referência é a versão 14.0: um caso com `version 14.1` para com erro na rodada; escreva `version 14` ou nenhuma linha `version`.
4. `ROADMAP.md` e `ROADMAP.pt-BR.md`: uma fase nova listando o que foi feito; em *Faltam para a 1.0*, só o que falta. Cada trecho de uma linha de *Faltam para a 1.0* que lista opções deve começar pelo comando a que elas pertencem (`` `cs` com a opção `pool` ``), para o mapa de cobertura não contar opções como comandos.
5. `tools/roadmap_coverage.py`: acrescentar o módulo novo em `MODULE_AREA`.
6. `docs/features.md` e `docs/features.pt-BR.md`: uma linha para o bloco.
7. A rodada no Stata, uma só, para tudo o que o bloco precisa (inclusive casos sem referência): peça ao usuário que rode `compat/gerar_esperados.do`, nunca um bloco improvisado de comandos. Esse arquivo já faz `set dp period` e `linesize 255`, apaga as sobras `odta_*` e roda `mata clear`. Antes de pedir, confira com o git que o branch está no GitHub e diga ao usuário que dê `git pull` nele. Se uma rodada falhar por motivo de ambiente (vírgula decimal, nome de arquivo, sobras de outra rodada), diga isso antes de pedir outra. Antes de pedir, confira também que os dados do caso estão definitivos (os modelos convergem, os testes não são degenerados): mudar os dados depois custa uma segunda rodada. O clone do usuário precisa estar no branch do bloco (uma rodada em outro branch não gera o log); se tiver alterações locais, `git stash -u` antes.
8. `python tools/compare.py` sem diferenças novas; abrir o PR; merge quando os testes passarem no Linux, no macOS e no Windows. A página de cobertura se atualiza sozinha depois do merge.

Notas para sessões que rodam na nuvem:

- O trabalho é feito num clone na nuvem, que alcança o GitHub; o computador do usuário não alcança. Depois de cada push, diga ao usuário que dê `git pull` no branch. Use `git push -u origin <branch>` no primeiro push, para o git conhecer o branch remoto.
- O `pytest` pode não instalar lá (sem acesso ao PyPI): os testes que dependem de pacotes opcionais (`pyreadstat`, `xlrd`, `xlwt`) só rodam no CI. Espere o CI e leia o resultado antes de pedir a rodada no Stata.
- Os logs do CI não baixam da nuvem, mas as anotações de cada job mostram as diferenças do `tools/compare.py` (`gh api repos/<dono>/<repo>/check-runs/<id>/annotations`). O `gh pr` falha porque o GraphQL está bloqueado: abra e mescle o pull request com `gh api` (REST), com merge commit.
- O log de referência da rodada no Stata nunca é editado à mão: quando um caso muda, a rodada no Stata é repetida.
- A anotação do CI de um caso com diferença conhecida mostra só os trechos que os limites `ignore()` e `tolerance()` não aceitam. Ajustes na fronteira (por exemplo, uma variância de fragilidade em 0) podem diferir entre Linux, macOS e Windows: espere os três jobs do CI antes de fechar os limites de uma diferença conhecida e mantenha-os tão estreitos quanto as linhas que realmente diferem.
- Depois da rodada no Stata, leia o log antes de comparar: copie-o do clone do usuário (`compat/stata/NNNN_<bloco>.log`, que o git ignora) para `compat/expected/`, sem alterar. Não edite o `.do` depois: ele deve continuar idêntico ao que o Stata executou. Compare em duas passagens: primeiro o layout (troque todo número por `#`, e aparecem títulos, rótulos, listas `e()` e a ordem delas); depois os números, ignorando as linhas `Iteration` e `Grid node`; essa segunda passagem acha diferenças reais, como uma constante na verossimilhança ou um parâmetro multiplicado por 2, que o ruído das iterações esconde.
- Não construa os dados do caso sobre a ordem de observações empatadas (`bysort g: gen ... if _n == 1` seguido de `runiform()`): a ordem dos empates no Stata é arbitrária, então os sorteios seguintes vão para outras observações e todas as estimativas diferem. Use `sort g, stable` ou gere por `_n` antes de ordenar.
- No shell da nuvem, não use `pkill -f` com um padrão que também casa com a linha de comando do próprio shell (mata a sessão); e não deixe um arquivo com nome de módulo padrão (`struct.py`) na pasta de um script: rode os scripts com `python -I`.
- Os conjuntos de dados de exemplo do Stata (`auto`, `nlsw88`...) não existem: um caso cria os próprios dados com `set seed`, `runiform()` e `invnormal()`, arredondados para inteiros ou centavos, de modo que o último dígito dos sorteios não importe.
- Um comando que roda outro sob `capture noisily` cancela um `quietly` interno: o código que roda um comando em silêncio (os ajustes dentro de `fp`, `nestreg`, `jackknife`) usa `Output.quietly()`.
- Sem pytest, rode as funções de teste com algumas linhas próprias (importe o módulo e chame cada `test_*` com o `Runner` de `tests/conftest.py`). Falhas que também ocorrem em `origin/main` vêm do ambiente (falta de `xlwt`, fixtures com `tmp_path`), não do bloco.
- Um caso com centenas de comandos deixa o `tools/compare.py` completo lento: ao corrigir, rode um trecho do caso e compare a saída de cada comando com o mesmo comando do log de referência; rode o `compare.py` uma vez no final.
- Números mostrados com `%18.0g` têm 16 dígitos significativos (`.81` aparece como `.8100000000000001`): use `format_value(x, "%18.0g")`, não `repr`. A contagem de iterações (`e(ic)`) de uma busca pode variar por plataforma: ignore a linha com uma faixa (`(3|1[0-9])`).
- Permutações e simulações: depois de `permute` ou `simulate` o gerador fica em outro estado que o do Stata, então tudo o que sorteia depois difere; liste essas linhas na diferença conhecida.
- `rnormal()` não dá os sorteios do Stata: crie dados normais com `round(invnormal(runiform()), .01)`.
- Uma primeira rodada no Stata mostra sobretudo layouts e numéricos, e um comando que falha com dados degenerados desperdiça a sua seção: dê a cada modelo dados que ele consiga estimar (caudas pesadas para `arch` com erros `t`, uma série com o que um modelo de componentes estima: ciclo, tendência, sazonalidade, correlação móvel) e mantenha as seções do caso independentes, cada uma sob `capture noisily`, para que uma rodada resolva o bloco.
- Para descobrir o que é uma diferença numérica, salve os dados do caso em csv e teste hipóteses fora do Stata: avalie a verossimilhança do OpenDTA nas estimativas de referência. Uma diferença constante (por exemplo `ln 4`) é uma constante da verossimilhança; uma diferença que some quando um valor pré-amostral vira zero ou a média da amostra é uma regra de pré-amostra; uma verossimilhança maior no OpenDTA nas mesmas estimativas, com estimativas diferentes, é o otimizador ou outro modelo, e vai para a diferença conhecida.
- Numa diferença conhecida, `ignore(/REGEX/)` tira linhas dos dois logs, então linhas que existem em só um deles (logs de iteração do Stata) não atrapalham; uma tabela que aparece várias vezes não se delimita por posição, então case os valores que diferem e deixe o `tolerance()` para os últimos dígitos.
- Faça o push sem `-q` e confirme com `git ls-remote --heads origin <branch>`. O aviso "commits não enviados" do hook de parada pode ser uma referência de rastreamento velha: `git fetch origin +refs/heads/<branch>:refs/remotes/origin/<branch>` resolve.
- Para descobrir por que o CI falha quando os logs não baixam, leia as anotações job a job, rode os arquivos de teste um por vez e defina `faulthandler_timeout`. O CI de um PR roda no merge com o `main`; um travamento em teste de GUI costuma ser um `QMessageBox` modal.
- Linhas que dependem da plataforma (uma verossimilhança plana numa fronteira, uma Hessiana numérica) diferem entre Linux, macOS e Windows: ignore-as com `ignore(/REGEX/)`. Quando a saída inteira de um comando é instável ou o comando ainda não existe, use `skip(/REGEX/)` na linha `. comando` dele, que tira a saída desse comando dos dois logs.
- Depois de uma edição por script que corta um arquivo, confira a lista de funções e o `git diff --stat`: uma âncora que casa com uma função anterior apaga código sem aviso.
- Achados de painel do caso 1905: os comandos do tipo `xtabond` estimam a constante por uma equação em nível separada (instrumento só `_cons`); o sigma2 de um passo é a soma dos quadrados dos resíduos diferenciados sobre `2(n_d - k)`, com `k` contando a constante; o `xtivreg, fd` tem constante, e as estatísticas dele vêm dos níveis das linhas que têm defasagem.
- Achados de multivariada do caso 1906: uma matriz que o comando guarda depois de um produto BLAS (correlação ou dissimilaridade) deve ficar exatamente simétrica (`np.tril(D) + np.tril(D, -1).T`): no macOS o ruído de arredondamento transformava a listagem `symmetric D[5,5]` do Stata numa listagem completa. `factor, ml` para depois de poucos passos de Newton, cerca de 1e-7 do máximo, e as cargas diferem na 4ª casa decimal (diferença conhecida com `ignore(/^Iteration \d+: /)` e `tolerance(2e-4)`). `procrustes`: df_m = ny*nx - nx(nx+1)/2 + ny + 1 (menos 1 com `norho`) e o RMSE de cada variável alvo usa df_r/ny. A ordem do `ereturn list` é a ordem das listas dadas a `post_e` (`commands/mvcommon.py`); o sinal dos autovetores precisa de uma regra explícita por comando (primeira coordenada não nula positiva no `mca`, maior coeficiente absoluto no `canon`).
- Logo depois de um push, `gh api repos/<dono>/<repo>/commits/<sha>/check-runs` vem vazio até os jobs serem criados: um laço que espera `status != completed` termina na hora. Espere até existirem os três jobs.
- Achados dos modelos (caso 1908): monte dados arredondados com coeficientes em décimos pares e evite comparações com zero exato (`round(soma, .01)` cria metades exatas que o Stata e o OpenDTA podem arredondar de modo diferente); evite `rpoisson()`, `rbinomial()` e `rnormal()` e parâmetros na fronteira (uma verossimilhança plana difere por plataforma). Confira que os dados do caso são idênticos nos dois programas antes de gastar a rodada do Stata, e veja o que o Stata rejeita (`nneighbors()` não é opção, `stteffects ipw` não aceita opção de modelo de sobrevivência e exige que toda observação tenha a falha). Primeiro o layout, depois os números: uma hipótese para uma fórmula (variância do `eivreg`, escala do `rreg`) se testa fora do programa num csv dos dados do caso contra o log de referência, e um algoritmo não documentado que não fecha em poucas tentativas vai para a diferença conhecida com o motivo. Uma falha sem anotação num job do CI pode ser uma exceção dentro de um caso: o `tools/compare.py` agora a mostra como `XX` com anotação.

Próximos números livres: casos de comparação a partir de **1910**, fases do roadmap a partir de **24**.

## Regra: a estatística no núcleo, o Stata no comando

O cálculo vai em funções simples em `src/opendta/core/` (ou em `src/opendta/estimation/`, para modelos): recebem arrays e números e devolvem números, arrays ou pequenos objetos de resultado, sem `Session`, sem ler variáveis e sem imprimir. O módulo do comando interpreta a sintaxe, lê os dados, chama essas funções, escreve a saída do Stata e guarda `r()`/`e()`. Exemplo: `core/epi.py` (fórmulas) e `commands/epitab.py` (tabelas). Os testes conferem as funções do núcleo contra fórmulas fechadas ou o SciPy, e o comando contra a saída.

## O que reaproveitar

| Necessidade | Onde |
|---|---|
| sintaxe: varlist, if/in, pesos, opções | `lang/syntax.py` (`parse_standard`, `match_options`, `find_top`), `commands/_util.py` (`touse`, `eval_vector`, `plural`), `commands/summarize.py` (`weights`, `_level`, `g9`, `name12`) |
| grupos de `by()` e estratos | `commands/epitab.py` (`_strata`), `core/grouping.py` (`group_ids`) |
| `r()` e matrizes | `s.r = {...}` (a ordem de inserção é a do `return list`, invertida onde indicado), `commands/matrix.py` (`Matrix`, `evaluate_matrix`) |
| resultados e tabela de estimação | `estimation/results.py` (`Estimates`, `post`, `coef_table`, `CoefRow`) |
| máxima verossimilhança | `estimation/ml.py` (`maximize`, `print_log`): os passos de Newton e o log de iterações do Stata |
| previsões depois de um modelo | `estimation/postest.py` (`xb_all`, `esample`), `estimation/estat.py` (`_binary`) |
| prefixos (como o bootstrap) | `estimation/bootstrap.py` |
| gráficos | `graphics/model.py` (`Graph`, `Layer`, `Panel`), `commands/graphcmd.py` (`finish`), `graphics/opts.py` (`graph_options`) |
| páginas de ajuda | `src/opendta/help/*.sthlp`, `ALIASES` em `commands/helpcmd.py`; `tests/test_i18n.py` confere que todo comando tem uma |
| erros | `core/errors.py` (`StataError(rc, mensagem)`); uma opção ainda não implementada dá `StataError(198, "option X is not available yet in OpenDTA")` em vez de ser ignorada |

Os testes usam a fixture `run` de `tests/conftest.py` (`run("comandos")` devolve a saída; `run.session.r`, `run.session.e`). O ambiente na nuvem pode não instalar o pytest: nesse caso, o CI roda os testes.

## Blocos

Estimativas a partir do trabalho feito até aqui (±50%): tempo de sessão e tokens gerados; uma sessão relê o contexto a cada passo, então o consumo total fica entre 100 e 150 vezes os tokens gerados. A ordem é uma sugestão: primeiro os blocos baratos que fecham uma área, depois os mais usados em pesquisa em saúde.

### 1. Epidemiologia: completar ✅

Feito na fase 12 (caso de comparação 1601), menos o `rocreg` com `ml` e `cluster()`, que ficam em *Faltam para a 1.0*. O plano original:

`rocreg` (regressão ROC não paramétrica e paramétrica; erros padrão por bootstrap); `cs`, `cc` e `ir` com `estandard`, `istandard`, `standard()`, `binomial` e `pool`; `mhodds` com `by()`; `roctab` com `binomial`, `bamber` e `lorenz`; `roccomp` com `binormal`. Reaproveitar `core/epi.py`, `commands/epitab.py`, `commands/roc.py`; tirar os erros "not available yet" correspondentes.

### 2. Manipulação de dados: opções ✅

Feito na fase 13 (caso de comparação 1701), exceto `import haver` (formato proprietário), `icd9p` (não foi possível obter uma cópia da lista de procedimentos da ICD-9-CM do CMS), `import sasxport, describe` e `member()`, os rótulos de valores do `formats.xpf`, `export sasxport` com `rename` e `vallabfile()`, e `list` com `string()` e `display`, que ficam em *Faltam para a 1.0*. O plano original:

`import sas`, `import spss` e `import sasxport` (pelo `pyreadstat`, opcional como para os `.dta` antigos; `io/`), `file` com leitura e gravação binária (`commands/filecmd.py`), `infile` com `byvariable()` (`commands/infile.py`), `putexcel` para `.xls`, as opções que faltam de `list` e `describe`. `import haver` precisa de um formato proprietário: deixar de fora e dizer isso. `icd9p` precisa da lista de procedimentos da ICD-9-CM, de domínio público (CMS); acrescentá-la em `src/opendta/data/icd/` com `tools/make_icd_tables.py` se uma cópia puder ser obtida.

### 3. Programação ✅

Feito na fase 14 (caso de comparação 1801), exceto `makecns` com lista de restrições e `displaycns` (dependem de `constraint`, do bloco 10), `matlist` com as opções de equações, as caixas de diálogo e os menus de `window` na interface gráfica e o carregamento de plugins, que ficam em *Faltam para a 1.0*. O plano original:

`unab`, `tsunab`, `fvunab`, `unabcmd`, `fvexpand`, `fvrevar`, `tsrevar`, `mark`, `_rmcoll`, `matlist`, `tabdisp`, `sortpreserve` (opção de program), `signestimationsample`, `makecns`, `postfile`/`post`/`postclose`, `pause`, `sleep`, `window` (comandos de diálogo, até onde a interface permitir), `plugin` (não carrega plugins do Stata: explicar), `viewsource` (mostrar o código do próprio OpenDTA, não o do Stata). Quase tudo em `lang/` e `commands/program.py`, `commands/progcmd.py`, `commands/syntaxcmd.py`; variáveis de fator em `estimation/fvars.py`.

### 4. Análise de sobrevivência ✅

Feito na fase 15 (caso de comparação 1901), exceto as opções que param com o erro 198 (`stcurve` com `unconditional`, `outfile()` e `kernel()`; `stcrreg` com `tvc()`; vce(robust) e pesos com `ggamma` ou `frailty()`; as opções gráficas de `ltable`, `strate` e `stci`). O plano original: `stsplit`, `stjoin`, `stfill`, `stgen`, `stbase`, `stvary`, `stci`, `strate`, `stptime`, `stir`, `stmh`, `stmc`, `ltable`, `sttocc`, `sttoct`, `snapspan`, `ct`, `stcurve` (gráfico depois de `stcox`/`streg`), `stcrreg` (riscos competitivos de Fine-Gray, ML); `streg, distribution(ggamma)` e `frailty()`; os empates `exactm` e `exactp` do `stcox`. Reaproveitar `estimation/st.py`, `km.py`, `stcox.py`, `streg.py`, `stfit.py`.

### 5. Modelos multinível (1–1,5 h, ~55 mil tokens)

Feito nas fases 16 e 30 (casos 1902 e 1916), exceto as peças listadas em Modelos multinível em [Faltam para a 1.0](../ROADMAP.pt-BR.md#faltam-para-a-10) (`meqrlogit` e `meqrpoisson` com mais de um nível, `intmethod(laplace)` com mais de um nível, efeitos cruzados nos comandos `me`). O plano original: `meglm` (a forma geral), `mecloglog`, `menbreg`, `meologit`, `meoprobit`, `mestreg`, `meqrlogit` e `meqrpoisson`; efeitos aleatórios cruzados no `mixed`. Reaproveitar `estimation/meglm.py` (quadratura adaptativa) e `mixed.py`.

### 6. Dados em painel ✅

Feito na fase 19 (caso 1905), exceto as peças listadas em Dados em painel em [Faltam para a 1.0](../ROADMAP.pt-BR.md#faltam-para-a-10) (`xtunitroot llc` e `ips`, que precisam dos momentos tabelados dos artigos originais, e algumas opções dos comandos de painel dinâmico, `xtregar`, `xtrc`, `xtivreg` e `xtdata`). O plano original: modelos de efeitos aleatórios por quadratura, `quadchk`, painéis dinâmicos por GMM, os estimadores lineares de painel, `xtfrontier`, `xtunitroot` e os comandos de descrição.

### 7. Estatística básica (3–4 h, ~135 mil tokens)

Feito na fase 17 (caso 1903), exceto as peças listadas em Estatística básica em [Faltam para a 1.0](../ROADMAP.pt-BR.md#faltam-para-a-10) (`ml` com os métodos `lf0`, `lf1`, `lf2` e `gf`, `boxcox` com outros modelos, estatísticas de teste do `cusum`…). O plano original: prefixos e reamostragem: `jackknife`, `permute`, `simulate` (reaproveitar `estimation/bootstrap.py`); `ml` (o comando em torno de `estimation/ml.py`), `mlexp`; construção de modelos: `nestreg`, `stepwise`, `xi`, `mkspline`, `fp`, `mfp`; depois da estimação: `contrast`, `pwcompare`, `pwmean`, `margins` depois de mais comandos (`estimation/margins.py`), `suest`, `predictnl`, `testnl`, `linktest`; transformações: `boxcox`, `ladder`, `gladder`, `lnskew0`; testes: `bitest`, `ksmirnov`, `median`, `nptrend`, `runtest`, `tetrachoric`, `ktau`, `pcorr`; concordância e tamanhos de efeito: `icc`, `loneway`, `esize`, `ameans`; exibições: `cumul`, `stem`, `dotplot`, `lv`, `lpoly`, `smooth`, `cusum`, `serrbar`, `sunflower`, `symplot`, `quantile`.

### 8. Estatística multivariada ✅

Feito na fase 20 (caso 1906), exceto as peças listadas em Estatística multivariada em [Faltam para a 1.0](../ROADMAP.pt-BR.md#faltam-para-a-10) (as medidas binárias de similaridade do `cluster`, pontos suplementares em `ca` e `mca`, `mds` com `method(nonmetric)`, `procrustes` com transformações oblíquas e irrestritas, e alguns subcomandos de `estat`). O plano original: `pca`, `factor`, `rotate`, `alpha`, `cluster` (ligações, dendrograma, `cluster generate`), `discrim` (lda, qda, knn, logistic), `manova`, `mvreg`, `mvtest`, `canon`, `ca`, `mca`, `mds`, `biplot`, `procrustes`, `hotelling`, `screeplot`, `scoreplot`, cada um com seus subcomandos de `estat`. Módulo novo `core/multivariate.py`.

### 9. Séries temporais (3–4 h, ~150 mil tokens)

Feito na fase 18 (caso 1904), exceto as peças listadas em Séries temporais em [Faltam para a 1.0](../ROADMAP.pt-BR.md#faltam-para-a-10) (restrições em `vec`, `svar`, `arch`, `ucm` e `sspace`, `mgarch vcc`, `predict` depois dos comandos de espaço de estados…). O plano original: `arch`, `arfima`, `mgarch`, `dfactor`, `sspace`, `ucm` (espaço de estados: compartilhar um único filtro de Kalman com o `arima` em `estimation/tsmodels.py`), `vec`, `vecrank`, `svar`, `fcast`, `varbasic`, `varlmar`, `varnorm`, `varstable`, `varwle`, `dfgls`, `pperron`, `prais`, `rolling`, `tsappend`, `tsfill`, `tsreport`, `tsfilter`, `tssmooth`, `pergram`, `cumsp`, `xcorr`, `wntestb`, `wntestq`, `psdensity`; `arima` com `sarima()`, `condition` e `diffuse`; `irf` depois de `svar` e `vec`, com erros padrão por bootstrap, e `irf ctable`, `cgraph`, `ograph`, `add` (`estimation/irf.py`).

### 10. Modelos (pronto)

Feito nas fases 21 e 22 (casos de compatibilidade 1907 e 1908), exceto o que consta em Modelos em [Faltam para a 1.0](../ROADMAP.pt-BR.md). A parte 10a (binários e contagem, censura e seleção, variáveis instrumentais) está em `estimation/binmodels.py`, `censmodels.py` e `ivmodels.py`, sobre o motor `estimation/mleq.py` e as verossimilhanças de `estimation/likelihoods.py`. A parte 10b: sistemas em `estimation/sysmodels.py` (álgebra em `core/syseq.py`), modelos de escolha em `estimation/choice.py` e `choiceprobit.py` (`core/choice.py`), métodos exatos em `estimation/exact.py` (`core/exact.py`), modelos não lineares em `estimation/nonlin.py` e `frontier.py`, e efeitos de tratamento em `estimation/teffects.py` e `treatment.py` (`core/treatfx.py`).

### 11. Outras suítes (pronto)

Feito na fase 23 (caso de compatibilidade 1909), exceto o que consta em Outras suítes em [Faltam para a 1.0](../ROADMAP.pt-BR.md). O ajuste do `sem` está em `estimation/semfit.py` (métodos `ml`, `mlmv`, `adf`; dados normais com valores ausentes em `core/mvnmiss.py`), com os grupos em `estimation/semgroup.py` e os subcomandos de `estat` em `estimation/semestat.py`; o `gsem` em `estimation/gsem.py` (quadratura adaptativa, respostas "custom"), com o `irt` sobre ele em `estimation/irt.py`; o `bayesmh` e os comandos `bayes` em `estimation/bayes.py`, sobre o amostrador de `core/mcmc.py`; os novos métodos do `power` e o gráfico em `commands/power2.py`.

## Antes da 1.0

Os blocos 1 a 11 saíram na versão 0.5.0. A versão 1.0 espera o que ainda consta em [Faltam para a 1.0](../ROADMAP.pt-BR.md#faltam-para-a-10): planeje como novos blocos (a partir do bloco 12), uma área por vez, executados como os blocos acima.

Quando essa lista estiver vazia: a rodada no Stata de todos os casos de comparação sem referência, as diferenças corrigidas, a versão 1.0.0 em `src/opendta/__init__.py` e `pyproject.toml`, as notas de lançamento escritas e a tag `v1.0.0` enviada (veja [development.pt-BR.md](development.pt-BR.md#publicar-uma-versão)).

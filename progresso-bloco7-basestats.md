# Progresso: bloco "Base statistics" (docs/next-steps.md, seção 7 no arquivo)

Arquivo local (não vai para o GitHub). Se você mandar "prossiga", o Claude lê este arquivo,
confere com o clone e continua da primeira etapa não concluída.
(Blocos anteriores arquivados: progresso-bloco1-compat.md ... progresso-bloco5-multilevel.md.)

Nota: você escreveu "bloco 2", mas "Base statistics" é o bloco 7 em docs/next-steps.md (o 2 é
"Data management", já feito). Segui o NOME do bloco. Corrija se quis outro.

Trabalho feito num clone na nuvem (o shell do seu Mac não alcança o GitHub); branch `block-base-stats`.
Seu clone local está em `block-multilevel` (já mesclado): depois do primeiro push,
`git fetch && git checkout block-base-stats` (se pedir, `git stash -u` antes).
Stata de referência: 14.0 (nunca `version 14.1` nos casos). Próximos números livres: caso 1903, fase 17.

## Etapas

- [x] 1. Preparar: ler next-steps/CLAUDE.md, clone na nuvem (main 4461371); levantamento: nenhum comando do bloco existe ainda
- [x] 2. Branch block-base-stats + explorar código a reutilizar (bootstrap.py, ml.py, margins.py, nptests.py, graphics)
- [ ] 3. Implementar, em subetapas, cada uma enviada ao GitHub:
  - [x] 3a. prefixos: jackknife, permute, simulate (estimation/resample.py; enviado)
  - [x] 3b. ml (lf, d0, d1, d2 + mleval/mlsum/mlvecsum/mlmatsum; init, search, maximize, display, score, clear) e mlexp (estimation/mlcmd.py; enviado). Não feitos (erro 198): lf0/lf1/lf2/gf*, ml check/query/report/graph/plot/trace/footnote/count, mlmatbysum, mlopts, constraints
  - [x] 3c. construção de modelos: nestreg, stepwise, xi, mkspline, fp, mfp (commands/modelbuild.py, estimation/fpoly.py; enviado)
  - [x] 3d. pós-estimação: contrast, pwcompare, pwmean, suest, predictnl, testnl, linktest; margins após tobit/ologit/oprobit (estimation/postest2.py; enviado). Limites: suest só após regress/logit/probit/poisson; contrast sem p./@ e só até 2 fatores; pwcompare sem groups/sort; pwmean sem unequal
  - [ ] 3e. transformações: boxcox, ladder, gladder, lnskew0
  - [ ] 3f. testes: bitest, ksmirnov, median, nptrend, runtest, tetrachoric, ktau, pcorr
  - [ ] 3g. concordância/efeito: icc, loneway, esize, ameans
  - [ ] 3h. gráficos/displays: cumul, stem, dotplot, lv, lpoly, smooth, cusum, serrbar, sunflower, symplot, quantile
  - [ ] 3i. páginas de ajuda
- [ ] 4. Testes pytest, caso compat 1903, ROADMAP fase 17 (en + pt-BR), MODULE_AREA, features (en + pt-BR), next-steps
- [ ] 5. pytest + compare.py, push e PR, CI verde
- [ ] 6. Pedir UMA rodada no Stata (git pull no branch + compat/gerar_esperados.do)
- [ ] 7. Comparar o log do 1903, corrigir, gravar em compat/expected/, known_differences (só algoritmo não documentado)
- [ ] 8. Merge do PR (gh api, merge commit)
- [ ] 9. PR pequeno em "How to run a block" (en + pt-BR), se houver aprendizado novo
- [ ] 10. Instrução do bloco 8 (Multivariate statistics)

## Feito

- 01h20: tarefa recebida. Progresso do bloco 5 arquivado em progresso-bloco5-multilevel.md. Clone na nuvem em /home/claude/opendta (main 4461371). Li "How to run a block" e a seção do bloco. Levantamento: nenhum dos ~45 comandos existe (só margins/jackknife em svy). Seu clone local está em block-multilevel, defasado do main.

- 01h35: 3a pronta e enviada (branch block-base-stats). Teste local sem pytest: executor próprio (do-files). Layouts de jackknife (df=N-1, F), permute (notas) e simulate com VERIFY para a rodada única. Próximo: 3b ml/mlexp.

- 02h10: 3b pronta e enviada. Conferido: ml lf/d0/d1 e mlexp coincidem com logit/MLE exato. Linhas initial/rescale/Wald/cabeçalho com VERIFY (um único log do Stata decidirá). Próximo: 3c nestreg, stepwise, xi, mkspline, fp, mfp.

- 03h00: 3c e 3d prontas e enviadas (testadas por do-file local; layouts com VERIFY para o log do Stata). Próximo: 3e boxcox/ladder/gladder/lnskew0, 3f testes, 3g icc/loneway/esize/ameans, 3h displays.

- 3e/3f/3g prontas e enviadas (commit c9b9fbd): bitest(i), ksmirnov, median, nptrend, runtest, tetrachoric, ktau, pcorr, ameans, icc, loneway, esize(i), lnskew0, ladder, boxcox (só model(lhsonly)). Layouts em VERIFY.
- Próximo: 3h (cumul, stem, lv, smooth, gladder, dotplot, lpoly, serrbar, sunflower, symplot, quantile, cusum), depois 3i (help), passo 4 (testes, caso 1903, docs).

- 3h, 3i, passo 4 prontos e enviados: gráficos/displays, páginas de ajuda, ROADMAP fase 17 (en+pt), next-steps, features, mapa de cobertura, caso compat/do/1903_base_stats.do (dados próprios; auto.dta não existe no OpenDTA), tests/test_basestats.py (~40 testes). Também corrigi: Output.quietly() (capture noisily anulava quietly interno), contrast agora guarda r().
- Em andamento: passo 5 (compare.py completo; depois PR e CI verde nos 3 sistemas).
- Depois: passo 6 (pedir UMA rodada do gerar_esperados.do no Stata 14.0), 7 (comparar o log e corrigir), 8 (merge), 9 (PR de aprendizados do processo), 10 (instrução do próximo bloco).

- Passo 5 concluído: PR #55 aberto, CI verde (Linux, macOS, Windows) no último commit.
- Passo 6 em andamento: pedir UMA rodada do compat/gerar_esperados.do no Stata 14.0 (clone na branch block-base-stats). Quando o usuário escrever "prossiga" com o log em compat/stata/1903_base_stats.log: passo 7 (copiar o log sem editar para compat/expected/, comparar em duas passadas, corrigir o código, known_differences só para algoritmo não documentado), depois 8 (merge), 9, 10.

- Passo 7 em andamento (log do Stata recebido e copiado para compat/expected/1903_base_stats.log). Corrigidos e iguais ao Stata: ameans, bitest/bitesti, ksmirnov, median, nptrend, runtest, tetrachoric, ktau, pcorr, icc, loneway, esize/esizei, ladder, reshape long (abrevia r1 r2 ... r4), stem, cumul, lv, smooth, lpoly, cusum. lnskew0 e boxcox: layout igual, números diferem por tolerância (known_differences).
- Falta no passo 7: sunflower (tabela de bins), bloco ml/mlexp, jackknife, permute, simulate, mkspline, xi, nestreg, fp/mfp, testnl, linktest, contrast, pwcompare, pwmean, margins, ologit/oprobit; depois known_differences.txt, compare.py completo, pytest, push, CI, merge (passo 8).
- Próximo bloco decidido pelo usuário: séries temporais (bloco 9 do next-steps), não Estatísticas multivariadas.
- Passo 7 quase concluído (nuvem, ainda sem push): ml/mlexp (initial/alternative/rescale, Log pseudolikelihood, e() na ordem do Stata, mlexp convergindo), jackknife, permute, mkspline, xi, nestreg, stepwise, fp/mfp, testnl, linktest, contrast, pwcompare, pwmean, margins após tobit, runtest iguais ao Stata. Restam só diferenças de algoritmo não documentado, registradas em known_differences.txt (1903): passo inicial do ml/ologit, busca do lnskew0, pseudosigma do lv (profundidade 1), sorteios de permute/simulate, tabela de caixas do sunflower, e(ic) do boxcox. pytest (mini-executor) sem regressões; falta compare.py completo, commit/push, CI e merge (passo 8).
- Passos 7, 8 e 9 concluídos: PR #55 (bloco Estatística básica) mergeado com CI verde nos 3 sistemas; PR #56 (aprendizados em "How to run a block", en + pt-BR) mergeado. Falta só o passo 10 (instrução do bloco de séries temporais). Faça `git checkout main && git pull` no seu clone.

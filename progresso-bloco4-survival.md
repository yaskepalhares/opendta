# Progresso: bloco 4 "Survival analysis" (docs/next-steps.md)

Arquivo local (não vai para o GitHub). Se você mandar "prossiga", o Claude lê este arquivo,
confere com o clone e continua da primeira etapa não concluída.
(Blocos anteriores: progresso-bloco1-compat.md, progresso-bloco2-datamgmt.md, progresso-bloco3-programming.md.)

Trabalho feito num clone na nuvem (o shell do seu Mac não alcança o GitHub); branch `block-survival`.
Seu clone local fica atrás: depois de cada etapa enviada, `git pull` no branch.
Stata de referência: 14.0 (nunca `version 14.1` nos casos).

## Etapas

- [x] 1. Preparar: clone na nuvem a partir do main (9922499), branch block-survival (criado, ainda sem commit)
- [x] 2. Explorar o código e planejar (manuais públicos do Stata 14 lidos para todos os comandos do bloco)
- [x] 3. Implementar os comandos (core + commands + páginas de ajuda), em subetapas, cada uma enviada ao GitHub:
  - [x] 3a. dados st: stsplit, stjoin, stfill, stgen, stvary, stbase, snapspan (commit enviado ao branch block-survival; ainda sem ajuda/testes)
  - [x] 3b. taxas: strate, stptime, stci, stir, stmh, stmc, ltable
  - [x] 3c. count-time e casos-controle: ct/ctset, cttost, sttoct, sttocc
  - [x] 3d. stcurve e stcrreg (Fine-Gray)
  - [x] 3e. streg: distribution(ggamma), frailty()/shared(); stcox: exactm, exactp
  - [x] 3f. páginas de ajuda de tudo isso
- [x] 4. Testes pytest, caso compat 1901, ROADMAP fase 15 (en + pt-BR), MODULE_AREA, features (en + pt-BR), next-steps
- [x] 5. pytest + compare.py, push e PR, CI verde
- [ ] 6. Pedir UMA rodada no Stata (git pull no branch + compat/gerar_esperados.do)
- [ ] 7. Comparar o log do 1901, corrigir, gravar em compat/expected/, known_differences (só algoritmo não documentado)
- [ ] 8. Merge do PR
- [ ] 9. PR pequeno em "How to run a block" (en + pt-BR), se houver aprendizado novo
- [ ] 10. Instrução do bloco 5 (Multilevel models)

## Feito

- 16h57: tarefa recebida. Progresso do bloco 3 arquivado em progresso-bloco3-programming.md. Clone na nuvem criado a partir do main (9922499, já com os PRs #46 e #47), branch block-survival. Seu clone local está no branch block-programming (já mesclado): antes de seguir com este bloco, depois do primeiro push, `git fetch` + `git checkout block-survival`.

- 17h10: subetapa 3a pronta e enviada (branch block-survival, primeiro push funcionou). Código em src/opendta/commands/stdata.py, core/survival.py e commands/stcmds.py (registro). Faltam ajuda, testes pytest e o caso compat. Próximo: 3b (strate, stptime, stci, stir, stmh, stmc, ltable).

## Decisões / achados do plano

- Numeração: o próximo caso compat é 1901 (fase 15 do ROADMAP): a linha "Next free numbers" do next-steps estava desatualizada (dizia 1801/14); corrijo no PR.
- Módulos novos: core/survival.py (contas puras), commands/stdata.py (stsplit/stjoin/stfill/stgen/stvary/stbase/snapspan), commands/strates.py (strate/stptime/stci/stir/stmh/stmc/ltable), commands/ctdata.py (ct/ctset/cttost/sttoct/sttocc), estimation/stcrreg.py (stcrreg + stcurve); streg.py (ggamma, frailty) e stcox.py (exactm/exactp) ganham as opções. Registro em commands/stcmds.py; MODULE_AREA ganha os módulos.
- stir reaproveita as tabelas de `ir` (commands/epitab.py); stmh/stmc reaproveitam mh_irr de core/epi.py.
- Os manuais públicos não dão o layout exato de stmh/stmc, stfill, stgen (várias funções), stbase sem at(), cttost, stjoin: implemento a melhor leitura, marco VERIFY e o caso 1901 leva sondas com `capture noisily`; a rodada no Stata decide.
- Sem pytest na nuvem se PyPI estiver bloqueado: runner próprio e o CI faz o pytest completo (como no bloco 3).
- Observação: o commit e56ca22 do bloco 3 incluiu odta_r2.dta e odta_results.dta na raiz do repositório (sobras de teste); não mexo nisso neste bloco, só aviso.
- 17h19: subetapa 3b pronta e enviada (commands/strates.py: strate, stptime, stci, stir, stmh, stmc, ltable; commit d07f8ef). Layouts não documentados marcados VERIFY. Próximo: 3c (commands/ctdata.py: ct/ctset, cttost, sttoct, sttocc).
- 17h22: subetapa 3c pronta e enviada (commands/ctdata.py: ctset/ct, cttost, sttoct, sttocc; commit d6381f8). Também corrigi o relatório do stset com fweight ("physical observations remaining, equal to / weighted observations, representing", conforme manual do cttost) em estimation/st.py. Próximo: 3d (estimation/stcrreg.py: stcrreg e stcurve).
- 17h28: subetapa 3d pronta e enviada (estimation/stcrreg.py Fine-Gray; estimation/stcurve.py; commit ecf464f). Variância robusta conferida contra bootstrap (0,084 vs 0,080). stcurve, hazard após stcox e tvc()/texp()/offset() do stcrreg ainda dão erro 198 "not available yet". Próximo: 3e (streg ggamma e frailty; stcox exactm/exactp).
- 17h49: subetapa 3e pronta e enviada (estimation/streg_ext.py ggamma, frailty(gamma|invgaussian), shared(); stcox exactm/exactp; stcurve após ggamma/frailty, condicional a alpha=1; commit 5e07ce7). Corrigido cancelamento numérico do termo gama para theta pequeno. weights/vce(robust) com ggamma/frailty e stcurve unconditional dão erro 198. Próximo: 3f (páginas de ajuda).
- 17h55: subetapa 3f pronta e enviada (help/stsplit, strate, ct, stcrreg, stcurve .sthlp + ALIASES). Próximo: etapa 4 (testes tests/test_survival2.py, caso compat 1901, ROADMAP/docs).
- 18h10: etapas 4 e 5 prontas. PR #49 (https://github.com/yaskepalhares/opendta/pull/49) aberto, CI verde nos 3 sistemas (commit 5c19546). Próximo: etapa 6, pedir UMA rodada no Stata: no seu Mac, `cd ~/opendta && git fetch && git checkout block-survival && git pull`, depois no Stata 14 `cd "~/opendta/compat"` e `do gerar_esperados.do`; me avise quando terminar (o log novo é compat/stata/1901_survival.log). Depois: etapa 7 (comparar) e 8 (merge).

## Entry (continuation)
- Step 7 (compare) advanced: layouts of strate/stptime/stci/stmh/stmc/ltable/st data tools, stcrreg e() list and header, exactm+robust message, sttocc header done; case 1901 data changed to informative times, old expected log removed.
- Pushed to block-survival. NEXT: user runs Stata ONCE (git pull on block-survival; cd "~/opendta/compat"; do gerar_esperados.do), then compare, save log in compat/expected/, known_differences, wait CI, merge #49, #50, "How to run a block" PR, stash question, block 5 instruction (steps 8-10 pending).

## Entry (final)
- Block 4 (Survival) DONE: Stata 14 log saved in compat/expected/1901_survival.log, restricted known difference registered, CI green on 3 OS, PR #49 merged; PR #50 (cleanup) merged; PR #52 (How to run a block learnings, en + pt-BR) merged.
- Remaining: the other stream's stash (roadmap-after-1-0, stash@{0}) - user decision; block 5 (Multilevel models) in a new chat.

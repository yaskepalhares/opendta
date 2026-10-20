# Progresso: bloco 5 "Multilevel models" (docs/next-steps.md)

Arquivo local (não vai para o GitHub). Se você mandar "prossiga", o Claude lê este arquivo,
confere com o clone e continua da primeira etapa não concluída.
(Blocos anteriores: progresso-bloco1-compat.md, progresso-bloco2-datamgmt.md, progresso-bloco3-programming.md, progresso-bloco4-survival.md.)

Trabalho feito num clone na nuvem (o shell do seu Mac não alcança o GitHub); branch `block-multilevel`.
Seu clone local está em `block-survival` (já mesclado): depois do primeiro push, `git fetch && git checkout block-multilevel`.
Stata de referência: 14.0 (nunca `version 14.1` nos casos).

## Etapas

- [x] 1. Preparar: clone na nuvem a partir do main (e3fba7f), branch block-multilevel (criado e enviado ao GitHub, sem commits ainda)
- [x] 2. Explorar estimation/meglm.py e mixed.py, comandos me*/mixed, ajuda e testes; ler os manuais públicos e planejar
- [ ] 3. Implementar os comandos (core/estimation + commands + ajuda), em subetapas, cada uma enviada ao GitHub:
  - [x] 3a. motor generalizado (famílias, links, parâmetros auxiliares, ordinal): meglm, mecloglog, menbreg, meologit, meoprobit (commit 83cc32c enviado)
  - [x] 3b. meqrlogit e meqrpoisson (commit 8934069 enviado; Laplace com derivadas numéricas e modo recalculado)
  - [x] 3c. mestreg (exponential, weibull, lognormal, loglogistic e gamma; time/nohr/tratio; sem pesos, sem predict além de xb/eta/reffects) — enviado
  - [ ] 3d. efeitos aleatórios cruzados em mixed (_all: R.var)
  - [ ] 3e. páginas de ajuda (meglm, mecloglog, menbreg, meologit, meoprobit, meqrlogit, mestreg)
- [ ] 4. Testes pytest, caso compat 1902, ROADMAP fase 16 (en + pt-BR), MODULE_AREA, features (en + pt-BR), next-steps
- [ ] 5. pytest + compare.py, push e PR, CI verde
- [ ] 6. Pedir UMA rodada no Stata (git pull no branch + compat/gerar_esperados.do)
- [ ] 7. Comparar o log do 1902, corrigir, gravar em compat/expected/, known_differences (só algoritmo não documentado)
- [ ] 8. Merge do PR (gh api, merge commit)
- [ ] 9. PR pequeno em "How to run a block" (en + pt-BR), se houver aprendizado novo
- [ ] 10. Instrução do bloco 6 (Panel data)

## Feito

- 22h40: tarefa recebida. Progresso do bloco 4 arquivado em progresso-bloco4-survival.md. Pasta local: ~/opendta (acesso concedido), clone local em block-survival. Clone na nuvem em /home/claude/opendta (main e3fba7f). Li CLAUDE.md e "How to run a block". Próximo: etapa 1 (branch) e 2 (explorar).

## Decisões / achados do plano

- Próximos números livres: caso compat 1902, fase do ROADMAP 16.
- 22h50: etapa 1 pronta. Sem pytest na nuvem (PyPI bloqueado): runner próprio para testes locais; o CI roda o pytest completo. Estado atual: estimation/meglm.py só faz melogit/meprobit/mepoisson, com UM nível de efeitos aleatórios (GLMM com famílias logit/probit/poisson); mixed.py faz níveis aninhados mas recusa efeitos cruzados (`_all:`). Próximo: ler os manuais públicos (meglm, mecloglog, menbreg, meologit, meoprobit, mestreg, meqrlogit, meqrpoisson, mixed) e fechar o escopo.
- 23h20: subetapa 3a pronta e enviada (commit 83cc32c, branch block-multilevel). Conferido contra a verossimilhança integrada exatamente (60 nós): gaussiana, cloglog, menbreg (mean e constant), gamma, binomial/probit, ologit, oprobit coincidem (diferença de ll ≤ 1e-4); caso 1002 (melogit) com saída idêntica à anterior. Layouts não confirmados pelo Stata (Family:/Link:, Overdispersion:, /lnphi, var(e.y), nomes em e(b)) ficam como VERIFY para a rodada única. Próximo: 3b (meqrlogit/meqrpoisson), depois 3c mestreg e 3d cruzados.
- 23h50: subetapa 3b pronta e enviada (commit 8934069). meqrlogit == melogit nos mesmos dados; layout estilo mixed (VERIFY na rodada). Próximo: 3c mestreg (lê estimation/streg.py para reaproveitar st e as distribuições).
- 00h20: subetapa 3c pronta e enviada. mestreg conferido contra integração exata (diferença de ll ≤ 2e-4) em todas as distribuições, com e sem entrada tardia. Refatorei cmd_me em preparo + _fit_me (dataclass Fit) para o mestreg reaproveitar. Próximo: 3d efeitos cruzados em mixed (`_all: R.var`), depois 3e ajuda, 4 testes/caso 1902/docs.

- [x] 3d: efeitos cruzados em `mixed` (`_all: R.var`) prontos e enviados (ll idêntica a ML denso independente). Limitação: `predict, reffects` com R. não suportado.
- Próximo: 3e ajuda (.sthlp), 4 testes/caso 1902/docs, 5 PR+CI, 6 rodada Stata.

- [x] 3e ajuda (.sthlp) pronta; passo 4 pronto (tests/test_multilevel2.py, compat/do/1902_multilevel.do, ROADMAP fase 16, docs, roadmap_coverage). Corrigido: gradiente da quadratura não adaptativa, stddeviations, conditional(fixedonly), estat group nos comandos me.
- [x] Passo 5: PR #53 aberto (branch block-multilevel). Próximo: esperar CI verde nos 3 SOs; depois passo 6 (pedir UMA rodada Stata do caso 1902), 7 (comparar/ajustar), 8 (merge), 9 (aprendizados), 10 (instrução do bloco 6).

- [x] CI do PR #53 verde (ubuntu, macOS, Windows) no commit da branch block-multilevel.
- [ ] Passo 6 (aguardando o usuário): UMA rodada do Stata 14.0 com compat/gerar_esperados.do na branch block-multilevel; log em compat/stata/1902_multilevel.log. Ao escrever "prossiga": comparar com tools/compare.py 1902, corrigir os VERIFY, salvar em compat/expected/, known_differences, merge, aprendizados, instrução do bloco 6.

## Bloco 5 (Multilevel) — CONCLUÍDO (2026-10-20)
- Comparação com o log do Stata feita; layouts alinhados; caso 1902 registrado como diferença conhecida restrita (iterações, e(ic), datasignature; dados cruzados após bysort; meglm gaussiana com link log; laplace ~1%).
- PR #53 (bloco) e PR #54 (notas do processo em "How to run a block") mergeados com CI verde nos 3 sistemas.
- Próximo: bloco 6 (Panel data) — iniciar novo chat com a instrução fornecida.

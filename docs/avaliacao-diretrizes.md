# Avaliação de aderência às diretrizes jurídicas

Data: 8 de outubro de 2026. Estado avaliado: `main` com a fase 0 e os ícones, mais o PR #3 (fase 1a).
Diretrizes de referência: [diretrizes-juridicas.md](diretrizes-juridicas.md).

> Análise técnica à luz das diretrizes do projeto, não parecer jurídico. Nenhuma mudança foi feita a partir dela; as decisões ficaram para depois.

## Resumo

No que é mais sensível o projeto está alinhado: não há código, arquivos ou ícones do Stata, nem descompilação, e o parser, o motor e a identidade são próprios. Os maiores desvios estão na **interface** e nas **mensagens de saída**, que hoje seguem de propósito o Stata 14 o mais fielmente possível (decisão anterior às diretrizes). Também faltam licença e registro de dependências.

## Já alinhado

| Diretriz | Situação |
|---|---|
| 4.1 Código (F) | Todo o código escrito do zero; nenhum `.ado` lido ou copiado; plugin de engenharia reversa recusado no início. Política registrada em `docs/politica-clean-room.md` e `CLAUDE.md`. |
| 6 e 11 Arquitetura | Fluxo próprio: lexer → parser → AST → avaliador → saída, em camadas separadas (`lang/`, `core/`, `commands/`, `gui/`). |
| 4.4 Marca (G) | Nome próprio (OpenDTA), "Stata" não usado como marca, ícone "Tabela" próprio (verde-azulado), aviso de não afiliação na janela About. |
| 4.3 Ícones | Os 12 ícones da barra foram desenhados do zero. |
| 4.2 Documentação (E) | Texto próprio; o manual é só citado por seção, sem reprodução; do-files de teste escritos por nós. |
| 12 Testes | O Stata é referência de resultado (`compat/`), não fonte de implementação. |
| 3.3 Sintaxe (B) | Nomes de comandos, sintaxe, regras de missing e códigos `r(#)` reproduzidos como compatibilidade. |

## Divergências

1. **Interface (4.3) — alta.** As janelas Review, Results, Command, Variables e Properties têm os mesmos nomes e posições do Stata 14. Os menus seguem a ordem do Stata, e vários itens copiam o texto literal ("Summaries, tables, and tests", "Twoway graph (scatter, line, etc.)", "Binary outcomes"). A barra de ferramentas segue a sequência e os nomes do Stata ("Clear --more-- condition", "Data Browser (Browse)"). Comentários no código registram a intenção ("Janela principal no layout do Stata 14").
2. **Mensagens de erro e saídas (5) — média.** O `CLAUDE.md` manda imitar o Stata, então textos e layouts (`list`, `describe`, "(2 real changes made)", "obs was 0, now 6") são literais. Os **códigos** `r(#)`/`_rc` são interoperabilidade (B) e fazem sentido manter. O **texto** e o **layout** são expressão (D) e, pela diretriz, deveriam ser próprios. Mudar isso exige que os testes de compatibilidade passem a comparar valores, não texto.
3. **Descrição do projeto (4.4, 18) — média.** O ROADMAP diz "interface o mais próxima possível da original", e o README e o `CLAUDE.md` falam em "interface no estilo do Stata 14". A formulação "compatibilidade de sintaxe para quem migra do Stata" atende melhor.
4. **Cores do Do-file Editor (4.3) — baixa a média.** Em `gui/theme.py`, as cores foram medidas na captura das preferências do Stata. São convencionais, mas copiadas de propósito. A fonte Menlo é do sistema e não traz problema.
5. **Licença e dependências (14, 15) — pendente.** Não há `LICENSE` nem registro de dependências. O PySide6 é LGPLv3 (ou GPL), com obrigações na distribuição. NumPy e SciPy usam BSD.
6. **Formato `.dta` (8) — pendente para a fase 1b.** A ideia era um leitor próprio a partir da especificação publicada. A diretriz pede análise separada e revisão jurídica antes de distribuir. A ReadStat (MIT), uma implementação independente já existente, é alternativa de menor exposição.
7. **Uso do Stata para gerar referências (12) — consultar advogado.** Verificar se a licença do Stata restringe o uso do programa para desenvolver ou testar outro software.

## Classificação (seção 17)

| Ponto | Categoria | Situação |
|---|---|---|
| Código, parser, motor | F / A | Alinhado |
| Nomes de comandos, sintaxe, `r(#)` | B | Alinhado |
| Ícones, logo, nome | G | Alinhado |
| Layout, menus e barra no padrão do Stata | D | Desalinhado |
| Texto das mensagens e layout das saídas | B / D | Revisar |
| Descrição "o mais próxima possível" | G / E | Revisar |
| Cores do editor copiadas | D | Baixo |
| Licença e registro de dependências | — | Pendente |
| `.dta` | C | Analisar antes da fase 1b |

## Ordem sugerida, se a decisão for alinhar

1. Redesenhar a interface: layout e menus com nomes próprios.
2. Separar os códigos `r(#)` (mantidos) do texto das mensagens (próprio) e passar os testes a comparar valores.
3. Ajustar a descrição no README, ROADMAP e `CLAUDE.md`.
4. Escolher a licença e criar o registro de dependências.
5. Decidir o caminho do `.dta` antes da fase 1b.

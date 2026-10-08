# OpenDTA

Interpretador de do-files com interface no estilo do Stata 14. É um projeto pessoal de estudo, com acesso restrito.

O objetivo é que um do-file escrito para o Stata 14 rode no OpenDTA sem alteração, com a mesma sintaxe e os mesmos comandos, abreviações, regras de missing e resultados guardados. A saída na janela Results deve ser a mais próxima possível da original. O plano completo está no [ROADMAP](ROADMAP.md).

> Stata é marca registrada da StataCorp LLC. O OpenDTA não tem relação com a StataCorp e não usa código dela. Veja a [política de desenvolvimento](docs/politica-clean-room.md).

## Estado atual: fase 0 (fundação)

Já funciona:

- leitura de do-files: comentários `*`, `//` e `/* */` (aninhados), continuação `///` e `#delimit ;`;
- macros locais e globais, incluindo aninhamento (``` ``i'' ```), `` `=exp' ``, `` `++i' ``, `local ++i` e funções estendidas (`word count`, `word # of`, `display`, `list`, `piece`, `env`);
- expressões com a precedência do Stata e as regras de missing (`.`, `.a`–`.z`);
- 106 nomes de função (contando sinônimos como `trim`/`strtrim`): matemáticas, de string, lógicas, distribuições e datas;
- `display` com formatos (`%9.2f`, `%10.3e`, `%-12s`, `%td`…) e diretivas (`_col`, `_skip`, `_dup`, `_n`, `_continue`, `as`…);
- `forvalues`, `foreach`, `while`, `if`/`else if`/`else`, `continue` e `continue, break`;
- os prefixos `quietly`, `noisily` e `capture` (também em blocos `{ }`);
- `do`, `run`, `include`, `args`, `tokenize`, `confirm`, `scalar`, `macro`, `set`, `version`, `cd`, `pwd`, `exit` e `error`;
- erros com as mensagens e os códigos `r(#)` do Stata;
- interface gráfica com as janelas Results, Command, Review, Variables e Properties, além de menus e barra de ferramentas;
- console de texto e modo batch (`opendta -b do arquivo.do`, que grava `arquivo.log` como o `stata -b`).

Ainda não há dados em memória (`use`, `generate`, `summarize`…). Isso começa na fase 1.

## Instalação

Requer Python 3.11 ou mais recente.

```bash
git clone https://github.com/yaskepalhares/opendta.git
cd opendta
python -m venv .venv
# Windows: .venv\Scripts\activate     macOS/Linux: source .venv/bin/activate
pip install -e ".[dev]"
```

## Uso

```bash
opendta                      # abre a interface gráfica
opendta arquivo.do           # abre a interface e executa o do-file
opendta --console            # console de texto
opendta -b do arquivo.do     # modo batch: grava arquivo.log no diretório atual
```

## Testes

```bash
pytest                       # testes automáticos
python tools/compare.py      # compara com os logs de referência do Stata 14
```

A comparação com o Stata está descrita em [docs/como-testar.md](docs/como-testar.md).

## Organização do código

```
src/opendta/
  core/       missing values, formatos de exibição, erros r(#), saída
  lang/       lexer de do-files, macros, expressões, funções, interpretador
  commands/   comandos (um módulo por grupo)
  gui/        interface gráfica (PySide6)
  session.py  estado de uma sessão (macros, scalars, r(), e(), c())
  cli.py      linha de comando: interface, console e batch
compat/
  do/         do-files de verificação
  expected/   logs de referência gerados no Stata 14
tools/
  compare.py  roda compat/do no OpenDTA e compara com compat/expected
```

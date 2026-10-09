# OpenDTA

**An open-source interpreter for statistical analysis built around do-files and `.dta` data files.**

OpenDTA runs analysis scripts (do-files), reads and writes `.dta` datasets and comes with a desktop interface: Results, Command, Review, Variables and Properties windows, plus a data browser. The goal is a free and transparent tool for reproducible statistical work, from data management to estimation, complex survey analysis and graphics, able to run existing scripts without changes.

It is a personal study project under active development. It will be released as open source; the license is still to be defined. The full plan is in the [ROADMAP](ROADMAP.md).

*A documentação abaixo está em português.*

## Estado atual: fase 4 (Mata)

Já funciona:

- leitura de do-files: comentários `*`, `//` e `/* */` (aninhados), continuação `///` e `#delimit ;`;
- macros locais e globais, incluindo aninhamento (``` ``i'' ```), `` `=exp' ``, `` `++i' ``, `local ++i` e funções estendidas (`word count`, `word # of`, `display`, `list`, `piece`, `env`);
- expressões com a precedência do Stata e as regras de missing (`.`, `.a`–`.z`);
- 123 nomes de função (contando sinônimos como `trim`/`strtrim`): matemáticas, de string, lógicas, distribuições, datas e números aleatórios;
- `display` com formatos (`%9.2f`, `%10.3e`, `%-12s`, `%td`…) e diretivas (`_col`, `_skip`, `_dup`, `_n`, `_continue`, `as`…);
- `forvalues`, `foreach`, `while`, `if`/`else if`/`else`, `continue` e `continue, break`;
- os prefixos `quietly`, `noisily` e `capture` (também em blocos `{ }`);
- `do`, `run`, `include`, `args`, `tokenize`, `confirm`, `scalar`, `macro`, `set`, `version`, `cd`, `pwd`, `exit` e `error`;
- erros com as mensagens e os códigos `r(#)` do Stata, seguidos de uma explicação do OpenDTA: o que deu errado no comando e como resolver (desligue com `set hints off`, ou pelo menu Help → Explain Errors);
- notação científica com expoente sobrescrito na exibição (`1.0000×10¹⁰` em vez de `1.00000e+10`); `set superscript off` volta à forma do Stata, e texto que pode ser lido de volta (`string()`, macros, `file write`, exportação) sempre usa a forma do Stata;
- interface gráfica com as janelas Results, Command, Review, Variables e Properties, além de menus, barra de ferramentas com ícones próprios (claro/escuro) e Preferences (ícone do app e fonte);
- console de texto e modo batch (`opendta -b do arquivo.do`, que grava `arquivo.log` como o `stata -b`);
- programação: `program`, `syntax`, `gettoken`, `marksample`, `return`/`ereturn`, `tempvar`/`tempfile`, `preserve`/`restore`, funções estendidas de macro, `matrix`, `.ado` de terceiros no adopath, `log`/`cmdlog`, `assert`, `set trace` e `timer`;
- Viewer com `help` em SMCL (páginas próprias e as de pacotes instalados);
- dados em memória com os tipos do Stata (inclusive a precisão de `float`), varlists, `if`/`in`, `by`/`bysort`, `_n`/`_N` e subscritos;
- `set obs`, `input`, `generate`, `replace`, `drop`, `keep`, `list`, `describe`, `count`, `sort`, `gsort`, `rename`, `order`, `label`, `format`, `compress` e `recast`;
- arquivos `.dta`: leitura própria dos formatos 113 a 119 (Stata 8 em diante) e gravação nos formatos 114, 115, 117, 118 e 119, com rótulos, notas e características; formatos ainda mais antigos são lidos pela ReadStat, se o pacote `pyreadstat` estiver instalado;
- `use` (com varlist, `if` e `in`), `save`, `saveold`, `notes`, `char`, `type`, `erase`;
- texto delimitado: `import delimited`, `export delimited`, `insheet` e `outsheet`;
- planilhas: `import excel` e `export excel`;
- texto em formato livre e em colunas fixas: `infile` (inclusive com dicionário `.dct`) e `infix`; `file open/write/read/close`;
- File > Open/Save/Save as, arrastar um `.dta` para a janela, Data Editor (`browse` e `edit`, com cada edição registrada como comando) e Variables Manager (`varmanage`);
- estatística descritiva: `summarize`, `tabstat`, `tabulate` (com `chi2`, `exact`, `gamma`, `taub`…), `tab1`, `tab2`, `table`, `ci`, `ttest`, `prtest`, `correlate`, `pwcorr`, `centile`, `pctile`, `xtile`, com `fweight`, `aweight` e `iweight`;
- manipulação: `egen` (funções de grupo e de linha), `collapse`, `contract`, `expand`, `fillin`, `merge`, `append`, `joinby`, `cross`, `reshape`, `duplicates`, `isid`, `levelsof`, `recode`, `encode`/`decode`, `destring`/`tostring`, `split`, `mvencode`/`mvdecode`;
- números aleatórios: Mersenne Twister de 64 bits, `set seed` e `runiform()`, `rnormal()`, `rbinomial()`, `rpoisson()` e outras distribuições.

- Mata: blocos `mata ... end`, funções do usuário, structs e ponteiros, cerca de 150 funções da biblioteca e a ligação com os dados (`st_data`, `st_store`, `st_local`, `st_matrix`…).

Ainda não há estimação (`regress`, `logit`…), que vem na fase 5.

## Instalação (macOS)

Requer Python 3.11 ou mais recente. No Terminal:

```bash
git clone https://github.com/yaskepalhares/opendta.git
cd opendta
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

No Linux os comandos são os mesmos. No Windows, o ambiente é ativado com `.venv\Scripts\activate`.

## Uso

| Comando | O que faz |
|---|---|
| `opendta` | abre a interface gráfica |
| `opendta arquivo.do` | abre a interface e executa o do-file |
| `opendta --console` | console de texto |
| `opendta -b do arquivo.do` | modo batch: grava `arquivo.log` no diretório atual |

## Ícone

O ícone do aplicativo ("Tabela") tem versão clara (padrão) e escura, escolhida em **Edit → Preferences**. Os conjuntos completos (Windows `.ico`, macOS `.icns`, Linux, web) estão em `assets/icons/opendta/` e `assets/icons/opendta-dark/`, e as fontes SVG ficam em `assets/icons/icon-work/`.

## Testes

Os testes automáticos rodam com `pytest`. A comparação com os logs de referência do Stata 14 roda com `python tools/compare.py`. O passo a passo completo no macOS está em [docs/como-testar.md](docs/como-testar.md).

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

## Aviso

Stata é marca registrada da StataCorp LLC. O OpenDTA não tem relação com a StataCorp e não usa código dela; a compatibilidade com do-files do Stata 14 é construída só a partir de documentação pública e observação de comportamento. Veja a [política de desenvolvimento](docs/politica-clean-room.md).

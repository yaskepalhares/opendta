# Como contribuir

[English](CONTRIBUTING.md) | **Português**

Obrigado pela ajuda. O OpenDTA é um projeto pessoal de estudo, mas relatos e pull requests são bem-vindos. As issues podem ser escritas em português ou em inglês.

## Relatar um problema

- **Um resultado difere do Stata**: abra uma issue com o modelo [Difference from Stata](https://github.com/yaskepalhares/opendta/issues/new?template=difference-from-stata.md). Um do-file curto que reproduza a diferença, com a saída dos dois programas, é o relato mais útil que existe.
- **Algo quebra ou a janela se comporta mal**: use o modelo [Bug report](https://github.com/yaskepalhares/opendta/issues/new?template=bug-report.md).
- **Falta um comando de que você precisa**: confira a lista no [ROADMAP](ROADMAP.pt-BR.md#faltam-para-a-10) e use o modelo [Command request](https://github.com/yaskepalhares/opendta/issues/new?template=command-request.md). Os pedidos ajudam a decidir a ordem.

## Regras do projeto

1. **Sala limpa.** O OpenDTA é escrito a partir dos manuais do Stata, da ajuda pública, da especificação do `.dta` e do comportamento observado. Nunca copie nem adapte código ou arquivos da StataCorp (`.ado`, `.mo`, `.sthlp`, `.dlg`, `.scheme`), e nunca descompile o Stata. Veja a [política de desenvolvimento](docs/clean-room-policy.pt-BR.md).
2. **Idioma.** Código, comentários, docstrings, nomes e mensagens de commit ficam em inglês. Cada página de documentação (README, ROADMAP, CONTRIBUTING, `docs/`) tem uma cópia em português com o sufixo `.pt-BR.md` e um seletor **English | Português** no topo; quando uma página muda, atualize as duas. A interface é traduzida em `src/opendta/i18n.py`.
3. **Saída como a do Stata.** Mensagens e resultados imitam os do Stata (em inglês). Depois de uma mensagem de erro, o OpenDTA acrescenta sua própria explicação (`src/opendta/core/hints/`).
4. **Detalhes não documentados** recebem um comentário `VERIFY` no código e um caso em `compat/do/`, para que uma execução no Stata os resolva.
5. **Testes.** Antes de abrir um pull request, o `pytest` precisa passar e o `python tools/compare.py` não pode mostrar diferenças novas.
6. **Pull requests.** Trabalhe numa branch e abra um pull request para o `main`; nada vai direto para o `main`.

## Preparar o ambiente

```bash
git clone https://github.com/yaskepalhares/opendta.git
cd opendta
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
pytest
```

A organização do código, os testes e a publicação de versões estão em [docs/development.pt-BR.md](docs/development.pt-BR.md); a comparação com o Stata, em [docs/testing.pt-BR.md](docs/testing.pt-BR.md).

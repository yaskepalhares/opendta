# OpenDTA

[English](README.md) | **Português**

**Do-files do Stata 14 e dados `.dta`, num programa livre e aberto.**

O OpenDTA roda do-files escritos para o Stata 14, lê e grava bases `.dta` e tem a interface de desktop que as pessoas já conhecem: as janelas Results, Command, Review, Variables e Properties, o Data Editor, o Do-file Editor, o Viewer e a janela de gráficos. Ele cobre manipulação de dados, estatística descritiva, estimação, amostras complexas, análise de sobrevivência, modelos multinível, séries temporais, imputação múltipla, modelos de equações estruturais, Mata e gráficos. Quando um comando falha, ele explica o erro para quem está começando, olhando os dados em memória, em português ou em inglês.

![OpenDTA com uma regressão e a explicação de um erro](docs/screenshot.png)

- [Primeiros passos](docs/getting-started.pt-BR.md)
- [O que o OpenDTA faz](docs/features.pt-BR.md)
- [Roadmap e os comandos do Stata 14 que ainda faltam](ROADMAP.pt-BR.md)
- [Testar e comparar com o Stata 14](docs/testing.pt-BR.md)
- [Desenvolvimento](docs/development.pt-BR.md)
- [Política de desenvolvimento (sala limpa)](docs/clean-room-policy.pt-BR.md)

## Instalar

- **macOS** 12 ou mais recente: [Apple Silicon](https://github.com/yaskepalhares/opendta/releases/latest/download/OpenDTA-macos-apple-silicon.dmg) · [Intel](https://github.com/yaskepalhares/opendta/releases/latest/download/OpenDTA-macos-intel.dmg)
- **Windows** 10 ou 11, 64 bits: [instalador](https://github.com/yaskepalhares/opendta/releases/latest/download/OpenDTA-windows-x64-setup.exe) · [.zip portátil](https://github.com/yaskepalhares/opendta/releases/latest/download/OpenDTA-windows-x64-portable.zip)
- **Linux e outros sistemas**: [a partir do código-fonte](docs/getting-started.pt-BR.md#a-partir-do-código-fonte)

Os aplicativos não têm assinatura digital, então a primeira abertura pede uma confirmação; veja [Primeiros passos](docs/getting-started.pt-BR.md#instalar). As cópias instaladas procuram atualizações sozinhas. Todas as versões estão na página de [Releases](https://github.com/yaskepalhares/opendta/releases).

## Como contribuir

Leia o [CONTRIBUTING.pt-BR.md](CONTRIBUTING.pt-BR.md) antes de abrir uma issue ou um pull request. O relato mais útil é uma diferença em relação ao Stata 14: abra uma issue com o modelo [Difference from Stata](https://github.com/yaskepalhares/opendta/issues/new?template=difference-from-stata.md) e um do-file curto. As traduções da interface ficam em `src/opendta/i18n.py`, e o empacotamento está descrito em [docs/development.pt-BR.md](docs/development.pt-BR.md#publicar-uma-versão).

## Agradecimentos

Feito com Python, NumPy, SciPy, Matplotlib, Qt (PySide6) e openpyxl; os arquivos `.dta` antigos são lidos pelo ReadStat. Os componentes de terceiros e suas licenças estão em [packaging/THIRD_PARTY_NOTICES.md](packaging/THIRD_PARTY_NOTICES.md).

O OpenDTA é um projeto pessoal de estudo, independente e sem relação com a StataCorp. Stata é marca registrada da StataCorp LLC; o OpenDTA não usa código dela, e a compatibilidade é construída só a partir de documentação pública e do comportamento observado. Distribuído sob a [licença MIT](LICENSE).

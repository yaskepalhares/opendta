# Primeiros passos

[English](getting-started.md) | **Português**

## Instalar

Baixe o instalador do seu sistema na [versão mais recente](https://github.com/yaskepalhares/opendta/releases/latest):

| Sistema | Arquivo |
|---|---|
| macOS com Apple Silicon (M1, M2, M3, M4…) | `OpenDTA-macos-apple-silicon.dmg` |
| macOS com processador Intel | `OpenDTA-macos-intel.dmg` |
| Windows 10 ou 11, 64 bits | `OpenDTA-windows-x64-setup.exe` (instalador) ou `OpenDTA-windows-x64-portable.zip` |

**macOS** (12 ou mais recente): abra o `.dmg` e arraste o OpenDTA para Aplicativos. Na dúvida sobre o arquivo, abra  → Sobre Este Mac: "Chip Apple M…" é Apple Silicon. Como o app não tem a assinatura paga da Apple, na primeira abertura vá em **Ajustes do Sistema → Privacidade e Segurança** e clique em **Abrir Mesmo Assim** (só uma vez).

**Windows**: o instalador funciona só para o usuário atual, sem precisar de administrador, e põe o OpenDTA no menu **Abrir com** dos arquivos `.do`, `.ado` e `.dta`. Como não tem assinatura digital, o SmartScreen pode mostrar "O Windows protegeu o computador": clique em **Mais informações → Executar assim mesmo**. O `.zip` portátil roda de qualquer pasta.

**Atualizações**: o OpenDTA procura uma versão nova ao abrir (no máximo uma vez por dia; dá para desligar em **Settings**, Ctrl+, / ⌘,) e em **Help → Check for Updates**. **Install and Restart** baixa o arquivo certo, instala e abre o OpenDTA de novo.

### A partir do código-fonte

Requer Python 3.11 ou mais recente:

```bash
git clone https://github.com/yaskepalhares/opendta.git
cd opendta
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

No Windows, ative o ambiente com `.venv\Scripts\activate`. Extras opcionais: `pip install -e ".[legacy]"` lê arquivos `.dta` do Stata 7 e anteriores (ReadStat), e `pip install -e ".[xls]"` lê e grava planilhas `.xls` do Excel 97-2003.

## Executar

| Comando | O que faz |
|---|---|
| `opendta` | abre a janela |
| `opendta arquivo.do` | abre a janela e roda o do-file |
| `opendta --console` | console de texto (prompt `. `) |
| `opendta -b do arquivo.do` | modo batch: grava `arquivo.log` na pasta atual, como o `stata -b` |

## Primeiros passos

A janela tem as partes conhecidas: Results, Command, Review, Variables e Properties, além do Data Editor (`browse`, `edit`), do Do-file Editor (`doedit`), do Viewer (`help`) e da janela de gráficos. Digite os comandos na janela Command, ou abra um do-file e rode com Ctrl+D (⌘+D).

```stata
use meusdados, clear
describe
summarize idade renda
regress renda idade i.sexo
```

Quando um comando falha, o OpenDTA mostra a mensagem e o código do Stata e, logo abaixo, uma explicação para quem está começando: o que deu errado com os dados em uso, por quê, como corrigir e um link para a página de ajuda. Para desligar, use `set hints off`.

A interface pode ficar em inglês ou em português (**Settings**, ou `set language pt`). Os resultados e as mensagens no formato do Stata continuam em inglês, a menos que você peça a tradução (`set messages translated`).

Digite `help` para ver a lista de comandos, ou veja [o que o OpenDTA faz](features.pt-BR.md).

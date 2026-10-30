**English** · [Português](#português)

## What's new in 0.9.0

New:

- **`power`, part 2**: sample size, power and effect size for `power twoway` (two-way ANOVA), `power repeated` (repeated-measures ANOVA), `power cmh`, `power mcc` (matched case-control), `power trend`, `power cox`, `power exponential` and `power logrank`, with lists of values in tables. The layout and the numbers were checked against the output of a reference run.

Interface:

- **User menu**: when no item has been added yet, it says so and points to the help page that explains how to add items.
- **Help → Contents** lists every help page with its title, and **Help → Search** asks for words and lists the pages that contain all of them.
- The **Graph** button says there is no graph yet when nothing has been drawn.
- **About** describes OpenDTA as an interpreter with a Stata-style interface.

What is still missing for version 1.0 is listed in the [roadmap](https://github.com/yaskepalhares/opendta/blob/main/ROADMAP.md#missing-for-10).

## Download

| Your Mac | File |
|---|---|
| Apple Silicon (M1, M2, M3, M4…) | **OpenDTA-macos-apple-silicon.dmg** |
| Intel | **OpenDTA-macos-intel.dmg** |

| Your PC | File |
|---|---|
| Windows 10 or 11, 64-bit (Intel/AMD) | **OpenDTA-windows-x64-setup.exe** (installer) |
| Same, without installing | **OpenDTA-windows-x64-portable.zip** |

On a Mac, not sure which one? Open  → **About This Mac**: "Chip Apple M…" means Apple Silicon; "Processor Intel…" means Intel. Requires macOS 12 or later.

## Installation

### macOS

1. Open the `.dmg` and drag **OpenDTA** to the **Applications** folder.
2. On first launch, macOS says it cannot verify the developer (the app does not have Apple's paid signature). Click **OK**, go to **System Settings → Privacy & Security** and, at the bottom of the page, click **Open Anyway**. This is needed only once.

From the Terminal instead: `xattr -dr com.apple.quarantine /Applications/OpenDTA.app`

### Windows

1. Run **OpenDTA-windows-x64-setup.exe**. It installs for the current user only (no administrator rights needed) and adds OpenDTA to the **Open with** menu of `.do`, `.ado` and `.dta` files.
2. The installer is not signed, so SmartScreen may say "Windows protected your PC". Click **More info → Run anyway**.

Portable version: unzip **OpenDTA-windows-x64-portable.zip** anywhere and run `OpenDTA.exe`.

## Updates

OpenDTA 0.2.0 and later look for new releases when they open (once a day) and in **Help → Check for Updates**, and can download and install them by themselves.

## Notice

Personal study project, released under the MIT License. Stata is a registered trademark of StataCorp LLC; OpenDTA is not affiliated with StataCorp.

---

## Português

### Novidades da 0.9.0

Novo:

- **`power`, parte 2**: tamanho de amostra, poder e tamanho de efeito para `power twoway` (ANOVA de dois fatores), `power repeated` (ANOVA de medidas repetidas), `power cmh`, `power mcc` (caso-controle pareado), `power trend`, `power cox`, `power exponential` e `power logrank`, com listas de valores em tabelas. O layout e os números foram conferidos com a saída de uma execução de referência.

Interface:

- **Menu User**: enquanto nenhum item foi adicionado, ele avisa e indica a página de ajuda que explica como adicionar itens.
- **Help → Contents** lista todas as páginas de ajuda com o título, e **Help → Search** pede palavras e lista as páginas que contêm todas elas.
- O botão **Graph** avisa que ainda não há gráfico quando nada foi desenhado.
- O **About** descreve o OpenDTA como um interpretador com interface no estilo do Stata.

### Download

| Seu Mac | Arquivo |
|---|---|
| Apple Silicon (M1, M2, M3, M4…) | **OpenDTA-macos-apple-silicon.dmg** |
| Intel | **OpenDTA-macos-intel.dmg** |

| Seu PC | Arquivo |
|---|---|
| Windows 10 ou 11, 64 bits (Intel/AMD) | **OpenDTA-windows-x64-setup.exe** (instalador) |
| O mesmo, sem instalar | **OpenDTA-windows-x64-portable.zip** |

No Mac, não sabe qual é o seu? Menu  → **Sobre Este Mac**: "Chip Apple M…" é Apple Silicon; "Processador Intel…" é Intel. Requer macOS 12 ou mais recente.

### Instalação

#### macOS

1. Abra o `.dmg` e arraste o **OpenDTA** para a pasta **Aplicativos**.
2. Na primeira abertura, o macOS avisa que não pode verificar o desenvolvedor (o app não tem a assinatura paga da Apple). Clique em **OK**, vá em **Ajustes do Sistema → Privacidade e Segurança** e, no fim da página, clique em **Abrir Mesmo Assim**. Isso só é preciso uma vez.

Se preferir o Terminal: `xattr -dr com.apple.quarantine /Applications/OpenDTA.app`

#### Windows

1. Execute **OpenDTA-windows-x64-setup.exe**. Ele instala só para o usuário atual (sem precisar de administrador) e põe o OpenDTA no menu **Abrir com** dos arquivos `.do`, `.ado` e `.dta`.
2. Como o instalador não tem assinatura digital, o SmartScreen pode avisar "O Windows protegeu o computador". Clique em **Mais informações → Executar assim mesmo**.

Versão portátil: descompacte **OpenDTA-windows-x64-portable.zip** em qualquer pasta e abra `OpenDTA.exe`.

### Atualizações

O OpenDTA 0.2.0 e as versões seguintes procuram novidades ao abrir (uma vez por dia) e em **Help → Check for Updates**, e conseguem baixar e instalar sozinhos.

### Aviso

Projeto pessoal de estudo, distribuído sob a licença MIT. Stata é marca registrada da StataCorp LLC; o OpenDTA não é afiliado à StataCorp.

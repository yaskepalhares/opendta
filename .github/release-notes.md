**English** · [Português](#português)

## What's new in 0.7.0

Graphs now look like Stata 14's, each change checked against Stata's own PNG output (case 1921 in `compat/`), and the interface got a few touches:

- **Graphs**: the `kdensity` bandwidth with `aweight` and `fweight` (the number of observations with `aweight`, the sum of the weights with `fweight`); `graph bar`, `hbar` and `box` with two `over()` (grouped axis); `sort()` and `descending` in `graph dot`; narrower bars; the y axis of bars starting at zero with about four ticks; `blabel(, format())`; light boxes with filled outliers, and a `graph box` axis that follows the data; `lfit` drawn over the scatter; the `histogram` x axis starting at zero and no border between many narrow bins.
- **Interface**: each toolbar button shows its name under the icon, in English and Portuguese; the update dialog follows the Settings style; the Results tables have solid lines; do-files dropped on the window open in the Do-file Editor instead of running.

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

### Novidades da 0.7.0

Os gráficos agora ficam como os do Stata 14, cada mudança conferida com os PNGs do próprio Stata (caso 1921 em `compat/`), e a interface ganhou alguns ajustes:

- **Gráficos**: a largura de banda do `kdensity` com `aweight` e `fweight` (nº de observações com `aweight`, soma dos pesos com `fweight`); `graph bar`, `hbar` e `box` com dois `over()` (eixo agrupado); `sort()` e `descending` no `graph dot`; barras mais estreitas; eixo y das barras a partir de zero com cerca de quatro marcas; `blabel(, format())`; caixas claras com outliers preenchidos, e eixo do `graph box` que acompanha os dados; `lfit` desenhado por cima do scatter; eixo x do `histogram` começando em zero e sem borda entre muitos bins estreitos.
- **Interface**: cada botão da barra de ferramentas mostra o nome sob o ícone, em inglês e português; a janela de atualização segue o estilo de Settings; as tabelas de Results têm linhas sólidas; do-files soltos na janela abrem no Editor de do-files em vez de rodar.

O que ainda falta para a versão 1.0 está no [roadmap](https://github.com/yaskepalhares/opendta/blob/main/ROADMAP.pt-BR.md#faltam-para-a-10).

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

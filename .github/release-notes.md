**English** · [Português](#português)

## What's new in 0.11.2

Changes:

- The Do-file Editor toolbar no longer has the **Do line** and **Run** buttons; **Do** stays.

## What's new in 0.11.1

New:

- **`memory`** command, with the layout of Stata 14's report.
- `$S_DATE`, `$S_TIME` and `c(flavor)`.

Fixes:

- `log using` in a folder that does not exist now fails with r(603) (and can be caught with `capture`) instead of crashing.
- `use varlist using file` no longer marks the data as changed, so `describe` does not print "Dataset has changed since last saved".
- A `.dta` file whose sort order repeats a variable now reads as unsorted.
- `describe` puts the variable label on the next line when the value label name is long.

## What's new in 0.11.0

New:

- **Class system**: `class` definitions, `.new`, member programs and variables, `class exit` and `classutil`.
- **Multivariate statistics, part 3**: more rotation criteria, `pca, vce(normal)`, `estat` after `canon`, `mvtest` and error terms in `manova`.
- **`estat dwatson`** after `regress`.

Fixes:

- **Excel 97-2003 `.xls` files** now open in the app: `import excel`, `export excel` and `putexcel` no longer ask for an extra package to be installed.
- `.dta` files (format 117 and later) whose map stores 0 for a section are read correctly.

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

### Novidades da 0.11.2

Mudanças:

- A barra do Editor de do-files não tem mais os botões **Executar linha** e **Executar sem eco**; o **Executar** continua.

### Novidades da 0.11.1

Novo:

- Comando **`memory`**, com o layout do relatório do Stata 14.
- `$S_DATE`, `$S_TIME` e `c(flavor)`.

Correções:

- `log using` numa pasta que não existe agora falha com r(603) (e pode ser capturado com `capture`) em vez de travar.
- `use varlist using arquivo` não marca mais os dados como alterados, então o `describe` não mostra "Dataset has changed since last saved".
- Um arquivo `.dta` cuja ordenação repete uma variável agora é lido como não ordenado.
- O `describe` põe o rótulo da variável na linha seguinte quando o nome do rótulo de valor é longo.

### Novidades da 0.11.0

Novo:

- **Sistema de classes**: definições com `class`, `.new`, programas e variáveis membros, `class exit` e `classutil`.
- **Estatística multivariada, parte 3**: mais critérios de rotação, `pca, vce(normal)`, `estat` depois de `canon`, `mvtest` e termos de erro no `manova`.
- **`estat dwatson`** depois de `regress`.

Correções:

- **Planilhas `.xls` do Excel 97-2003** agora abrem no app: `import excel`, `export excel` e `putexcel` não pedem mais a instalação de um pacote extra.
- Arquivos `.dta` (formato 117 em diante) cujo mapa guarda 0 para uma seção são lidos corretamente.

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

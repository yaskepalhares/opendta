**English** · [Português](#português)

## What's new in 0.3

**If your copy says it is 0.2.0 and keeps offering 0.2.2:** versions 0.2.1 and 0.2.2 were built with the number 0.2.0 inside, so they kept offering the same update. This release fixes it: install it once and the updates work again.

New since 0.2:

- epidemiology: `cs`, `cc`, `ir` and `mcc` and their immediate forms, with Mantel-Haenszel stratification; `tabodds`, `mhodds`, `symmetry`, `dstdize` and `istdize`; interrater agreement with `kap` and `kappa`; ROC analysis with `roctab`, `roccomp`, `rocfit`, `lroc` and `lsens`;
- data management: `clonevar`, `ds`, `lookfor`, `inspect`, `labelbook`, `compare`, `cf`, `separate`, `ipolate`, `stack`, `xpose`, `expandcl`, `drawnorm`, `corr2data`, `statsby`, `putmata`, `snapshot`, `outfile`; files and folders (`copy`, `mkdir`, `dir`, `shell`, `filefilter`, `zipfile`…); example datasets (`sysuse`, `webuse`); ICD-9-CM and ICD-10-CM codes (`icd9`, `icd10`); `odbc`; business calendars (`bcal`);
- `sem` (path models, confirmatory factor analysis and latent variables, `estat gof`), `power` (means, proportions and correlations) and `irf` after `var`;
- `xtgee` and `xtreg, pa`/`be`/`mle`; `anova` with sequential sums of squares, nested terms, error terms and repeated measures; two-way `tabulate, summarize()`; `ci variances`; two-sample `prtest`; `misstable tree` and `nested`; `reshape error`; the old `merge` syntax; `export excel` to `.xls`;
- error explanations for beginners, with a diagnosis from the data and a clickable link to the help;
- `replace y = 2*L.y` is now sequential, as in Stata; `sample #, count` with # at least the number of observations no longer draws random numbers, as in Stata.

The commands still missing are listed in the [ROADMAP](https://github.com/yaskepalhares/opendta/blob/main/ROADMAP.md#missing-for-10), and the [coverage map](https://yaskepalhares.github.io/opendta/) shows the progress.

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

### Novidades da 0.3

**Se a sua cópia diz que é a 0.2.0 e fica oferecendo a 0.2.2:** as versões 0.2.1 e 0.2.2 foram geradas com o número 0.2.0 por dentro, por isso ofereciam sempre a mesma atualização. Esta versão corrige isso: instale-a uma vez e as atualizações voltam a funcionar.

Novo desde a 0.2:

- epidemiologia: `cs`, `cc`, `ir` e `mcc` e suas formas imediatas, com estratificação de Mantel-Haenszel; `tabodds`, `mhodds`, `symmetry`, `dstdize` e `istdize`; concordância entre avaliadores com `kap` e `kappa`; análise ROC com `roctab`, `roccomp`, `rocfit`, `lroc` e `lsens`;
- manipulação de dados: `clonevar`, `ds`, `lookfor`, `inspect`, `labelbook`, `compare`, `cf`, `separate`, `ipolate`, `stack`, `xpose`, `expandcl`, `drawnorm`, `corr2data`, `statsby`, `putmata`, `snapshot`, `outfile`; arquivos e pastas (`copy`, `mkdir`, `dir`, `shell`, `filefilter`, `zipfile`…); dados de exemplo (`sysuse`, `webuse`); códigos ICD-9-CM e ICD-10-CM (`icd9`, `icd10`); `odbc`; calendários de dias úteis (`bcal`);
- `sem` (modelos de caminhos, análise fatorial confirmatória e variáveis latentes, `estat gof`), `power` (médias, proporções e correlações) e `irf` depois do `var`;
- `xtgee` e `xtreg, pa`/`be`/`mle`; `anova` com somas de quadrados sequenciais, termos aninhados, termos de erro e medidas repetidas; `tabulate, summarize()` de duas entradas; `ci variances`; `prtest` de duas amostras; `misstable tree` e `nested`; `reshape error`; a sintaxe antiga do `merge`; `export excel` para `.xls`;
- explicações de erro para iniciantes, com diagnóstico a partir dos dados e link clicável para a ajuda;
- `replace y = 2*L.y` agora é sequencial, como no Stata; `sample #, count` com # maior ou igual ao número de observações não sorteia mais números aleatórios, como no Stata.

Os comandos que ainda faltam estão no [ROADMAP](https://github.com/yaskepalhares/opendta/blob/main/ROADMAP.pt-BR.md#faltam-para-a-10), e o [mapa de cobertura](https://yaskepalhares.github.io/opendta/?lang=pt) mostra o progresso.

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

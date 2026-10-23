**English** · [Português](#português)

## What's new in 0.5.0

Four new groups of commands, each checked against Stata 14 output (the remaining differences, and why, are in `compat/known_differences.txt`):

- **Panel data**: random-effects models by quadrature, `quadchk`, dynamic panels by GMM (`xtabond`, `xtdpdsys`, `xtdpd`, with `estat abond` and `estat sargan`), the linear panel estimators, `xtfrontier`, `xtunitroot` (`ht`, `breitung`, `hadri`, `fisher`), `xtline`, `xttab` and `xtdata`.
- **Multivariate statistics**: `pca`, `factor` and `rotate`, `alpha`, `cluster` (k-means, k-medians and the hierarchical linkages, with the dendrogram), `discrim` and `candisc`, `manova`, `mvreg`, `mvtest`, `canon`, `hotelling`, `ca`, `mca`, `mds`, `procrustes` and `biplot`, with their `estat` subcommands.
- **Models**: `cloglog`, `scobit`, `hetprobit`, `fracreg`, `binreg`, `glogit`; `zip`, `zinb`, `tpoisson`, `tnbreg`, `gnbreg`; `intreg`, `truncreg`, `heckman`, `heckprobit`, `heckoprobit`, `biprobit`; `ivregress liml` and `gmm`, `ivprobit`, `ivtobit`, `ivpoisson`; `constraint`, `cnsreg`, `sureg`, `reg3`, `eivreg`, `vwls`, `rreg`; `clogit`, `rologit`, `asclogit`, `slogit`, `mprobit`, `asmprobit`, `asroprobit`; `exlogistic`, `expoisson`; `nl`, `nlsur`, `gmm`, `frontier`; `teffects`, `eteffects`, `stteffects`, `etregress` and `etpoisson`.
- **Other suites**: `sem` with `mlmv`, `adf`, groups and the `estat` subcommands; `gsem`; `irt` (`1pl` to `hybrid`) with `irtgraph`; `bayesmh`, the `bayes` prefix for `regress`, `logit`, `probit` and `poisson`, `bayesstats`, `bayestest` and `bayesgraph`; new `power` methods and the `graph` option; `mi impute mvn`, `monotone`, `truncreg`, `intreg` and `nbreg`; `svy:` with sixteen more estimators.
- **Fixes**: `ivregress` now uses weights; `matrix list` shows equation names (`y:x`, `cut1:_cons`) and lists a 1×1 matrix as symmetric, and `` `:colfullnames' ``, `` `:coleq' `` and `matrix coleq`/`roweq` work with equations.

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

### Novidades da 0.5.0

Quatro novos grupos de comandos, cada um conferido com a saída do Stata 14 (as diferenças que restam, e o motivo, estão em `compat/known_differences.txt`):

- **Dados em painel**: modelos de efeitos aleatórios por quadratura, `quadchk`, painéis dinâmicos por GMM (`xtabond`, `xtdpdsys`, `xtdpd`, com `estat abond` e `estat sargan`), os estimadores lineares de painel, `xtfrontier`, `xtunitroot` (`ht`, `breitung`, `hadri`, `fisher`), `xtline`, `xttab` e `xtdata`.
- **Estatística multivariada**: `pca`, `factor` e `rotate`, `alpha`, `cluster` (k-médias, k-medianas e as ligações hierárquicas, com o dendrograma), `discrim` e `candisc`, `manova`, `mvreg`, `mvtest`, `canon`, `hotelling`, `ca`, `mca`, `mds`, `procrustes` e `biplot`, com seus subcomandos de `estat`.
- **Modelos**: `cloglog`, `scobit`, `hetprobit`, `fracreg`, `binreg`, `glogit`; `zip`, `zinb`, `tpoisson`, `tnbreg`, `gnbreg`; `intreg`, `truncreg`, `heckman`, `heckprobit`, `heckoprobit`, `biprobit`; `ivregress liml` e `gmm`, `ivprobit`, `ivtobit`, `ivpoisson`; `constraint`, `cnsreg`, `sureg`, `reg3`, `eivreg`, `vwls`, `rreg`; `clogit`, `rologit`, `asclogit`, `slogit`, `mprobit`, `asmprobit`, `asroprobit`; `exlogistic`, `expoisson`; `nl`, `nlsur`, `gmm`, `frontier`; `teffects`, `eteffects`, `stteffects`, `etregress` e `etpoisson`.
- **Outras suítes**: `sem` com `mlmv`, `adf`, grupos e os subcomandos de `estat`; `gsem`; `irt` (do `1pl` ao `hybrid`) com `irtgraph`; `bayesmh`, o prefixo `bayes` para `regress`, `logit`, `probit` e `poisson`, `bayesstats`, `bayestest` e `bayesgraph`; novos métodos do `power` e a opção `graph`; `mi impute mvn`, `monotone`, `truncreg`, `intreg` e `nbreg`; `svy:` com mais dezesseis estimadores.
- **Correções**: o `ivregress` agora usa os pesos; o `matrix list` mostra os nomes de equação (`y:x`, `cut1:_cons`) e lista uma matriz 1×1 como simétrica, e `` `:colfullnames' ``, `` `:coleq' `` e `matrix coleq`/`roweq` funcionam com equações.

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

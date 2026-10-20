**English** · [Português](#português)

## What's new in 0.4

New since 0.3:

- time series: `arima` with seasonal and multiplicative terms, `condition` and `diffuse`; `arch` (ARCH, GARCH and the other variance terms, ARCH-M, t and GED errors), `arfima`, `mgarch`, `ucm`, `sspace` and `dfactor`; `var` with `varlmar`, `varnorm`, `varstable`, `varwle`, `varbasic` and `fcast`; `vecrank`, `vec` and `svar`; `irf` with bootstrap errors, `irf ctable`, `cgraph` and `ograph`; `dfgls`, `pperron`, `prais`, the spectral and white-noise tools, `tsfilter`, `tssmooth`, `tsappend`, `tsfill`, `tsreport` and `rolling`;
- survival: `strate`, `stptime`, `stci`, `stir`, `stmh`, `stmc`, `ltable`, `stsplit`, `stjoin`, `stfill`, `stcrreg` (Fine-Gray), `stcurve`, `streg` with generalized gamma and frailty, and count-time data (`ctset`, `ct`, `cttost`, `sttoct`, `sttocc`);
- multilevel models: `meglm` and its family (`mecloglog`, `menbreg`, `meologit`, `meoprobit`, `meqrlogit`, `meqrpoisson`), `mestreg` and crossed random effects in `mixed`;
- base statistics: `jackknife`, `permute` and `simulate`; `ml` with its evaluators and `mlexp`; `mkspline`, `xi`, `nestreg`, `stepwise` and `fp`/`mfp`; `testnl`, `predictnl`, `linktest`, `contrast`, `pwcompare`, `pwmean` and `suest`; `boxcox`, `ladder`, `bitest`, `ksmirnov`, `icc`, `loneway`, `esize`, `lv`, `stem`, `lpoly` and more;
- programming: `unab` and the other variable-list tools, `mark`/`markout`, `postfile`, `matlist`, `tabdisp`, `sortpreserve`;
- epidemiology and data management: pooled and standardized `cs`/`cc`/`ir`, `rocreg`, binary `file` reading and writing, `import sas`/`spss`, and more options of `describe` and `list`.

Results of the new commands were checked against Stata 14 reference logs; the few remaining differences are listed in the repository.

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

### Novidades da 0.4

Novo desde a 0.3:

- séries temporais: `arima` com termos sazonais e multiplicativos, `condition` e `diffuse`; `arch` (ARCH, GARCH e os outros termos de variância, ARCH-M, erros t e GED), `arfima`, `mgarch`, `ucm`, `sspace` e `dfactor`; `var` com `varlmar`, `varnorm`, `varstable`, `varwle`, `varbasic` e `fcast`; `vecrank`, `vec` e `svar`; `irf` com erros bootstrap, `irf ctable`, `cgraph` e `ograph`; `dfgls`, `pperron`, `prais`, as ferramentas espectrais e de ruído branco, `tsfilter`, `tssmooth`, `tsappend`, `tsfill`, `tsreport` e `rolling`;
- sobrevivência: `strate`, `stptime`, `stci`, `stir`, `stmh`, `stmc`, `ltable`, `stsplit`, `stjoin`, `stfill`, `stcrreg` (Fine-Gray), `stcurve`, `streg` com gama generalizada e fragilidade, e dados de contagem no tempo (`ctset`, `ct`, `cttost`, `sttoct`, `sttocc`);
- modelos multinível: `meglm` e sua família (`mecloglog`, `menbreg`, `meologit`, `meoprobit`, `meqrlogit`, `meqrpoisson`), `mestreg` e efeitos aleatórios cruzados no `mixed`;
- estatística básica: `jackknife`, `permute` e `simulate`; `ml` com seus avaliadores e `mlexp`; `mkspline`, `xi`, `nestreg`, `stepwise` e `fp`/`mfp`; `testnl`, `predictnl`, `linktest`, `contrast`, `pwcompare`, `pwmean` e `suest`; `boxcox`, `ladder`, `bitest`, `ksmirnov`, `icc`, `loneway`, `esize`, `lv`, `stem`, `lpoly` e outros;
- programação: `unab` e as demais ferramentas de listas de variáveis, `mark`/`markout`, `postfile`, `matlist`, `tabdisp`, `sortpreserve`;
- epidemiologia e manipulação de dados: `cs`/`cc`/`ir` com combinação e padronização, `rocreg`, leitura e escrita binária com `file`, `import sas`/`spss` e mais opções de `describe` e `list`.

Os resultados dos comandos novos foram conferidos com logs de referência do Stata 14; as poucas diferenças restantes estão listadas no repositório.

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

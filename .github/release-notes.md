**English** · [Português](#português)

## What's new in 0.6.0

New options across nine areas, each checked against Stata 14 output (the remaining differences, and why, are in `compat/known_differences.txt`):

- **Programming**: `makecns` with constraint lists, `displaycns` and `nocnsnotes`, and `matrix dispCns`; `constraint #` without `define`; `matlist` with equation names (`showcoleq()`, `roweqonly`, `coleqonly`, `lines(eq)`…); `window` dialogs, message boxes and the User menu in the interface.
- **Data management**: `table` with a supercolumn variable and `by()`; `list` with `string()` and `display`; `import sasxport` with `describe` and `member()` and the value labels of `formats.xpf`; `export sasxport` with `rename` and `vallabfile()`.
- **Epidemiology**: `rocreg` with `probit ml` (binormal ROC model by maximum likelihood) and with `cluster()`.
- **Base statistics**: `ml` with the `lf0`–`lf2` and `gf0`–`gf2` methods, `constraints()` and `svy`; more options of `boxcox`, `icc`, `loneway`, `esize`, `pwmean`, `ameans`, `lv` and `stem`.
- **Panel data**: `xtivreg` with `first` and `regress`; `xtdata, re`; `xtrc, betas`; more options of `xtunitroot fisher` and `hadri`; `artests()` in the dynamic-panel commands.
- **Time series**: robust variances in `prais` and `arfima`; constraints in `arch`, `svar`, `sspace` and `dfactor`; `predict` after `ucm`, `sspace` and `dfactor`.
- **Multilevel models**: nested levels of random effects in `meglm` and the other `me` commands, weights at every level and robust and cluster variances; `predict` after `mestreg`.
- **Multivariate statistics**: `clustermat`, the binary similarity measures, `manovatest`, `camat`, `cabiplot`, `caprojection`, `mcaplot`, `mcaprojection`, `mdslong`, `mdsconfig`, `mdsshepard`, `procoverlay`, and `procrustes` with the unrestricted and oblique transformations.
- **Models**: `tnbreg` and `nbreg` with `dispersion(constant)`; `glogit, or`; `glm` and `binreg` with `mu()`; `offset()` in `intreg`, `truncreg`, `heckman`, `heckprobit`, `heckoprobit` and `biprobit`; `heckman` with `mills()`; `sureg` and `reg3` with `constraints()`; `ivregress gmm` with HAC weight matrices, `center` and `igmm`; `frontier` with `uhet()`, `vhet()` and `cm()`. `tpoisson` and `tnbreg` also take `ul()`, an option that Stata 14 does not have.
- **Interface**: the Review window can list the commands that a do-file ran, under its `do` line.
- **Fixes** found while checking against Stata: hierarchical cluster variables follow the dendrogram; negative binomial likelihoods no longer lose digits when the dispersion goes to zero; nested multilevel quadrature gives Stata's estimates and standard errors; and many `e()` and `r()` lists now have Stata's contents and order.

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

### Novidades da 0.6.0

Novas opções em nove áreas, cada uma conferida com a saída do Stata 14 (as diferenças que restam, e o motivo, estão em `compat/known_differences.txt`):

- **Programação**: `makecns` com listas de restrições, `displaycns` e `nocnsnotes`, e `matrix dispCns`; `constraint #` sem `define`; `matlist` com nomes de equação (`showcoleq()`, `roweqonly`, `coleqonly`, `lines(eq)`…); diálogos, caixas de mensagem e o menu User do `window` na interface.
- **Manipulação de dados**: `table` com variável de supercoluna e `by()`; `list` com `string()` e `display`; `import sasxport` com `describe` e `member()` e os rótulos de valores do `formats.xpf`; `export sasxport` com `rename` e `vallabfile()`.
- **Epidemiologia**: `rocreg` com `probit ml` (modelo ROC binormal por máxima verossimilhança) e com `cluster()`.
- **Estatística básica**: `ml` com os métodos `lf0`–`lf2` e `gf0`–`gf2`, `constraints()` e `svy`; mais opções do `boxcox`, `icc`, `loneway`, `esize`, `pwmean`, `ameans`, `lv` e `stem`.
- **Dados em painel**: `xtivreg` com `first` e `regress`; `xtdata, re`; `xtrc, betas`; mais opções do `xtunitroot fisher` e `hadri`; `artests()` nos comandos de painel dinâmico.
- **Séries temporais**: variâncias robustas no `prais` e no `arfima`; restrições no `arch`, `svar`, `sspace` e `dfactor`; `predict` depois do `ucm`, `sspace` e `dfactor`.
- **Modelos multinível**: níveis aninhados de efeitos aleatórios no `meglm` e nos outros comandos `me`, pesos em todos os níveis e variâncias robusta e por cluster; `predict` depois do `mestreg`.
- **Estatística multivariada**: `clustermat`, as medidas binárias de similaridade, `manovatest`, `camat`, `cabiplot`, `caprojection`, `mcaplot`, `mcaprojection`, `mdslong`, `mdsconfig`, `mdsshepard`, `procoverlay`, e `procrustes` com as transformações irrestrita e oblíqua.
- **Modelos**: `tnbreg` e `nbreg` com `dispersion(constant)`; `glogit, or`; `glm` e `binreg` com `mu()`; `offset()` no `intreg`, `truncreg`, `heckman`, `heckprobit`, `heckoprobit` e `biprobit`; `heckman` com `mills()`; `sureg` e `reg3` com `constraints()`; `ivregress gmm` com matrizes de pesos HAC, `center` e `igmm`; `frontier` com `uhet()`, `vhet()` e `cm()`. O `tpoisson` e o `tnbreg` também aceitam `ul()`, opção que o Stata 14 não tem.
- **Interface**: a janela Review pode listar os comandos que um do-file executou, abaixo da linha do `do`.
- **Correções** encontradas na conferência com o Stata: as variáveis da análise de cluster hierárquica seguem o dendrograma; as verossimilhanças binomiais negativas não perdem mais dígitos quando a dispersão vai a zero; a quadratura multinível aninhada dá as estimativas e os erros-padrão do Stata; e muitas listas de `e()` e `r()` agora têm o conteúdo e a ordem do Stata.

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

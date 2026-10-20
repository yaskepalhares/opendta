**English** · [Português](#português)

## What's new in 0.4.1

Interface improvements:

- **Do-file Editor**: a modern toolbar with icons (new, open, save, save as, undo and redo, find, comment/uncomment, indent and outdent), a **+** button for a new tab, and **Do**, **Run** and **Do current line** (Ctrl+Enter) on the right, as in Stata 14.
- **Color schemes** for the Results window and for the Do-file Editor, chosen separately in **Settings**: Classic Stata, High contrast, Soft, Solarized, Monokai, Dracula, One Dark, Nord, Gruvbox and GitHub.
- **Results**: `(N observations deleted)` is highlighted in red.
- **Settings** button at the right end of the toolbar, with new options: an initial folder, an automatic log of each session, and reopening the Do-file Editor's files from the last session (both off by default).
- **Confirmation before risky commands** typed in the Command window (`clear` with unsaved changes, `drop` or `keep` removing many observations or variables, `save, replace` over a much larger file), with a switch for each command in **Settings › General › Safety**.
- **Review, Variables and Properties** open in their default place every time OpenDTA starts, and go back there when chosen again in the **Window** menu.
- Portuguese interface: the standard buttons (Save / Don't Save / Cancel) and the right-click menus of text fields are now translated.

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

### Novidades da 0.4.1

Melhorias na interface:

- **Do-file Editor**: barra de botões moderna, com ícones (novo, abrir, salvar, salvar como, desfazer e refazer, localizar, comentar/descomentar, aumentar e diminuir recuo), um botão **+** para nova aba e **Do**, **Run** e **Do na linha atual** (Ctrl+Enter) à direita, como no Stata 14.
- **Esquemas de cores** para a janela Results e para o Do-file Editor, escolhidos separadamente em **Settings**: Clássico Stata, Alto contraste, Suave, Solarized, Monokai, Dracula, One Dark, Nord, Gruvbox e GitHub.
- **Results**: `(N observations deleted)` aparece em vermelho.
- Botão **Settings** na ponta direita da barra, com novas opções: pasta inicial, log automático de cada sessão e reabertura dos arquivos do Do-file Editor da última sessão (ambos desligados por padrão).
- **Confirmação antes de comandos arriscados** digitados na janela Command (`clear` com alterações não salvas, `drop` ou `keep` removendo muitas observações ou variáveis, `save, replace` sobre um arquivo bem maior), com um interruptor para cada comando em **Settings › General › Safety**.
- **Review, Variables e Properties** abrem na posição padrão a cada vez que o OpenDTA inicia e voltam para ela quando escolhidas de novo no menu **Window**.
- Interface em português: os botões padrão (Salvar / Não salvar / Cancelar) e os menus de clique direito dos campos de texto agora estão traduzidos.

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

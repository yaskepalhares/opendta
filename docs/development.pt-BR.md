# Desenvolvimento

[English](development.md) | **Português**

Como o código está organizado, como é testado e como uma versão é publicada. Para propor uma mudança, veja o [CONTRIBUTING.pt-BR.md](../CONTRIBUTING.pt-BR.md).

## Organização do código

```
src/opendta/
  core/       missing values, formatos de exibição, erros r(#), saída, explicações de erro (hints/)
  lang/       lexer de do-files, macros, expressões, funções, interpretador
  commands/   comandos (um módulo por grupo)
  io/         leitura e gravação de .dta, texto e planilhas
  estimation/ modelos, pós-estimação, svy, sem, irf e bootstrap
  graphics/   gráficos (matplotlib, esquemas próprios)
  mata/       a linguagem Mata
  gui/        interface gráfica (PySide6), diálogos e Graph Editor
  help/       páginas de ajuda (.sthlp)
  i18n.py     idioma da interface (inglês ou português)
  session.py  estado de uma sessão (macros, scalars, r(), e(), c())
  cli.py      linha de comando: interface, console e batch
compat/
  do/         do-files de verificação
  expected/   logs de referência gerados no Stata 14
packaging/    especificação do PyInstaller, script do instalador do Windows e avisos de terceiros
tools/
  compare.py         roda compat/do no OpenDTA e compara com compat/expected
  build_macos.sh     gera OpenDTA.app e o .dmg
  build_windows.ps1  gera o instalador do Windows e o .zip portátil
```

## Testes

`pip install -e ".[dev]"` instala o pytest, o ReadStat (`pyreadstat`) e o pandas. Depois:

```bash
pytest
python tools/compare.py
```

O `pytest` confere o comportamento interno e roda no GitHub a cada mudança (Linux, Windows e macOS). O `tools/compare.py` roda cada do-file de `compat/do/` e compara o log com a referência do Stata 14 em `compat/expected/`; as diferenças conhecidas estão em `compat/known_differences.txt`. O passo a passo, inclusive como gerar novas referências no Stata, está em [testing.pt-BR.md](testing.pt-BR.md).

## Ícone

O ícone do aplicativo ("Tabela") tem versão clara (padrão) e escura, escolhida em **Settings**. Os conjuntos completos (Windows `.ico`, macOS `.icns`, Linux, web) estão em `assets/icons/opendta/` e `assets/icons/opendta-dark/`, e as fontes SVG em `assets/icons/icon-work/`.

## Publicar uma versão

1. Atualize a versão em `src/opendta/__init__.py` e em `pyproject.toml` (as duas precisam coincidir com a tag).
2. Faça o merge do pull request no `main`.
3. Envie uma tag:
   ```bash
   git tag v1.0.0
   git push origin v1.0.0
   ```

O workflow [`.github/workflows/release.yml`](../.github/workflows/release.yml) confere se a tag coincide com a versão, gera os `.dmg` para Apple Silicon e Intel, o instalador do Windows e o `.zip` portátil, e cria a Release com as notas de `.github/release-notes.md`. Os links de download apontam sempre para a Release mais recente, e as cópias instaladas oferecem a atualização.

Para gerar localmente: `tools/build_macos.sh` num Mac, ou `powershell -ExecutionPolicy Bypass -File tools\build_windows.ps1` no Windows (precisa de Python 3.11+ de 64 bits e do [Inno Setup 6](https://jrsoftware.org/isinfo.php)). Em outros sistemas, `pyinstaller packaging/opendta.spec` gera a pasta `dist/OpenDTA`.

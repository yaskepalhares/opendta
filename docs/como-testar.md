# Como testar e comparar com o Stata 14 (macOS)

As instruções valem para o macOS, com o Terminal (zsh) e o Stata 14 for Mac. O repositório fica em `~/opendta`. Se estiver em outra pasta, troque o caminho nos comandos.

> **Atenção ao colar comandos no Terminal.** No zsh interativo, `#` não marca comentário: uma linha como `pytest  # roda tudo` dá erro. Os blocos abaixo não têm comentários e podem ser colados como estão.

São dois níveis de teste: os testes automáticos e a comparação com o Stata.

## 0. Preparar o ambiente (uma vez)

1. Instale o Python 3.11 ou mais recente, pelo instalador de [python.org](https://www.python.org/downloads/macos/) ou pelo Homebrew (`brew install python@3.12`).
2. No Terminal, clone o repositório e crie o ambiente virtual:
   ```bash
   cd ~
   git clone https://github.com/yaskepalhares/opendta.git
   cd opendta
   python3 -m venv .venv
   source .venv/bin/activate
   pip install -e ".[dev]"
   ```

O extra `[dev]` instala o pytest, a ReadStat (`pyreadstat`) e o pandas, usados nos testes.

## Antes de cada rodada de testes

Abra o Terminal, entre na pasta e ative o ambiente:

```bash
cd ~/opendta
source .venv/bin/activate
```

O prompt passa a mostrar `(.venv)`. Para testar um pull request antes de ele entrar na `main`, troque para o branch dele (o nome aparece no PR):

```bash
git fetch
git checkout fase-1b-arquivos
git pull
pip install -e ".[dev]"
```

Para voltar à versão principal:

```bash
git checkout main
git pull
```

Repetir o `pip install -e ".[dev]"` só é necessário quando as dependências mudam, mas não faz mal rodar sempre.

## 1. Testes automáticos (`pytest`)

Verificam o comportamento interno e rodam a cada mudança no GitHub (Linux, Windows e macOS). Para rodar no seu Mac:

```bash
pytest
```

O resultado esperado é `N passed`, sem `failed`. Se algo falhar, copie a saída do Terminal e me envie.

## 2. Abrir o programa

```bash
opendta
```

Outras formas de uso:

```bash
opendta --console
opendta -b do arquivo.do
```

A primeira abre o console de texto (prompt `. `). A segunda roda um do-file em modo batch e grava `arquivo.log` na pasta atual, como o `stata -b`.

## 3. Comparação com o Stata 14

A pasta `compat/` guarda do-files de verificação (`compat/do/`) e as saídas de referência geradas no Stata 14 (`compat/expected/`). É essa comparação que define se o OpenDTA está fiel.

### Gerar as referências no Stata (a cada lote novo de do-files)

1. Abra o Stata 14.
2. Na janela Command, vá até a pasta `compat` do repositório e rode o gerador:
   ```stata
   cd "~/opendta/compat"
   do gerar_esperados.do
   ```
3. Os arquivos `compat/expected/*.log` aparecem. O gerador força `set dp period` durante a execução e restaura a sua configuração no fim.
4. Envie os logs: anexe os arquivos novos aqui na conversa, ou faça commit pelo Terminal:
   ```bash
   cd ~/opendta
   git add compat/expected
   git commit -m "Referências do Stata 14"
   git push
   ```

### Comparar

Com o ambiente ativado (`source .venv/bin/activate`), rode a comparação de todos os do-files:

```bash
python tools/compare.py
```

Para mostrar as linhas diferentes:

```bash
python tools/compare.py --show-diff
```

Para comparar só um arquivo, pelo começo do nome:

```bash
python tools/compare.py 0104
```

A saída indica, para cada do-file, se as saídas são iguais (`ok`), diferentes (`XX`) ou se ainda falta a referência (`?`). Cabeçalhos de log, linhas em branco, espaços no fim das linhas e o carimbo de data do `describe` são ignorados.

### Arquivos `.dta` (fase 1b)

O leitor e o gravador de `.dta` são próprios. Nos testes automáticos, cada arquivo gravado é relido pela ReadStat, uma implementação independente da especificação. O juiz final continua sendo o Stata.

1. **Stata → OpenDTA.** No Stata, grave um arquivo com rótulos, notas e missing estendidos:
   ```stata
   sysuse auto, clear
   notes: teste de notas
   replace rep78 = .a in 1
   save "~/Desktop/auto14.dta", replace
   saveold "~/Desktop/auto13.dta", replace
   ```
   No OpenDTA, abra o arquivo e confira:
   ```stata
   use "~/Desktop/auto14.dta", clear
   describe
   notes
   list in 1/5
   ```
   Repita com `auto13.dta`. Também dá para abrir pelo menu File > Open ou arrastando o arquivo para a janela.
2. **OpenDTA → Stata.** No OpenDTA:
   ```stata
   save "~/Desktop/do_opendta.dta", replace
   saveold "~/Desktop/do_opendta13.dta", replace
   ```
   No Stata, abra cada um e rode `describe`, `notes`, `char list` e `list`. Qualquer aviso do Stata ao abrir o arquivo é uma divergência.

Os do-files `compat/do/0104_arquivos.do` a `0107_infile_infix.do` criam e apagam os próprios arquivos (`.dta`, CSV, Excel e texto) na pasta atual, com nomes começando por `odta_`.

## Reportar uma divergência encontrada no uso

Abra uma *issue* com o modelo **Divergência com o Stata**: um do-file mínimo, a saída do Stata 14 e a saída do OpenDTA. Cada divergência confirmada vira um caso em `compat/do/`, para não voltar a acontecer.

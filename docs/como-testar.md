# Como testar e comparar com o Stata 14

São dois níveis de teste.

## 1. Testes automáticos (`pytest`)

Verificam o comportamento interno e rodam a cada mudança no GitHub. Para rodar localmente:

```bash
pytest
```

## 2. Comparação com o Stata 14

A pasta `compat/` guarda do-files de verificação (`compat/do/`) e as saídas de referência geradas no Stata 14 (`compat/expected/`). É essa comparação que define se o OpenDTA está fiel.

### Gerar as referências no Stata (uma vez por lote novo de do-files)

1. Abra o Stata 14.
2. Vá até a pasta `compat` do repositório:
   ```stata
   cd "C:/caminho/para/opendta/compat"
   ```
3. Rode:
   ```stata
   do gerar_esperados.do
   ```
4. Os arquivos `compat/expected/*.log` aparecem. Faça commit deles, ou me envie.

### Comparar

```bash
python tools/compare.py               # todos
python tools/compare.py 0002          # só um arquivo
python tools/compare.py --show-diff   # mostra as linhas diferentes
```

A saída indica, para cada do-file, se as saídas são iguais (`ok`), diferentes (`XX`) ou se ainda falta a referência (`?`). Cabeçalhos de log, linhas em branco e espaços no fim das linhas são ignorados.

### Reportar uma divergência encontrada no uso

Abra uma *issue* com o modelo **Divergência com o Stata**: um do-file mínimo, a saída do Stata 14 e a saída do OpenDTA. Cada divergência confirmada vira um caso em `compat/do/`, para não voltar a acontecer.

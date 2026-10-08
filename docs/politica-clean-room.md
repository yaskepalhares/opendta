# Política de desenvolvimento ("clean room")

O OpenDTA reimplementa o comportamento do Stata a partir de fontes públicas. Não usa código da StataCorp.

## Pode ser usado como fonte

- Os manuais do Stata (PDFs), a ajuda (`help`) e os FAQs públicos
- A especificação pública do formato `.dta`
- Artigos, livros e materiais didáticos sobre o Stata e sobre os métodos estatísticos
- A observação do comportamento: rodar um do-file no Stata 14 e comparar a saída com a do OpenDTA
- Bibliotecas abertas com licença compatível (por exemplo, ReadStat para `.dta`)

## Não pode ser usado

- Arquivos `.ado`, `.mo`, `.mlib`, `.sthlp`, `.dlg` ou `.scheme` distribuídos com o Stata, nem para copiar trechos nem como base de adaptação
- Descompilação, desmontagem ou inspeção dos executáveis e bibliotecas do Stata
- Ícones, logotipos e outros elementos gráficos do Stata

## Na prática

- Cada comando é escrito a partir do manual. Quando o manual não define um detalhe (por exemplo, quando o `%g` troca a notação fixa pela exponencial), o código recebe a marca `VERIFICAR` e um caso em `compat/do/`. O comportamento real é fixado pelo log gerado no Stata.
- Pacotes do SSC e do-files de usuários podem ser usados como casos de teste, respeitando suas licenças.
- A interface reproduz o layout e o fluxo de trabalho do Stata 14, mas os ícones e a identidade visual são próprios.

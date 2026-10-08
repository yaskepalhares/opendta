# Instruções para sessões de desenvolvimento

- Projeto: OpenDTA, um interpretador de do-files com interface no estilo do Stata 14. A meta está em ROADMAP.md.
- Siga docs/politica-clean-room.md: nada de código, arquivos ou elementos gráficos da StataCorp, nem descompilação. Use só manuais, ajuda pública, especificação do .dta e observação de comportamento.
- Documentação, comentários e mensagens de commit em português. Nomes de código em inglês. Mensagens de erro e saídas imitam o Stata (em inglês).
- Detalhes não documentados recebem a marca `VERIFICAR` no código e um caso em compat/do/.
- Antes de abrir PR: `pytest` passando e `python tools/compare.py` sem regressões.
- Trabalhe em branch e abra PR para revisão. Não faça push direto na main.

"""Interface language (phase 9).

Two languages: English (default, as in Stata) and Portuguese. The choice applies to
menus, windows, dialogs, OpenDTA's error explanations (core/hints.py)
and OpenDTA's own messages. Messages in Stata's format (errors
such as "variable x not found") stay in English unless the user asks
for the translation: those who compare with Stata logs and manuals need them identical.

    tr("Open...")            -> "Abrir..." with Portuguese on
    L("en text", "pt text")  -> picks by language (texts with variables)
    translate_error(msg)     -> translated error message, if requested
"""

from __future__ import annotations

import re

_LANG = "en"
_MESSAGES = "stata"      # "stata" (English, as in Stata) or "translated"

LANGUAGES = {"en": "English", "pt": "Português"}


def set_language(lang: str) -> None:
    global _LANG
    _LANG = "pt" if str(lang).lower().startswith("pt") else "en"


def language() -> str:
    return _LANG


def set_messages(mode: str) -> None:
    global _MESSAGES
    _MESSAGES = "translated" if str(mode).lower().startswith("t") else "stata"


def messages() -> str:
    return _MESSAGES


def L(en: str, pt: str) -> str:
    return pt if _LANG == "pt" else en


def tr(text: str) -> str:
    """Translation of a fixed interface text (the text itself if there is none)."""
    if _LANG != "pt" or not text:
        return text
    if text in PT:
        return PT[text]
    # "&File" and the like: try without the accelerator
    bare = text.replace("&", "")
    if bare in PT:
        t = PT[bare]
        return "&" + t if text.startswith("&") else t
    return text


# ---------------------------------------------------------------------------
# interface texts
# ---------------------------------------------------------------------------

PT: dict[str, str] = {
    # menus principais
    "File": "Arquivo", "Edit": "Editar", "Data": "Dados", "Graphics": "Gráficos",
    "Statistics": "Estatística", "User": "Usuário", "Window": "Janela", "Help": "Ajuda",
    # File
    "Open...": "Abrir...", "Open": "Abrir", "Save": "Salvar", "Save as...": "Salvar como...",
    "Save as": "Salvar como", "Do...": "Executar do-file...", "Do": "Executar",
    "Change working directory...": "Mudar a pasta de trabalho...",
    "Change working directory": "Mudar a pasta de trabalho", "Log": "Log", "Begin...": "Iniciar...",
    "Begin log": "Iniciar log", "Close": "Fechar", "Suspend": "Suspender", "Resume": "Retomar",
    "View...": "Ver...", "View": "Ver", "Print...": "Imprimir...", "Print": "Imprimir", "Exit": "Sair",
    # Edit
    "Copy": "Copiar", "Copy table": "Copiar tabela", "Find...": "Localizar...", "Find": "Localizar",
    "Clear Results": "Limpar Results", "Preferences...": "Preferências...", "Preferences": "Preferências",
    "Find in Results:": "Localizar em Results:",
    "Input requested by display:": "Entrada pedida pelo display:",
    # Data
    "Describe data": "Descrever os dados", "Data Editor": "Editor de dados",
    "Create or change data": "Criar ou mudar dados", "Create new variable": "Criar variável",
    "Change contents of variable": "Mudar o conteúdo da variável",
    "Create new variable (extended)": "Criar variável (estendido)", "Recode": "Recodificar",
    "Encode string": "Codificar texto", "Destring": "Texto para número", "Variables": "Variáveis",
    "Rename": "Renomear", "Label variable": "Rotular variável", "Define value label": "Definir rótulo de valores",
    "Keep variables": "Manter variáveis", "Drop variables": "Apagar variáveis", "Order": "Ordenar colunas",
    "Sort": "Ordenar", "Sort (ascending/descending)": "Ordenar (crescente/decrescente)",
    "Data utilities": "Utilitários de dados", "Duplicates report": "Relatório de duplicatas",
    "Codebook": "Codebook", "Missing values": "Valores ausentes", "Compress": "Compactar", "Count": "Contar",
    "Combine datasets": "Combinar bases", "Merge": "Juntar (merge)", "Append": "Empilhar (append)",
    "Reshape": "Reformatar (reshape)", "Collapse": "Resumir (collapse)",
    # Graphics
    "Twoway graph (scatter, line, etc.)": "Gráfico twoway (dispersão, linha etc.)", "Line plot": "Gráfico de linha",
    "Bar chart": "Gráfico de barras", "Dot chart": "Gráfico de pontos", "Pie chart": "Gráfico de setores",
    "Histogram": "Histograma", "Box plot": "Box plot", "Kernel density": "Densidade kernel",
    "Scatterplot matrix": "Matriz de dispersão", "Combine graphs": "Combinar gráficos",
    "Export graph": "Exportar gráfico", "Manage graphs": "Gerenciar gráficos",
    "Margins plot": "Gráfico de margens", "Residuals vs. fitted": "Resíduos x ajustados",
    "Normal quantile plot": "Gráfico de quantis normais", "Lowess smoother": "Suavização lowess",
    # Statistics
    "Summaries, tables, and tests": "Resumos, tabelas e testes", "Summary statistics": "Estatísticas descritivas",
    "One-way table": "Tabela simples", "Two-way table with chi2": "Tabela cruzada com qui-quadrado",
    "Table of statistics": "Tabela de estatísticas", "Correlations": "Correlações", "t test": "Teste t",
    "Proportion test": "Teste de proporção", "Confidence intervals": "Intervalos de confiança",
    "Centiles": "Percentis", "Nonparametric tests": "Testes não paramétricos",
    "Normality tests": "Testes de normalidade", "ANOVA": "ANOVA",
    "Linear models and related": "Modelos lineares e afins", "Linear regression": "Regressão linear",
    "Instrumental variables": "Variáveis instrumentais", "Absorbed indicators": "Indicadores absorvidos",
    "Quantile regression": "Regressão quantílica", "Tobit": "Tobit", "Binary outcomes": "Desfechos binários",
    "Logistic (odds ratios)": "Logística (razões de chances)", "Logit": "Logit", "Probit": "Probit",
    "Ordinal outcomes": "Desfechos ordinais", "Ordered logit": "Logit ordinal", "Ordered probit": "Probit ordinal",
    "Categorical outcomes": "Desfechos categóricos", "Multinomial logit": "Logit multinomial",
    "Count outcomes": "Desfechos de contagem", "Poisson": "Poisson", "Negative binomial": "Binomial negativa",
    "Generalized linear models": "Modelos lineares generalizados", "GLM": "GLM",
    "Longitudinal/panel data": "Dados longitudinais/painel", "Declare panel": "Declarar painel",
    "Fixed effects": "Efeitos fixos", "Random effects": "Efeitos aleatórios", "Panel summary": "Resumo do painel",
    "Time series": "Séries temporais", "Declare time variable": "Declarar variável de tempo",
    "ARIMA models": "Modelos ARIMA", "Vector autoregression": "Autorregressão vetorial (VAR)",
    "Lag-order selection": "Seleção do número de defasagens",
    "Dickey-Fuller unit-root test": "Teste de raiz unitária de Dickey-Fuller", "Correlogram": "Correlograma",
    "Newey-West regression": "Regressão com erros de Newey-West",
    "Check for Updates...": "Procurar atualizações...", "Check for Updates": "Procurar atualizações",
    "Update Available": "Atualização disponível", "Skip This Version": "Pular esta versão",
    "Remind Me Later": "Lembrar depois", "Install and Restart": "Instalar e reiniciar",
    "Open Download Page": "Abrir a página de download", "Update": "Atualização", "Updates": "Atualizações",
    "Check for updates when OpenDTA opens (once a day)": "Procurar atualizações ao abrir o OpenDTA (uma vez por dia)",
    "Settings": "Configurações",
    "Settings...": "Configurações...",
    "Search": "Buscar",
    "General": "Geral",
    "Menus, windows, dialogs and OpenDTA's own messages.": "Menus, janelas, caixas de diálogo e as mensagens próprias do OpenDTA.",
    "Translate Stata-style messages": "Traduzir as mensagens no estilo do Stata",
    "Error messages in Portuguese. Keep them in English to compare with Stata logs and manuals.": "Mensagens de erro em português. Deixe em inglês para comparar com os logs e os manuais do Stata.",
    "Errors": "Erros",
    "Explain errors": "Explicar os erros",
    "After an error message, explain what went wrong and suggest a fix. Same as set hints on.": "Depois da mensagem de erro, explica o que deu errado e sugere como corrigir. O mesmo que set hints on.",
    "Appearance": "Aparência",
    "Theme": "Tema",
    "Light, dark, or the same as the operating system.": "Claro, escuro ou o mesmo do sistema operacional.",
    "Light": "Claro",
    "Follow the system": "Seguir o sistema",
    "The icon in the Dock or taskbar and in the windows.": "O ícone no Dock ou na barra de tarefas e nas janelas.",
    "Font": "Fonte",
    "Font family": "Família da fonte",
    "Results, Command and Review windows.": "Janelas Results, Command e Review.",
    "Font size": "Tamanho da fonte",
    "Numbers": "Números",
    "Superscript exponents": "Expoentes sobrescritos",
    "Show scientific notation as 1.0000×10¹⁰ instead of 1.00000e+10. Text that can be read back (exports, macros) always uses Stata's form.": "Mostra a notação científica como 1.0000×10¹⁰ em vez de 1.00000e+10. Textos que podem ser lidos de volta (exportações, macros) usam sempre a forma do Stata.",
    "Stata": "Stata",
    "Precise": "Precisa",
    "Arithmetic": "Aritmética",
    "Stata reproduces Stata's rounding. Precise gives the same estimates with more exact arithmetic (set numerics).": "Stata reproduz os arredondamentos do Stata. Precisa dá as mesmas estimativas com aritmética mais exata (set numerics).",
    "Check automatically": "Verificar automaticamente",
    "Look for a new version on GitHub when OpenDTA opens, at most once a day.": "Procura uma versão nova no GitHub quando o OpenDTA abre, no máximo uma vez por dia.",
    "Check Now": "Verificar agora",
    "Version": "Versão",
    "About": "Sobre",
    "Open-source interpreter for do-files and .dta data, compatible with Stata 14.": "Interpretador de código aberto para do-files e dados .dta, compatível com o Stata 14.",
    "Project": "Projeto",
    "Source code, issues and releases.": "Código-fonte, problemas relatados e versões.",
    "Report a problem": "Relatar um problema",
    "Differences from Stata and bugs.": "Diferenças em relação ao Stata e erros.",
    "Release notes": "Notas das versões",
    "What changed in each version.": "O que mudou em cada versão.",
    "License": "Licença",
    "Free to use, copy, modify and distribute. Stata is a registered trademark of StataCorp LLC; OpenDTA is not affiliated with StataCorp.": "Livre para usar, copiar, modificar e distribuir. Stata é marca registrada da StataCorp LLC; o OpenDTA não é afiliado à StataCorp.",
    "Reset": "Redefinir",
    "Restore Defaults": "Restaurar padrões",
    "Restore default settings": "Restaurar as configurações padrão",
    "Theme, icon, font, language and the other settings on these pages.": "Tema, ícone, fonte, idioma e as demais configurações destas páginas.",
    "Restore all settings to their defaults?": "Restaurar todas as configurações para o padrão?",
    "Multiple imputation": "Imputação múltipla", "Declare data (mlong)": "Declarar os dados (mlong)",
    "Register imputed variables": "Registrar as variáveis a imputar",
    "Impute by linear regression": "Imputar por regressão linear", "Chained equations": "Equações encadeadas",
    "Estimate and combine": "Estimar e combinar",
    "Multilevel mixed-effects models": "Modelos multinível de efeitos mistos",
    "Linear regression (mixed)": "Regressão linear (mixed)", "Logistic regression": "Regressão logística",
    "Probit regression": "Regressão probit", "Poisson regression": "Regressão de Poisson",
    "Intraclass correlation": "Correlação intraclasse",
    "Survival analysis": "Análise de sobrevivência", "Declare survival data": "Declarar dados de sobrevivência",
    "Describe": "Descrever", "Summarize": "Resumir", "List survivor function": "Listar a função de sobrevivência",
    "Kaplan-Meier graph": "Gráfico de Kaplan-Meier", "Compare survivor functions": "Comparar funções de sobrevivência",
    "Cox model": "Modelo de Cox", "Parametric models": "Modelos paramétricos",
    "Test proportional hazards": "Testar riscos proporcionais",
    "Failure variable:": "Variável de falha (evento):",
    "Subject ID variable:": "Variável de identificação do sujeito:",
    "Origin (time expression):": "Origem (expressão de tempo):", "Scale:": "Escala:",
    "Efron method for ties": "Método de Efron para empates",
    "Report coefficients, not hazard ratios": "Mostrar coeficientes, não razões de risco",
    "Strata variables:": "Variáveis de estrato:", "Distribution:": "Distribuição:",
    "Accelerated failure-time metric": "Métrica de tempo de falha acelerado",
    "Failure function (1 - S)": "Função de falha (1 - S)",
    "Confidence interval": "Intervalo de confiança", "Test:": "Teste:",
    "Exponential": "Exponencial", "Lognormal": "Lognormal", "Loglogistic": "Loglogística",
    "Survey data analysis": "Amostras complexas (survey)", "Declare design": "Declarar o desenho",
    "Describe design": "Descrever o desenho", "Means": "Médias", "Proportions": "Proporções",
    "Tables": "Tabelas", "Design effects": "Efeitos de desenho",
    "Postestimation": "Pós-estimação", "Predictions": "Predições", "Marginal effects": "Efeitos marginais",
    "Wald tests": "Testes de Wald", "Linear combinations": "Combinações lineares",
    "Nonlinear combinations": "Combinações não lineares", "Compare estimates": "Comparar estimativas",
    "Likelihood-ratio test": "Teste da razão de verossimilhança", "Hausman test": "Teste de Hausman",
    "Bootstrap": "Bootstrap",
    # Window / Help
    "Command": "Comando", "Results": "Resultados", "Graph": "Gráfico", "Viewer": "Visualizador",
    "Do-file Editor": "Editor de do-files", "Review": "Histórico", "Properties": "Propriedades",
    "Search...": "Buscar...", "Contents": "Conteúdo", "Explain Errors (set hints)": "Explicar erros (set hints)",
    "About Error Explanations": "Sobre as explicações de erro",
    "Superscript Exponents (set superscript)": "Expoentes sobrescritos (set superscript)",
    "About OpenDTA": "Sobre o OpenDTA", "Error explanations": "Explicações de erro",
    # toolbar and panels
    "Log (begin or close)": "Log (iniciar ou fechar)", "Data Editor (Edit)": "Editor de dados (editar)",
    "Data Browser (Browse)": "Navegador de dados (ver)", "Variables Manager": "Gerenciador de variáveis",
    "Clear --more-- condition": "Limpar a condição --more--", "Break": "Interromper", "Toolbar": "Barra de ferramentas",
    "Filter commands here": "Filtrar comandos aqui", "Filter variables here": "Filtrar variáveis aqui",
    "Name": "Nome", "Label": "Rótulo", "Type": "Tipo", "Format": "Formato", "Value label": "Rótulo de valores",
    "Notes": "Notas", "Filename": "Arquivo", "Observations": "Observações", "Size": "Tamanho",
    "Memory": "Memória", "Sorted by": "Ordenado por", "Property": "Propriedade", "Value": "Valor",
    "Data in memory have changed.": "Os dados na memória foram alterados.",
    "Do you want to save the changes before opening another dataset?":
        "Deseja salvar as alterações antes de abrir outra base?",
    "Stata data (*.dta);;All files (*)": "Dados Stata (*.dta);;Todos os arquivos (*)",
    "Stata data (*.dta)": "Dados Stata (*.dta)",
    "Do-files (*.do *.ado);;All files (*)": "Do-files (*.do *.ado);;Todos os arquivos (*)",
    "Logs and help (*.smcl *.log *.sthlp *.txt);;All files (*)":
        "Logs e ajuda (*.smcl *.log *.sthlp *.txt);;Todos os arquivos (*)",
    "SMCL log (*.smcl);;Text log (*.log)": "Log SMCL (*.smcl);;Log de texto (*.log)",
    'Write "{cmd}" in the Command window': 'Escreve "{cmd}" na janela Command',
    "Open a dialog to build the command": "Abre uma caixa de diálogo para montar o comando",
    # Preferences
    "Interface appearance": "Aparência da interface", "Light (default)": "Claro (padrão)", "Dark": "Escuro",
    "Tables": "Tabelas", "Solid table lines": "Linhas contínuas nas tabelas",
    "Draws the lines of the tables in the Results window as continuous lines. Off: dashes, as in Stata "
    "and in the logs.":
        "Desenha as linhas das tabelas na janela Results de forma contínua. Desligado: tracejadas, como "
        "no Stata e nos logs.",
    "Application icon": "Ícone do aplicativo", "Font of the Results, Command and Review windows":
        "Fonte das janelas Results, Command e Review", "Family": "Família",
    "Default on this system": "Padrão neste sistema", "Language": "Idioma",
    "Interface language": "Idioma da interface",
    "Also translate Stata-style messages (errors)": "Traduzir também as mensagens no formato do Stata (erros)",
    "Keep them in English to compare with Stata logs and manuals.":
        "Deixe em inglês para comparar com os logs e os manuais do Stata.",
    "Table — light (default)": "Tabela — claro (padrão)", "Table — dark": "Tabela — escuro",
    # Do-file Editor
    "New": "Novo", "Run": "Executar sem eco", "Next": "Próximo", "Replace": "Substituir",
    "Replace all": "Substituir tudo", "Find:": "Localizar:", "Replace:": "Substituir:",
    "Execute the selection, or the whole file (do)": "Executa a seleção, ou o arquivo inteiro (do)",
    "Execute the selection, or the whole file, without showing commands (run)":
        "Executa a seleção, ou o arquivo inteiro, sem mostrar os comandos (run)",
    "Do-files (*.do);;Ado-files (*.ado)": "Do-files (*.do);;Ado-files (*.ado)",
    "Untitled": "Sem título",
    "Session": "Sessão", "Safety": "Segurança", "Initial folder": "Pasta inicial", "Choose...": "Escolher...",
    "Default": "Padrão", "Not set": "Não definida", "Color scheme": "Esquema de cores",
    "The working directory OpenDTA starts in. Takes effect the next time OpenDTA opens.":
        "A pasta de trabalho em que o OpenDTA começa. Vale na próxima vez que o OpenDTA abrir.",
    "Save a log of each session": "Salvar um log de cada sessão",
    "Show the commands of do-files in Review": "Mostrar os comandos dos do-files no Review",
    "After do, run or the Do button of the Do-file Editor, Review lists the commands the do-files ran, under "
    "the do line (click the arrow to open them; commands with errors in red). Off: one line per do, as in Stata.":
        "Depois de do, run ou do botão Do do Editor de do-files, o Review lista os comandos que os do-files "
        "executaram, sob a linha do do (clique na seta para abri-los; comandos com erro em vermelho). Desligado: "
        "uma linha por do, como no Stata.",
    "selection": "seleção", "line": "linha",
    "Starts a text log when OpenDTA opens, in a logs folder inside the initial folder (or your home folder). "
    "Takes effect the next time OpenDTA opens.":
        "Inicia um log em texto ao abrir o OpenDTA, numa pasta logs dentro da pasta inicial (ou da sua pasta "
        "pessoal). Vale na próxima vez que o OpenDTA abrir.",
    "Confirm clear": "Confirmar clear", "Confirm drop": "Confirmar drop", "Confirm keep": "Confirmar keep",
    "Confirm save, replace": "Confirmar save, replace",
    "In the Command window, ask before clear when the data in memory have unsaved changes.":
        "Na janela Command, pergunta antes do clear quando os dados na memória têm alterações não salvas.",
    "In the Command window, ask before drop removes many observations or most of the variables.":
        "Na janela Command, pergunta antes do drop remover muitas observações ou a maior parte das variáveis.",
    "In the Command window, ask before keep removes many observations or most of the variables.":
        "Na janela Command, pergunta antes do keep remover muitas observações ou a maior parte das variáveis.",
    "In the Command window, ask before save, replace overwrites a file much larger than the data in memory.":
        "Na janela Command, pergunta antes do save, replace sobrescrever um arquivo bem maior que os dados "
        "na memória.",
    "Colors of the Results window and the Review list: the classic look, high contrast, a soft one, or a "
    "scheme common in code editors.":
        "Cores da janela Results e da lista Review: o visual clássico, alto contraste, um suave, ou um esquema "
        "comum em editores de código.",
    "Syntax colors and background of the Do-file Editor, independent of the Results window.":
        "Cores de sintaxe e fundo do Editor de do-files, independentes da janela Results.",
    "Reopen last session's files": "Reabrir os arquivos da última sessão",
    "When the Do-file Editor opens, reopen the files that were open when it was last used.":
        "Ao abrir o Editor de do-files, reabre os arquivos que estavam abertos na última vez que foi usado.",
    "Confirm command": "Confirmar comando", "Run anyway": "Executar mesmo assim",
    "Cancel and open Settings": "Cancelar e abrir Configurações",
    "You can turn these confirmations off in Settings > General > Safety.":
        "Você pode desligar estas confirmações em Configurações > Geral > Segurança.",
    "New tab": "Nova aba",
    "Undo": "Desfazer", "Redo": "Refazer", "Comment/uncomment": "Comentar/descomentar",
    "Indent": "Aumentar recuo", "Outdent": "Diminuir recuo", "Do line": "Executar linha",
    "Execute the current line (do)": "Executa a linha atual (do)",
    # Graph window / editor
    "Save...": "Salvar...", "Export...": "Exportar...", "Save graph": "Salvar gráfico",
    "OpenDTA graph (*.gph)": "Gráfico do OpenDTA (*.gph)", "Edit graph": "Editar gráfico",
    "Graph Editor": "Editor de gráfico", "Title": "Título", "Subtitle": "Subtítulo", "Note": "Nota",
    "Caption": "Legenda inferior", "X axis title": "Título do eixo x", "Y axis title": "Título do eixo y",
    "Show legend": "Mostrar legenda", "Scheme": "Esquema", "Apply": "Aplicar", "Width (in)": "Largura (pol.)",
    "Height (in)": "Altura (pol.)", "Grid lines on the y axis": "Linhas de grade no eixo y",
    "Log scale (x)": "Escala logarítmica (x)", "Log scale (y)": "Escala logarítmica (y)",
    "Revert": "Desfazer edições",
    # Data Browser / Variables Manager / Viewer
    "Browse": "Ver", "Mode": "Modo", "Back": "Voltar", "Forward": "Avançar",
    "help topic, or a file to view": "tópico de ajuda, ou um arquivo para ver",
    # command dialogs
    "OK": "OK", "Submit": "Enviar", "Cancel": "Cancelar", "Copy to Command window": "Copiar para a janela Command",
    "Variables:": "Variáveis:", "Dependent variable:": "Variável dependente:",
    "Independent variables:": "Variáveis independentes:", "Variable:": "Variável:",
    "Row variable:": "Variável da linha:", "Column variable:": "Variável da coluna:",
    "If (expression):": "If (expressão):", "In (range):": "In (intervalo):",
    "Weights:": "Pesos:", "Options:": "Opções:", "Additional options:": "Outras opções:",
    "Group variable:": "Variável de grupo:", "Over:": "Por grupo (over):", "Statistic:": "Estatística:",
    "Detail": "Detalhado", "Chi-squared test": "Teste qui-quadrado", "Row percentages": "Porcentagens na linha",
    "Column percentages": "Porcentagens na coluna", "Robust standard errors": "Erros-padrão robustos",
    "Cluster variable:": "Variável de cluster:", "Odds ratios": "Razões de chances",
    "Incidence-rate ratios": "Razões de taxas de incidência", "Confidence level:": "Nível de confiança:",
    "Number of bins:": "Número de barras:", "Scale:": "Escala:", "Normal density": "Densidade normal",
    "Y variable:": "Variável y:", "X variable:": "Variável x:", "Add linear fit": "Acrescentar ajuste linear",
    "Graph title:": "Título do gráfico:", "PSU variable:": "Variável de UPA:", "Strata variable:": "Variável de estrato:",
    "Sampling weight:": "Peso amostral:", "FPC variable:": "Variável de FPC:",
    "Variance estimation:": "Estimação da variância:", "Subpopulation (variable):": "Subpopulação (variável):",
    "Estimator:": "Estimador:", "Panel variable:": "Variável de painel:", "Time variable:": "Variável de tempo:",
    "Model:": "Modelo:", "Family:": "Família:", "Link:": "Ligação:", "Exponentiated coefficients":
        "Coeficientes exponenciados", "Two-sample test by group:": "Teste de duas amostras por grupo:",
    "Hypothesized mean:": "Média sob H0:", "Unequal variances": "Variâncias diferentes",
    "Dialog": "Diálogo", "The command will be:": "O comando será:",
    "(default)": "(padrão)", "(missing)": "(missing)", "Base outcome:": "Categoria de referência:",
    "Between effects": "Efeitos entre grupos", "Maximum likelihood": "Máxima verossimilhança",
    "Binomial": "Binomial", "Gaussian": "Gaussiana", "Inverse Gaussian": "Gaussiana inversa",
    "Gamma": "Gama", "Identity": "Identidade", "Complementary log-log": "Complemento log-log",
    "Power -1": "Potência -1", "Bonferroni comparisons": "Comparações de Bonferroni",
    "Covariances instead of correlations": "Covariâncias em vez de correlações",
    "Display numeric codes instead of value labels": "Mostrar os códigos em vez dos rótulos",
    "Treat missing values like other values": "Tratar valores ausentes como os demais",
    "Exposure variable:": "Variável de exposição:", "Fisher's exact test": "Teste exato de Fisher",
    "Density": "Densidade", "Fraction": "Fração", "Frequency": "Frequência", "Percent": "Porcentagem",
    "Hypothesized proportion:": "Proporção sob H0:", "Response variable:": "Variável resposta:",
    "Linearized": "Linearização", "Jackknife": "Jackknife", "BRR": "BRR", "SDR": "SDR",
    "Certainty": "Certeza (certainty)", "Scaled": "Escalonado (scaled)", "Centered": "Centrado (centered)",
    "Strata with a single PSU:": "Estratos com uma só UPA:",
    "Percentages instead of proportions": "Porcentagens em vez de proporções",
    "Quantile:": "Quantil:", "Relative-risk ratios": "Razões de risco relativo",
    "Standardized beta coefficients": "Coeficientes beta padronizados", "Summary table": "Tabela-resumo",
    "Data Editor (Browse)": "Editor de dados (ver)",
    "Save changes to {name}?": "Salvar as alterações em {name}?",
    "No effect: OpenDTA always runs with set more off":
        "Sem efeito: o OpenDTA roda sempre com set more off",
    "Not available yet: commands run in the interface thread":
        "Ainda não disponível: os comandos rodam na linha de execução da interface",
    "After an error message, explain what went wrong and how to fix it. Same as typing "
    "\"set hints on\" or \"set hints off\".":
        "Depois da mensagem de erro, explica o que deu errado e como corrigir. O mesmo que "
        "digitar \"set hints on\" ou \"set hints off\".",
    "Show scientific notation as 1.0000\u00d710\u00b9\u2070 instead of 1.00000e+10. Same as typing "
    "\"set superscript on\" or \"set superscript off\".":
        "Mostra a notação científica como 1.0000\u00d710\u00b9\u2070 em vez de 1.00000e+10. O mesmo "
        "que digitar \"set superscript on\" ou \"set superscript off\".",
}


# ---------------------------------------------------------------------------
# error messages in Stata's format (only when translation is requested)
# ---------------------------------------------------------------------------

_ERRORS: list[tuple[str, str]] = [
    (r"variable (\S+) not found", r"variável \1 não encontrada"),
    (r"variable (\S+) already defined", r"a variável \1 já existe"),
    (r"(\S+) ambiguous abbreviation", r"\1: abreviação ambígua"),
    (r"(\S+) not found", r"\1 não encontrado"),
    (r"command (\S+) is unrecognized", r"comando \1 não reconhecido"),
    (r"option (\S+) not allowed", r"opção \1 não permitida"),
    (r"unknown function (\S+)\(\)", r"função desconhecida \1()"),
    (r"invalid syntax", "sintaxe inválida"),
    (r"type mismatch", "tipos incompatíveis"),
    (r"varlist required", "lista de variáveis obrigatória"),
    (r"varlist not allowed", "lista de variáveis não permitida"),
    (r"too few variables specified", "poucas variáveis especificadas"),
    (r"too many variables specified", "variáveis demais especificadas"),
    (r"no observations", "nenhuma observação"),
    (r"no variables defined", "nenhuma variável definida"),
    (r"Obs\. nos\. out of range", "números de observação fora do intervalo"),
    (r"no; data in memory would be lost", "não; os dados na memória seriam perdidos"),
    (r"file (.+) not found", r"arquivo \1 não encontrado"),
    (r"file (.+) already exists", r"o arquivo \1 já existe"),
    (r"file (.+) could not be opened", r"não foi possível abrir o arquivo \1"),
    (r"assertion is false", "a afirmação é falsa"),
    (r"last estimates not found", "últimas estimativas não encontradas"),
    (r"not sorted", "não ordenado"),
    (r"= ?exp(ression)? required", "expressão (=exp) obrigatória"),
    (r"unmatched quote", "aspas sem par"),
    (r"invalid numlist", "numlist inválida"),
    (r"(\S+) invalid name", r"\1: nome inválido"),
    (r"data not set up for svy, use svyset", "dados sem desenho amostral; use svyset"),
    (r"time variable not set.*", "variável de tempo não declarada; use tsset"),
    (r"panel variable not set.*", "variável de painel não declarada; use xtset"),
    (r"observations must be at least (\d+)", r"o número de observações deve ser pelo menos \1"),
    (r"conformability error", "erro de conformidade (dimensões)"),
    (r"matrix not positive definite", "matriz não positiva definida"),
    (r"insufficient observations", "observações insuficientes"),
    (r"graph (\S+) not found", r"gráfico \1 não encontrado"),
    (r"graph (\S+) already exists", r"o gráfico \1 já existe"),
    (r"no graphs in memory", "nenhum gráfico na memória"),
    (r"convergence not achieved", "convergência não alcançada"),
    (r"outcome does not vary", "o desfecho não varia"),
]


def translate_error(msg: str) -> str:
    if _LANG != "pt" or _MESSAGES != "translated" or not msg:
        return msg
    lines = msg.split("\n")
    out = []
    for ln in lines:
        for pat, rep in _ERRORS:
            if re.fullmatch(pat, ln):
                ln = re.sub(pat, rep, ln)
                break
        out.append(ln)
    return "\n".join(out)

# Instruções do projeto: software estatístico compatível com Stata

## 1. Objetivo

Desenvolver um software estatístico open source, independente do Stata, com interface própria e implementação própria das funcionalidades estatísticas.

O principal objetivo é oferecer uma experiência familiar para usuários habituados ao Stata, permitindo, quando tecnicamente e juridicamente apropriado, o uso de comandos com sintaxe compatível ou semelhante à do Stata.

A intenção **não é criar uma cópia ou fork do Stata**, mas desenvolver uma implementação independente de um ambiente estatístico baseado em funcionalidades e metodologias estatísticas conhecidas.

---

## 2. Princípio jurídico e técnico fundamental

O projeto deve seguir uma abordagem de **implementação independente (clean-room / clean implementation)**.

A regra geral é:

> Reproduzir funcionalidades e comportamentos necessários à interoperabilidade, mas não reproduzir código, documentação, assets, identidade visual ou outros elementos protegidos do Stata.

Sempre que houver dúvida entre:

- implementar uma funcionalidade de forma independente; ou
- reproduzir diretamente algum elemento do Stata;

preferir a primeira alternativa.

Este documento não constitui parecer jurídico. Questões jurídicas específicas devem ser submetidas posteriormente a um advogado especializado em propriedade intelectual e software.

---

# 3. O que pode ser reproduzido

## 3.1 Funcionalidades estatísticas

É permitido como objetivo do projeto implementar funcionalidades estatísticas equivalentes às existentes em softwares estatísticos consolidados.

Exemplos:

- regressão linear;
- regressão logística;
- regressão de Poisson;
- estatística descritiva;
- tabelas de frequência;
- testes de hipótese;
- análise de sobrevivência;
- modelos multiníveis;
- métodos para amostras complexas;
- pesos amostrais;
- estratificação;
- conglomerados;
- estimativas de variância;
- marginal effects;
- padronização;
- análise de tendência temporal.

A implementação deve ser própria.

---

## 3.2 Algoritmos e métodos estatísticos

Métodos estatísticos conhecidos podem ser implementados independentemente.

Por exemplo:

```text
OLS
maximum likelihood
logistic regression
Poisson regression
Taylor linearization
bootstrap
jackknife
```

O algoritmo/método estatístico não deve ser confundido com a implementação proprietária específica de um software.

Sempre que possível, implementar os métodos a partir de:

- literatura científica;
- livros-texto;
- documentação acadêmica;
- especificações matemáticas;
- padrões estatísticos;
- outras fontes independentes.

---

## 3.3 Sintaxe compatível

O software pode ter como objetivo reconhecer comandos semelhantes ou compatíveis com Stata.

Exemplo:

```stata
regress renda idade sexo
```

poderia ser interpretado pelo nosso parser como:

```text
regress
    ↓
regressão linear
    ↓
dependent variable = renda
independent variables = idade, sexo
```

Da mesma forma, podemos buscar compatibilidade com comandos como:

```text
summarize
tabulate
regress
logit
poisson
generate
replace
drop
keep
sort
merge
reshape
collapse
svy
margins
```

A compatibilidade de sintaxe deve ser tratada como uma camada de interoperabilidade.

---

# 4. O que NÃO deve ser copiado

## 4.1 Código-fonte do Stata

Nunca:

- copiar código-fonte;
- copiar arquivos `.ado`;
- copiar componentes internos;
- descompilar o software para reutilizar código;
- modificar código proprietário e incorporá-lo ao projeto;
- copiar implementação interna de algoritmos.

Mesmo quando a funcionalidade desejada for idêntica, a implementação deve ser escrita do zero.

---

## 4.2 Documentação

Não copiar:

- textos do manual;
- descrições de comandos;
- exemplos;
- tabelas;
- figuras;
- explicações;
- notas técnicas;
- páginas da documentação;
- estrutura extensa de conteúdo reproduzida sistematicamente.

Podemos documentar nossas próprias funcionalidades, utilizando linguagem própria.

Exemplo aceitável:

```text
regress estimates the coefficients of a linear regression
model using ordinary least squares.
```

A documentação deve explicar como **nosso software** funciona.

---

## 4.3 Interface visual

A interface deve ser própria.

É aceitável possuir elementos convencionais de softwares estatísticos, como:

```text
File
Edit
Data
Statistics
Help
```

ou:

- editor de comandos;
- painel de variáveis;
- console;
- janela de resultados;
- gráficos;
- tabelas;
- explorador de dados.

Entretanto, evitar reprodução deliberada da:

- disposição específica da interface do Stata;
- identidade visual;
- ícones;
- imagens;
- elementos gráficos;
- aparência característica;
- combinação de elementos destinada a fazer o programa parecer o Stata.

---

## 4.4 Marca e identidade

Não utilizar:

- logo do Stata;
- elementos de identidade visual da StataCorp;
- nome "Stata" como nome do produto;
- branding que possa sugerir que o projeto é oficial ou afiliado à StataCorp.

O projeto deve possuir:

- nome próprio;
- logo próprio;
- identidade visual própria;
- documentação própria.

Quando necessário mencionar Stata, utilizar referências descritivas e factuais, especialmente para indicar compatibilidade ou migração.

Exemplo:

> "Command syntax compatibility is provided for users migrating from Stata."

Evitar afirmações que possam sugerir endosso, associação ou origem comum.

---

# 5. Mensagens de erro

Não há necessidade de reproduzir literalmente as mensagens de erro do Stata.

Por exemplo, em vez de copiar:

```text
variable x not found
r(111);
```

podemos utilizar:

```text
Error: variable 'x' was not found.
```

ou:

```text
Unknown variable: x
```

A finalidade funcional da mensagem pode ser equivalente, mas a implementação e apresentação devem ser próprias.

---

# 6. Comandos

O objetivo é permitir que usuários habituados ao Stata tenham uma curva de aprendizado mínima.

Exemplo:

```stata
use dados.dta

summarize idade

regress renda idade sexo

logit diabetes idade sexo
```

Idealmente, o software deverá conseguir interpretar comandos dessa natureza.

Entretanto, a implementação deverá seguir:

```text
User command
      ↓
Parser próprio
      ↓
AST / representação interna
      ↓
Statistical engine próprio
      ↓
Resultados
      ↓
Output/UI própria
```

Nunca implementar comandos simplesmente copiando a implementação interna do Stata.

---

# 7. Exemplo: `svy`

Um caso prioritário para usuários de epidemiologia é o suporte a amostras complexas.

Por exemplo:

```stata
svyset setor [pweight=peso], strata(estrato)

svy: regress renda idade sexo
```

O software pode implementar independentemente:

```text
svyset
    ↓
configuração do desenho amostral

svy:
    ↓
execução do modelo considerando o desenho amostral

peso
estrato
conglomerado/UPA
    ↓
estimativas e variâncias apropriadas
```

A metodologia pode ser baseada em literatura estatística e epidemiológica, sem copiar a implementação do Stata.

---

# 8. Compatibilidade com arquivos `.dta`

Investigar separadamente a possibilidade de:

- ler arquivos `.dta`;
- escrever arquivos `.dta`;
- preservar labels;
- preservar value labels;
- preservar formatos;
- preservar missing values;
- preservar metadados.

A implementação do suporte ao formato deve ser independente.

Não incorporar componentes proprietários do Stata.

Antes de distribuir oficialmente suporte completo a `.dta`, realizar revisão jurídica específica sobre:

- especificação do formato;
- engenharia reversa;
- interoperabilidade;
- eventuais direitos sobre componentes do formato;
- licença e documentação utilizada para implementar o parser.

---

# 9. Classificação de risco para decisões de desenvolvimento

Utilizar a seguinte classificação preliminar:

| Elemento | Tratamento |
|---|---|
| Implementar regressão linear | OK |
| Implementar OLS | OK |
| Implementar logística | OK |
| Implementar `svy` | OK, implementação própria |
| Implementar `margins` | OK, implementação própria |
| Usar metodologia estatística conhecida | OK |
| Criar parser próprio | OK |
| Criar sintaxe compatível | Geralmente OK, revisar casos específicos |
| Usar comandos com nomes iguais | Avaliar caso a caso |
| Ler formato `.dta` | Avaliar separadamente |
| Reproduzir resultados estatísticos | Objetivo de compatibilidade |
| Copiar código `.ado` | NÃO |
| Copiar código do Stata | NÃO |
| Descompilar Stata para reutilização | NÃO |
| Copiar manual | NÃO |
| Copiar exemplos da documentação | Evitar |
| Copiar ícones | NÃO |
| Copiar logo | NÃO |
| Copiar identidade visual | NÃO |
| Fazer interface deliberadamente semelhante | Evitar |
| Utilizar o nome Stata como marca do projeto | NÃO |

Esta tabela é uma orientação de desenvolvimento, não uma conclusão jurídica definitiva.

---

# 10. Princípio de implementação

Sempre que houver uma funcionalidade que também exista no Stata, perguntar:

### Pergunta 1
A funcionalidade é uma ideia, método estatístico, algoritmo ou comportamento geral?

Se sim:

> Implementar independentemente.

### Pergunta 2
Estamos copiando código ou uma implementação proprietária específica?

Se sim:

> Parar e implementar novamente do zero.

### Pergunta 3
Estamos copiando texto, imagens, ícones, documentação ou identidade visual?

Se sim:

> Não utilizar.

### Pergunta 4
A compatibilidade exige determinado comportamento para que scripts possam ser migrados?

Se sim:

> Implementar o comportamento de forma independente e documentar nossa própria implementação.

---

# 11. Arquitetura recomendada

O projeto deve separar claramente:

```text
                    ┌─────────────────────┐
                    │       GUI           │
                    │     própria         │
                    └──────────┬──────────┘
                               │
                    ┌──────────▼──────────┐
                    │   Command Parser    │
                    │      próprio        │
                    └──────────┬──────────┘
                               │
                    ┌──────────▼──────────┐
                    │ Internal AST / IR   │
                    └──────────┬──────────┘
                               │
             ┌─────────────────┼─────────────────┐
             │                 │                 │
       Data Engine       Statistical Engine   Graphics
             │                 │                 │
             └─────────────────┼─────────────────┘
                               │
                    ┌──────────▼──────────┐
                    │    Output Engine    │
                    │       próprio       │
                    └─────────────────────┘
```

Essa separação é importante tanto para a arquitetura do software quanto para manter uma implementação claramente independente.

---

# 12. Testes de compatibilidade

Podemos comparar nosso software com o Stata para verificar se os resultados são equivalentes.

Exemplo:

```text
Input:
regress y x1 x2
```

Executar independentemente:

```text
Stata → resultados de referência
Nosso software → resultados próprios
```

Comparar:

- coeficientes;
- erros-padrão;
- estatísticas de teste;
- p-values;
- intervalos de confiança;
- graus de liberdade;
- resultados de previsão;
- resultados marginais.

A comparação de resultados deve ser usada como **teste de compatibilidade**, não como mecanismo para copiar a implementação.

Sempre que possível, utilizar também resultados matematicamente conhecidos ou implementações independentes como referência.

---

# 13. Testes e validação estatística

Cada procedimento estatístico deve possuir:

1. testes unitários;
2. testes numéricos;
3. casos extremos;
4. comparação com resultados publicados;
5. comparação com referências matemáticas;
6. comparação com outros softwares independentes;
7. testes de regressão para evitar alterações acidentais.

Para métodos epidemiológicos, priorizar validação em:

- dados reais anonimizados;
- dados simulados;
- exemplos de artigos científicos;
- datasets públicos;
- cenários de amostras complexas.

---

# 14. Licenciamento do projeto

O projeto deverá utilizar uma licença open source apropriada.

Avaliar opções como:

```text
MIT
Apache-2.0
GPLv3
AGPLv3
```

A escolha deve considerar:

- possibilidade de uso comercial;
- obrigação ou não de disponibilizar modificações;
- bibliotecas utilizadas;
- compatibilidade entre licenças;
- possibilidade de distribuição proprietária de versões derivadas;
- objetivo de construir um ecossistema de extensões.

Não escolher a licença apenas por preferência. Fazer uma análise das dependências antes.

---

# 15. Dependências

Toda biblioteca externa deve ser registrada.

Para cada dependência, verificar:

```text
Nome
Versão
Licença
URL oficial
Uso no projeto
Compatibilidade da licença
```

Nunca incorporar código de terceiros sem verificar sua licença.

---

# 16. Política de desenvolvimento

O projeto deve manter uma separação clara entre:

### Referência

```text
Stata
```

e:

### Implementação

```text
Nosso software
```

O Stata pode ser utilizado como referência de:

- sintaxe;
- comportamento esperado;
- interoperabilidade;
- experiência de usuários;
- terminologia necessária para compatibilidade.

Mas o código, documentação, interface e identidade do projeto devem ser desenvolvidos independentemente.

---

# 17. Regra prática para o Claude

Durante o desenvolvimento, sempre que uma solicitação envolver algo existente no Stata, classificar primeiro em uma destas categorias:

```text
A. Funcionalidade estatística
B. Sintaxe/interoperabilidade
C. Formato de dados
D. Interface/UX
E. Documentação
F. Código/implementação
G. Marca/identidade
```

Aplicar as seguintes regras:

```text
A → implementar independentemente
B → buscar compatibilidade, evitando copiar expressão protegida
C → analisar separadamente
D → criar design próprio
E → escrever documentação própria
F → código 100% próprio
G → não reproduzir
```

Se houver dúvida jurídica relevante, **não assumir que algo é permitido**. Sinalizar a questão e propor uma alternativa tecnicamente equivalente e de menor risco.

---

# 18. Objetivo final

O produto final deve poder ser descrito como:

> Um software estatístico open source, independente, com interface própria e motor estatístico próprio, que busca oferecer compatibilidade de sintaxe e interoperabilidade para usuários familiarizados com ambientes estatísticos baseados em comandos.

O projeto deve evitar ser descrito ou desenvolvido como:

> "um clone do Stata"

ou:

> "uma cópia open source do Stata".

A distinção não é apenas de nomenclatura. Ela deve estar refletida na arquitetura, no código, na documentação, na interface e no processo de desenvolvimento.
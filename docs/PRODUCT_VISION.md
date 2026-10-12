# Delta Chemistry e Grimperium: problema, mecanismo e visão do produto

> **Nomes.** **Delta Chemistry** é o hub de rotinas de cálculo de propriedades (a visão descrita abaixo). **Grimperium** é o módulo de entalpia de formação no estado padrão, o único implementado hoje. Nas seções seguintes, “Grimperium” designa o que existe no código; “hub” e “Delta Chemistry”, a visão.

## O problema, sem pressupor conhecimento de química computacional

Para calcular uma propriedade termodinâmica de uma molécula, existem métodos computacionais com diferentes custos e graus de precisão. Na prática, muitas vezes é necessário escolher entre **obter uma estimativa rapidamente, com maior erro**, e **esperar muito mais tempo por um resultado mais preciso**.

O Grimperium investiga uma alternativa a essa escolha: **usar resultados de métodos mais precisos para ensinar um modelo de machine learning a corrigir os resultados de métodos mais rápidos**. A meta é reduzir o erro mantendo parte importante da vantagem de custo computacional do método rápido; isso é um objetivo a ser validado em cada domínio, não uma garantia universal.

## Como a alternativa funciona

1. **Escolher métodos comparáveis:** um método rápido e outro de referência, mais preciso, para a **mesma propriedade** e dentro de um domínio químico definido.
2. **Construir dados em pares:** executar ou reunir resultados rastreáveis dos dois métodos para moléculas compatíveis. Essa é uma etapa cara; por isso a automação do pipeline importa.
3. **Treinar a correção:** calcular a diferença entre referência e método rápido; o modelo de machine learning aprende a estimar essa diferença usando características da molécula e os dados pertinentes ao método.
4. **Aplicar:** para uma nova molécula, executar o método rápido e adicionar a correção prevista pelo modelo compatível.
5. **Verificar se compensou:** medir erros antes/depois da correção em dados separados do treinamento, inspecionar falhas e comparar o custo/tempo do novo percurso.

```text
delta de treinamento = valor_de_referência − valor_do_método_rápido
valor_corrigido = valor_do_método_rápido + delta_previsto_pelo_modelo
```

Essa estratégia é o **delta-learning**. Ela *não* transforma o método rápido em um método preciso por si só, nem dispensa a verificação científica da correção.

## Por que um hub?

A finalidade maior do Grimperium não é oferecer apenas um script de cálculo isolado. A **visão** é um hub para reunir **diferentes propriedades, métodos rápidos, referências mais precisas e modelos de correção**, com interfaces consistentes para configurar experimentos, gerar dados, treinar, comparar e usar os métodos corrigidos.

O hub precisa conservar a identidade da propriedade e do método, compatibilidade entre modelos, unidades, domínio de aplicação, proveniência dos dados e evidências de validação. Não se deve misturar alvos de treinamento de propriedades diferentes ou aplicar um modelo fora do domínio para o qual foi criado.

## O que existe no repositório hoje

**O único tipo de propriedade suportado atualmente no catálogo é a entalpia padrão de formação (`standard_enthalpy_of_formation`, ΔHf°).** Isso está explícito em `src/grimperium/calculation/methods/registry.py`.

- Exemplo de par para treinamento: cálculo semiempírico **PM7** (rápido) e referência **CBS-QB3** (mais precisa), com diferenças aprendidas por ML.
- O projeto possui pipelines científicos, carregamento/tratamento de dados, treinamento e avaliação, definições de métodos e interfaces de acompanhamento. A existência de um modelo treinado ou de uma definição de método **não significa que todo fluxo de execução já aplique ML automaticamente**.
- Há um componente independente, **Semi-Imperium**, para calcular/catalogar entalpia padrão de formação com AM1, PM3 e PM7. Ele **não** é um pipeline de correção por ML. Ver [Semi-Imperium](SEMI_IMPERIUM.md).
- A aplicação de dashboard permite inspecionar o progresso dos cálculos e comparações de dados; detalhes em [dashboard local](../apps/dashboard/README.md).

## O que ainda é visão

Expandir o catálogo para outras propriedades e novos pares método rápido/referência, treinar correções específicas e avaliar cada combinação com critérios científicos próprios. **Não anunciar suporte a uma propriedade nova até existir implementação, testes e validação correspondentes.**

## Como descrever o projeto para outras pessoas

> O Grimperium investiga uma alternativa à escolha entre cálculos rápidos e pouco precisos ou cálculos precisos e demorados. Em vez de abandonar os métodos rápidos, usa dados de métodos mais precisos para treinar modelos de Machine Learning que aprendem a corrigir as estimativas dos métodos rápidos. O projeto automatiza a geração e o tratamento desses dados e, hoje, aplica essa abordagem à entalpia padrão de formação. A visão é expandir a ideia para um hub de diferentes propriedades e métodos.

Essa é a descrição **do produto e de sua visão**. Em currículo, a descrição pode ser adaptada à vaga, mas não deve inverter o mecanismo: **o ML aprende correções a partir de cálculos mais precisos; não apenas “integra ferramentas” nem substitui o método rápido**.

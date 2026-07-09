# Módulo 0 — Fundamentos e Desenho Experimental

## Conceitos-chave

### O que é um transcriptoma

O transcriptoma é o conjunto completo de moléculas de RNA transcritas a partir do genoma de um organismo, em um tecido, estágio de desenvolvimento ou condição experimental específicos, em um dado momento. Diferente do genoma, que é (praticamente) estático, o transcriptoma é dinâmico: muda entre tecidos, entre condições (estresse, infecção, tratamento) e ao longo do tempo. RNA-Seq é a técnica de sequenciamento de alto rendimento usada para capturar esse conjunto de transcritos, convertendo RNA em cDNA e sequenciando-o em plataformas como Illumina (short-read), PacBio Iso-Seq ou Oxford Nanopore (ONT, long-read).

### O que é um organismo não-modelo

Organismo-modelo é aquele com genoma de referência publicado, bem anotado e mantido pela comunidade (ex.: *Drosophila melanogaster*, *Arabidopsis thaliana*, camundongo). Organismo não-modelo é qualquer espécie sem esse recurso, o que descreve a esmagadora maioria da biodiversidade eucariótica, incluindo praticamente todos os insetos de interesse agronômico fora de meia dúzia de espécies-praga extensivamente estudadas. Os dois organismos-caso deste curso, *Anticarsia gemmatalis* (lagarta-da-soja, Lepidoptera) e *Mahanarva spectabilis* (cigarrinha-das-pastagens, Hemiptera), são exemplos típicos: relevância econômica e ecológica claras, mas nenhum genoma de referência de qualidade disponível publicamente no momento de escrita deste material.

### Genome-guided vs. montagem de novo

Existem dois caminhos para reconstruir transcritos a partir de reads de RNA-Seq:

- **Genome-guided (referência-guiada):** os reads são alinhados a um genoma de referência (com um alinhador ciente de splicing, como STAR ou HISAT2), e os transcritos são reconstruídos a partir desse alinhamento (ex.: StringTie, Cufflinks). É o caminho preferido sempre que existe uma referência de qualidade razoável, mesmo que de uma espécie próxima.
- **De novo:** os transcritos são reconstruídos diretamente a partir dos reads, sem nenhuma referência, usando grafos de de Bruijn (Trinity, rnaSPAdes) ou, no caso de long-read, algoritmos de clustering e consenso (RNA-Bloom2, IsoSeq3/pbccs+Isoseq, isONform). É o caminho obrigatório quando não existe genoma utilizável.

### Principais desafios de montagens de novo

Montagem de novo é substancialmente mais difícil que montagem genômica de novo ou montagem genome-guided, por razões específicas do RNA:

1. **Heterozigosidade:** indivíduos diploides (ou pools de indivíduos) carregam variantes alélicas que o assemblador pode interpretar erroneamente como transcritos distintos, inflando o número de "genes" reconstruídos.
2. **Isoformas (splicing alternativo):** um mesmo gene pode gerar múltiplos transcritos com sobreposição parcial de sequência, o que confunde os grafos de de Bruijn baseados em k-mer. O grafo não "sabe" se uma bifurcação é uma variante alélica, um erro de sequenciamento ou uma isoforma real.
3. **Redundância:** o mesmo motivo acima gera contigs quase idênticos e fragmentados no output bruto do assemblador, exigindo uma etapa posterior de clustering/redução de redundância (CD-HIT, por exemplo) antes de a montagem ser usada como referência.
4. **Ausência de referência para QC:** sem genoma, não é possível calcular métricas como taxa de mapeamento contra o genoma verdadeiro ou completude por comparação direta. O controle de qualidade precisa recorrer a métricas indiretas: BUSCO (completude via genes ortólogos de cópia única esperados no clado), N50, ExN50, taxa de remapeamento dos próprios reads contra a montagem (Bowtie2/Salmon) e métricas de contiguidade/redundância (Transrate, DETONATE).

Esses quatro desafios explicam por que a montagem de novo é tratada, no fluxo deste curso, como um processo iterativo de "montar → avaliar → filtrar/reduzir redundância → reavaliar", em vez de um passo único.

### Por que montagem de novo? Árvore de decisão

Antes de decidir por uma montagem de novo, vale percorrer explicitamente os critérios abaixo, nessa ordem:

1. **Existe genoma de referência publicado para a espécie-alvo?**
   - Sim, com boa qualidade (contiguidade e anotação adequadas) → use pipeline genome-guided (STAR/HISAT2 + StringTie). Não monte de novo desnecessariamente: você perde precisão ao ignorar informação genômica disponível.
   - Não → vá para o critério 2.

2. **Existe genoma de referência de uma espécie filogeneticamente próxima (mesmo gênero ou família), com sintenia/identidade esperada alta?**
   - Sim, e o objetivo é apenas quantificação de expressão gênica em genes conservados → um pipeline genome-guided cross-species pode funcionar como aproximação, mas com perda de sensibilidade para genes específicos da linhagem e para eventos de splicing alternativo divergentes. Avalie criticamente se essa perda é aceitável para a pergunta biológica.
   - Não, ou o objetivo inclui descoberta de genes/isoformas específicos da espécie (ex.: genes de resistência a inseticida, toxinas, proteínas de secreção) → monte de novo.

3. **O objetivo da pesquisa exige um catálogo completo e específico de transcritos da espécie (anotação funcional, busca de famílias gênicas, iscas para desenho de primers/CRISPR, caracterização estrutural de proteínas)?**
   - Sim → montagem de novo é o caminho correto, independentemente de existir ou não referência próxima, pois só ela captura a sequência real (não uma projeção) da espécie de interesse.

Na prática, para *A. gemmatalis* e *M. spectabilis*, os critérios 1 e 2 falham: nenhuma das duas tem genoma de referência de qualidade disponível, e os genomas mais próximos publicados — de outros lepidópteros e hemípteros — divergem demais para uma abordagem cross-species confiável em nível de transcrito individual. O critério 3, por outro lado, é positivo: ambos os projetos têm como objetivo caracterização funcional/estrutural específica da espécie. A montagem de novo é, portanto, a escolha correta, não uma segunda opção por falta de recurso.

### Desenho experimental

**Número de réplicas biológicas.** Para a etapa de *montagem* em si (gerar um catálogo de transcritos), réplicas biológicas não são estritamente necessárias; o que importa é maximizar a diversidade de transcritos capturados, geralmente combinando múltiplos tecidos/estágios/condições em uma única montagem de referência ("pooled assembly"). Já para a etapa subsequente de *expressão diferencial* sobre essa montagem, réplicas biológicas (não técnicas) são indispensáveis: o mínimo defensável na literatura de RNA-Seq é n=3 por condição, com n=4-6 preferível quando o orçamento permite, para dar poder estatístico ao DESeq2/edgeR. É uma distinção que quem está começando costuma confundir: "quantas amostras preciso sequenciar" tem respostas diferentes se a pergunta é "montar um transcriptoma" ou "comparar duas condições".

**Profundidade de sequenciamento.** A referência clássica mais citada para o número de reads necessário em montagem de novo é Francis et al. (2013, BMC Genomics, DOI 10.1186/1471-2164-14-167), não MacManes 2014 como às vezes se atribui de memória; o artigo de MacManes de 2014 (Frontiers in Genetics) trata de outro aspecto correlato, otimização de trimming, e não fixa esse número de reads. Francis et al. compararam remontagens de várias espécies animais não-modelo em profundidades crescentes de subamostragem e concluíram que montagens representativas já são obtidas com cerca de 20 milhões de reads para uma única biblioteca de tecido e cerca de 30 milhões para bibliotecas de organismo inteiro (misturando muitos tipos celulares), com ganho marginal de novos genes além de ~60 milhões de reads e risco crescente de acúmulo de erros de sequenciamento em transcritos muito expressados nessa faixa. Esse é o número (20-40M reads/amostra) citado no material do curso como "prática usual".

**A literatura mais recente concorda, com ressalvas.** Patterson et al. (2019, BMC Genomics, DOI 10.1186/s12864-019-5965-x) revisitaram a questão testando profundidades e tecnologias de sequenciamento mais amplas, concluindo que a quantidade de sequência exônica reconstruída satura em torno de 2-8 Gbp de dados (o que, a depender do tamanho de read, corresponde a uma faixa próxima da de Francis et al.), e que sequenciamento adicional além disso recupera majoritariamente transcritos de único éxon não anotados, de valor biológico questionável. Ou seja: reforça, e não substitui, a recomendação de Francis et al. de 2013. O guia mais recente de referência geral sobre o tema, Raghavan et al. (2022, Briefings in Bioinformatics, DOI 10.1093/bib/bbab563), não fixa um número específico de reads; os autores observam que estudos atuais rotineiramente sequenciam "centenas de milhões de reads" no total do projeto (não necessariamente por amostra) e remetem a Conesa et al. para recomendações de desenho experimental mais amplas.

Não encontrei, na busca realizada para este módulo, nenhuma publicação de 2024-2026 que revise ou substitua numericamente a faixa de 20-40M reads/amostra como consenso para de novo com short-read. A discussão que de fato mudou nesse período é sobre tecnologia (short-read vs. long-read), não sobre profundidade de short-read em si. Vale dizer isso de forma explícita aos alunos: a faixa de 20-40M é uma boa regra prática consolidada há mais de uma década, mas o estado da arte hoje está deslocando parte do problema (resolução de isoforma) para long-read complementar, não necessariamente pedindo mais reads Illumina.

**Single-end vs. paired-end.** Para montagem de novo, paired-end é fortemente preferível a single-end: a informação de distância de inserto entre as duas pontas ajuda o assemblador a resolver repetições, distinguir isoformas com éxons compartilhados e reduzir fragmentação. Single-end só se justifica em cenários de orçamento muito restrito, aceitando qualidade de montagem inferior.

**Tamanho de leitura.** 2x100 a 2x150 pb (paired-end) é o padrão atual em plataformas Illumina (NovaSeq/NextSeq) e equilibra custo, cobertura de éxons curtos e capacidade do assemblador de resolver junções.

**Stranded vs. unstranded.** Bibliotecas *stranded* (com preservação de informação de fita, via kits dUTP ou similares) permitem determinar a orientação transcricional de cada read. O impacto downstream é direto: sem essa informação, o assemblador e as ferramentas de quantificação não conseguem distinguir transcritos antisenso sobrepostos ou genes vizinhos em fitas opostas, o que infla artificialmente a contagem de "genes" fundidos incorretamente e reduz a precisão da anotação de UTRs. Para montagem de novo de organismos não-modelo, onde não há anotação prévia para desambiguar depois, biblioteca stranded é a recomendação padrão do curso: o custo adicional é pequeno frente ao ganho de precisão.

### Illumina (short-read) vs. long-read (PacBio Iso-Seq / ONT)

Reads curtos (Illumina, ~100-150 pb) têm taxa de erro baixa (<0,1%) e custo por base muito baixo, mas raramente cobrem um transcrito inteiro: o assemblador precisa *inferir* a estrutura completa de éxons e isoformas a partir de fragmentos sobrepostos, o que é exatamente a origem da maior parte das dificuldades de heterozigosidade/isoforma/redundância descritas acima. Reads longos (PacBio Iso-Seq, ou ONT com protocolos de cDNA/RNA direta) sequenciam a molécula de cDNA de ponta a ponta em uma única leitura (full-length transcript), eliminando a necessidade de inferência de conectividade entre éxons distantes: cada read já É, por definição, uma isoforma candidata completa.

O que muda no fluxo de trabalho com long-read:

- A etapa de "montagem" propriamente dita (grafo de de Bruijn) é substituída por *clustering* de reads full-length seguido de geração de consenso (pipelines como IsoSeq3/Cupcake para PacBio, ou RNA-Bloom2/isONform/RATTLE para ONT e cenários de novo sem referência).
- O problema de resolução de isoforma, que é o ponto mais frágil de assemblers short-read como Trinity, é resolvido de forma muito mais direta.
- Em compensação, a taxa de erro por read é maior (embora tecnologias mais recentes, como PacBio HiFi e ONT com basecalling atualizado, tenham reduzido bastante essa diferença), e o custo por Gb ainda costuma ser mais alto que Illumina, o que limita a profundidade típica de projetos long-read puros.
- Uma avaliação sistemática recente (Genome Biology, 2026, DOI 10.1186/s13059-026-04001-5) comparou os montadores de novo long-read RATTLE, RNA-Bloom2 e isONform contra Trinity (short-read) em dados simulados, spike-ins de referência conhecida e dados reais de humano e ervilha (sem uso de genoma de referência). RNA-Bloom2 combinado com Corset para clustering de transcritos foi a combinação com melhor desempenho em acurácia e eficiência computacional entre as opções long-read testadas, um indicativo forte de que o ecossistema de ferramentas de novo para long-read amadureceu o suficiente para uso em projetos reais de organismos sem referência, embora os próprios autores destaquem a falta histórica de benchmarks e protocolos estabelecidos para esse cenário específico (reference-free).
- Em paralelo, o Nature Methods (2024, DOI 10.1038/s41592-024-02298-3) publicou um benchmark sistemático de métodos de RNA-Seq long-read para identificação e quantificação de transcritos, no contexto do consórcio LRGASP (Long-read RNA-Seq Genome Annotation Assessment Project), um esforço comunitário relevante para quem for adotar long-read no futuro, ainda que focado majoritariamente em cenários com genoma de referência disponível.
- Uma abordagem híbrida (short-read para profundidade/quantificação + long-read para resolução de isoforma, com sequenciamento da mesma amostra ou pool nas duas tecnologias) é a direção que a literatura mais recente aponta como estado da arte para organismos não-modelo com orçamento que comporte as duas plataformas, mas não é obrigatória. Um projeto sólido em short-read puro (como as montagens Trinity de *A. gemmatalis* e *M. spectabilis* usadas neste curso) continua sendo cientificamente válido e é o que a maioria dos laboratórios com orçamento típico de pós-graduação no Brasil consegue executar.

### Fontes de dados públicos

Para localizar dados de RNA-Seq de insetos não-modelo já publicados (útil para desenho experimental comparativo, ou para complementar dados próprios com dados públicos da mesma espécie ou de espécies próximas):

- **NCBI SRA (Sequence Read Archive):** busca por nome científico do organismo + "RNA-Seq" em https://www.ncbi.nlm.nih.gov/sra; o BioProject/BioSample associado geralmente documenta desenho experimental, número de réplicas e protocolo de biblioteca usados por outros grupos, útil como referência de desenho antes de submeter seu próprio projeto.
- **ENA (European Nucleotide Archive):** espelha a maior parte do conteúdo do SRA e frequentemente é mais rápido para download de FASTQ diretamente (sem precisar do `sra-tools`/`fasterq-dump`), via https://www.ebi.ac.uk/ena.
- Para insetos de importância agrícola especificamente, vale também checar bancos genômicos especializados (i5k Workspace, InsectBase) que agregam metadados de projetos de sequenciamento de insetos não-modelo, mesmo quando o dado bruto está hospedado no SRA/ENA.
- **Nota de infraestrutura para o público deste curso:** como mencionado em outros módulos práticos, o acesso a FTP/HTTPS de bancos externos pode estar bloqueado em alguns ambientes de servidor institucional (exceto GitHub). Nesse caso, a estratégia é baixar os dados públicos em uma máquina sem essa restrição e transferir por outro meio, o que não é um problema do desenho experimental em si, mas afeta a logística de obtenção de dados de comparação.

## Prática usual vs. estado da arte (2024-2026)

| Aspecto | Prática usual/consolidada | Estado da arte 2024-2026 | Fonte (DOI) |
|---|---|---|---|
| Profundidade de sequenciamento (short-read) | 20-40 milhões de reads paired-end por amostra; platô de descoberta de genes além de ~60M | Nenhuma revisão numérica publicada em 2024-2026 encontrada na busca realizada; recomendação permanece a mesma, discussão deslocou-se para tecnologia (long-read complementar), não para volume de short-read | Francis et al. 2013, BMC Genomics, 10.1186/1471-2164-14-167 |
| Confirmação independente da faixa de profundidade | — | Estudo com profundidades/tecnologias mais amplas confirma platô de sequência exônica em 2-8 Gbp; ganho adicional é majoritariamente transcritos de éxon único não anotados | Patterson et al. 2019, BMC Genomics, 10.1186/s12864-019-5965-x |
| Guia de referência geral para desenho/montagem de novo | Recomendações dispersas em múltiplos artigos-fonte | Guia consolidado 2022 (dentro da janela de literatura "recente"); não fixa número de reads, remete a Conesa et al. para desenho experimental | Raghavan et al. 2022, Briefings in Bioinformatics, 10.1093/bib/bbab563 |
| Resolução de isoformas | Inferência indireta via grafo de de Bruijn a partir de reads curtos (Trinity), sensível a erro em regiões de splicing complexo | Montadores de novo dedicados a long-read (RNA-Bloom2+Corset, isONform, RATTLE) alcançam melhor acurácia/eficiência que Trinity em benchmark controlado sem uso de genoma de referência | Evaluation, Genome Biology 2026, 10.1186/s13059-026-04001-5 |
| Benchmarking comunitário de long-read | Ferramentas avaliadas isoladamente por cada grupo desenvolvedor | Esforço comunitário sistemático (LRGASP) de benchmarking de identificação/quantificação de transcritos via long-read | Nature Methods 2024, 10.1038/s41592-024-02298-3 |
| Estratégia de sequenciamento recomendada para novos projetos não-modelo | Illumina short-read puro (2x100-150bp, stranded, paired-end) | Híbrido short-read (quantificação/profundidade) + long-read (resolução de isoforma) quando orçamento permite; short-read puro continua válido e é o cenário mais viável para a maioria dos laboratórios | Síntese das fontes acima; nenhuma fonte isolada recomenda abandonar short-read puro |

## Recomendação para o curso

Para o público deste minicurso (pós-graduandos/pesquisadores com orçamento típico de projeto de pós-graduação/pesquisa no Brasil, sem acesso rotineiro a long-read), a recomendação é:

1. **Manter Illumina paired-end, 2x100-150pb, biblioteca stranded, 20-40M reads/amostra** como padrão de desenho experimental. Essa recomendação está bem estabelecida há mais de uma década (Francis et al. 2013) e foi reconfirmada por trabalho independente em 2019 (Patterson et al.), sem nenhuma revisão numérica publicada em 2024-2026 que a contradiga. É também o desenho compatível com os dois casos reais usados no curso.
2. **Combinar amostras de múltiplos tecidos/condições em uma única montagem de referência (pooled assembly)** para maximizar a diversidade de transcritos capturados, separando réplicas biológicas apenas na etapa posterior de quantificação/expressão diferencial, não na montagem.
3. **Apresentar long-read (Iso-Seq/ONT) como complemento de estado da arte, não como pré-requisito.** É importante que os alunos saibam que essa opção existe e para que problema específico ela resolve (resolução de isoforma), mas adotar isso como exigência do curso excluiria a maioria dos cenários reais de bancada de pós-graduação. O material deve deixar claro que um projeto Trinity short-read bem executado (como os dois organismos-caso) é cientificamente sólido e publicável, não uma solução de segunda classe.
4. **Ser transparente sobre a lacuna de evidência.** Ao apresentar a faixa de 20-40M reads, o instrutor deve dizer explicitamente aos alunos que essa é uma prática consolidada e não revisada por evidência publicada recente. Isso evita que o número seja repassado como se fosse dogma imutável, e modela a postura correta de checar a literatura antes de fixar um parâmetro de desenho experimental.

## Aplicação aos organismos-caso

Nenhum dos dois organismos-caso deste curso tem um SRA accession registrado na memória do projeto: as montagens Trinity existentes (42.372 transcritos para *A. gemmatalis*; 90.344 genes/103.560 transcritos para *M. spectabilis*) partiram de FASTQ já processados/fasta já montado, sem que os detalhes de desenho experimental (número exato de reads brutas, réplicas usadas na montagem) estejam documentados nos materiais do curso. Por isso, a etapa de desenho experimental é tratada aqui de forma didática e geral, não como uma reconstrução retroativa de como esses dois projetos específicos foram desenhados.

Dito isso, os números reais de saída servem como referência de ordem de grandeza esperada para organismos de tamanho genômico/complexidade transcricional semelhante:

- **Mahanarva spectabilis** (Hemiptera, glândula salivar): 90.344 genes / 103.560 transcritos brutos do Trinity, N50 = 738 pb, BUSCO C=96,2% em insecta_odb10, reduzindo para 12.445 proteínas após TransDecoder. A proporção transcritos:genes (~1,15) é relativamente baixa, compatível com uma montagem de um único tecido (glândula salivar) e não de um pool multi-tecido, o que reforça o ponto de que a escolha de quais amostras compor na montagem (critério de desenho discutido acima) tem efeito direto e mensurável no perfil final de isoformas reconstruídas. O BUSCO de 96,2% é um exemplo concreto de como a completude é avaliada na ausência de genoma de referência (usando o conjunto ortólogo esperado para Insecta), ilustrando diretamente o desafio nº 4 descrito acima ("ausência de referência para QC").
- **Anticarsia gemmatalis** (Lepidoptera): 42.372 transcritos, número bruto Trinity sem redução de redundância documentada no material disponível. Serve em módulos posteriores (Módulo 3, redução de redundância e avaliação) como caso para discutir quanto desse total é redundância alélica/isoforma versus genes distintos reais.

Para ambos os casos, ao ensinar desenho experimental de forma prospectiva ("se você fosse desenhar um novo experimento de RNA-Seq para uma espécie de inseto não-modelo hoje"), a recomendação da seção anterior (Illumina paired-end stranded, 20-40M reads/amostra, múltiplos tecidos/condições combinados na montagem) é o que o curso orienta os alunos a aplicar, inclusive como exercício prático de desenho para uma hipotética terceira réplica ou tecido desses mesmos organismos.

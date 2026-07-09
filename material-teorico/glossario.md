# Glossário

> ~100 termos usados ao longo dos 8 módulos do minicurso, em ordem alfabética. As definições foram escritas pensando no público do curso — pós-graduandos e pesquisadores que ainda não têm domínio prévio das ferramentas específicas.

**Adaptador (adapter)** — sequência sintética ligada às pontas de um fragmento de DNA/cDNA durante a preparação da biblioteca. O sequenciador precisa dela para ler o fragmento, mas ela deve sair no trimming.

**Anotação funcional** — atribuir a uma sequência montada uma função biológica provável (nome de gene, domínio, via metabólica, termo GO), geralmente por similaridade com bancos de referência.

**BLAST** — sigla de Basic Local Alignment Search Tool. Algoritmo clássico de busca por similaridade de sequência, hoje majoritariamente substituído por DIAMOND em escala de transcriptoma pela diferença de velocidade.

**BUSCO** — Benchmarking Universal Single-Copy Orthologs. Avalia a completude de uma montagem verificando se um conjunto de genes ortólogos de cópia única, esperado para o clado, está presente e íntegro.

**BUSCO C/S/D/F/M** — notação de saída do BUSCO: Complete, Single-copy, Duplicated, Fragmented, Missing.

**Backbone (do curso)** — pipeline conceitual central (QC → montagem → avaliação → anotação → quantificação → interpretação) que estrutura os módulos.

**Back-mapping** — remapeamento das reads originais contra a própria montagem, para medir que fração delas está representada; um proxy útil de completude.

**BioProject / BioSample** — identificadores do NCBI: o primeiro agrupa um projeto de sequenciamento inteiro, o segundo identifica uma amostra biológica específica dentro dele.

**BUSCO odb10 / odb12** — versões das linhagens de referência do OrthoDB usadas pelo BUSCO. odb12 é a geração mais recente (2024-2025); odb10 continua funcional, mas está congelada.

**Butterfly** — módulo final do Trinity, resolve isoformas/variantes dentro de cada componente do grafo e produz os transcritos finais.

**Chrysalis** — módulo intermediário do Trinity que agrupa contigs do Inchworm em componentes e reconstrói o grafo de De Bruijn completo por componente.

**cDNA** — DNA complementar, sintetizado a partir de um molde de RNA pela transcriptase reversa. É isso, e não o RNA em si, que o RNA-Seq de fato sequencia — exceto nos protocolos de RNA direta da ONT.

**CD-HIT-EST** — agrupa (clustering) sequências de nucleotídeos por identidade; usada para reduzir a redundância de transcritos.

**Componente (De Bruijn)** — subgrafo conectado de k-mers gerado pelo Chrysalis do Trinity, correspondente a um cluster de contigs relacionados.

**Contaminação (transcriptômica)** — presença de sequências não pertencentes ao organismo-alvo (bactérias, fungos, DNA humano, endossimbiontes) misturadas ao transcriptoma montado.

**Contig** — sequência contígua montada a partir da sobreposição de reads ou k-mers.

**Corset** — agrupa transcritos em "genes" com base em reads compartilhadas e correlação de expressão; não depende do montador de origem.

**Cobertura (coverage)** — número médio de vezes que uma dada posição da sequência é representada pelas reads sequenciadas.

**De Bruijn, grafo de** — estrutura de dados em que k-mers são nós conectados por sobreposição de k-1 bases; base algorítmica de Trinity e rnaSPAdes.

**De novo (montagem)** — reconstrução de transcritos/sequências diretamente a partir de reads, sem uso de genoma de referência.

**DESeq2** — pacote R/Bioconductor para análise de expressão diferencial, baseado em modelo binomial negativo com shrinkage bayesiano.

**DIAMOND** — busca por similaridade de sequência proteica, ordens de magnitude mais rápida que o BLAST; hoje é o padrão em anotação de transcriptomas.

**Endossimbionte** — microrganismo (geralmente bactéria) que vive dentro das células de um hospedeiro em relação simbiótica, comum em insetos sugadores de seiva (ex.: Buchnera, Sulcia, Wolbachia).

**EnTAP** — pipeline de anotação funcional dedicado a eucariotos não-modelo, com filtro de contaminantes/expressão/frame integrado.

**ENA** — European Nucleotide Archive; repositório europeu equivalente ao NCBI SRA para dados de sequenciamento brutos.

**Enriquecimento funcional (ORA)** — teste estatístico (em geral hipergeométrico) que checa se um conjunto de genes de interesse (ex.: DEGs) concentra mais do que o esperado de um termo GO ou via KEGG, em comparação a um universo de referência.

**EvidentialGene (tr2aacds)** — pipeline que classifica transcritos em primário/alternativo/redundante com base em estrutura de ORF/CDS, usado para consolidar montagens de múltiplos montadores/k-mers.

**ExN50** — métrica do Trinity que calcula N50 ponderado pela contribuição de expressão de cada transcrito, mais informativa que N50 bruto para transcriptomas.

**FASTQ** — formato de arquivo padrão para armazenar reads de sequenciamento junto com seus escores de qualidade (Phred).

**FastQC** — inspeciona a qualidade de reads brutas e gera um relatório por amostra (Phred, adaptador, GC, duplicação).

**fastp** — combina filtragem de qualidade, remoção de adaptador e geração de relatório num único passo de trimming/QC.

**FCS-GX (Foreign Contamination Screen)** — usada pelo NCBI para triagem de contaminação em montagens, antes do depósito no TSA/GenBank.

**FAIR (princípios)** — Findable, Accessible, Interoperable, Reusable; conjunto de princípios para tornar dados (e, mais recentemente, software) cientificamente reutilizáveis.

**Filogenia** — reconstrução das relações evolutivas entre sequências ou espécies, tipicamente representada como uma árvore.

**Frame (quadro de leitura)** — um dos três possíveis registros de leitura de códons de três nucleotídeos numa sequência, relevante para tradução correta de uma ORF.

**Genome-guided (montagem)** — reconstrução de transcritos a partir do alinhamento de reads contra um genoma de referência.

**Gene Ontology (GO)** — vocabulário controlado hierárquico que descreve função molecular, processo biológico e componente celular de genes/proteínas.

**GSEA (Gene Set Enrichment Analysis)** — método de enriquecimento funcional que usa o ranking completo de genes (não apenas um corte de significância) para detectar tendências sutis em conjuntos de genes.

**Grafo de De Bruijn** — ver "De Bruijn, grafo de".

**Heterozigosidade** — presença de variantes alélicas distintas num locus, que pode confundir montadores de novo ao gerar bifurcações espúrias no grafo.

**HMMER / hmmsearch** — busca domínios de proteína usando modelos ocultos de Markov (HMM); usada contra o banco Pfam.

**IMRaD** — estrutura padrão de artigo científico: Introduction, Methods, Results, and Discussion.

**Inchworm** — primeiro módulo do Trinity, monta contigs lineares gulosos a partir de k-mers.

**InterProScan** — integra múltiplos bancos de assinaturas de proteína (Pfam, PROSITE, PRINTS etc.) num único relatório de domínios/GO/pathways.

**IQ-TREE** — software de inferência filogenética por máxima verossimilhança com seleção automática de modelo (ModelFinder) e suporte de ramo via ultrafast bootstrap.

**Isoforma** — variante de transcrito de um mesmo gene, geralmente produzida por splicing alternativo.

**Iso-Seq** — protocolo/pipeline da PacBio para sequenciamento e processamento de transcritos completos (full-length) via long-read.

**k-mer** — subsequência de comprimento fixo *k* extraída de uma read, unidade básica dos grafos de De Bruijn.

**KEGG (Kyoto Encyclopedia of Genes and Genomes)** — banco de dados de vias metabólicas e de sinalização, usado para anotação funcional e enriquecimento.

**KO (KEGG Orthology)** — identificador de ortologia funcional dentro do sistema KEGG, atribuído por ferramentas como eggNOG-mapper.

**Log2FoldChange (LFC)** — medida de magnitude de mudança de expressão entre duas condições, em escala logarítmica de base 2.

**LRGASP** — Long-read RNA-Seq Genome Annotation Assessment Project; consórcio internacional de benchmarking de métodos de RNA-Seq long-read.

**MAFFT** — alinha múltiplas sequências; etapa intermediária comum antes da inferência filogenética.

**MAplot** — gráfico que relaciona expressão média (eixo x) e log2FoldChange (eixo y) de genes num experimento de expressão diferencial.

**MultiQC** — agrega relatórios de múltiplas amostras/ferramentas (FastQC, fastp, Salmon etc.) num único relatório HTML navegável.

**N50** — métrica de contiguidade de montagem: comprimento tal que 50% do total de bases montadas estão em contigs desse tamanho ou maior; enganosa se usada isoladamente.

**NCBI SRA (Sequence Read Archive)** — repositório do NCBI para dados brutos de sequenciamento de alto rendimento.

**Não-modelo (organismo)** — espécie sem genoma de referência publicado/bem anotado.

**Normalização in silico** — redução computacional da profundidade de cobertura de regiões muito expressas antes da montagem, para economizar recursos sem perder complexidade de k-mer.

**Ortogrupo** — conjunto de genes de duas ou mais espécies descendentes de um único gene no ancestral comum dessas espécies.

**Ortólogo** — gene relacionado a outro por um evento de especiação (divergência entre espécies).

**Ortologia** — a relação entre esses genes: vínculo evolutivo entre genes de espécies diferentes, derivados de um ancestral comum por especiação.

**ORF (Open Reading Frame)** — trecho de sequência com potencial para codificar uma proteína, delimitado por um possível início e fim de tradução.

**OrthoFinder** — infere ortogrupos e a árvore de espécies a partir dos proteomas de múltiplas espécies.

**OrgDb** — associa, para uma dada espécie, identificadores de genes a anotações (GO, KEGG etc.); pacote Bioconductor tipicamente inexistente para organismos não-modelo.

**Parálogo** — gene relacionado a outro por um evento de duplicação gênica dentro da mesma linhagem.

**Paired-end** — modo de sequenciamento em que ambas as pontas de um fragmento são lidas, gerando reads R1/R2 com informação de distância de inserto.

**Pathview** — sobrepõe dados de expressão diferencial em diagramas de vias do KEGG (pacote R).

**Pfam** — banco de dados de famílias de domínios de proteína, baseado em modelos ocultos de Markov.

**Phred score (Q score)** — medida logarítmica de confiança na chamada de base de sequenciamento; Q30 corresponde a 1 erro em 1000 bases.

**Platô (Ex50)** — no gráfico de ExN50, ponto em que o N50 ponderado por expressão se estabiliza, indicando que os transcritos mais expressos já estão bem representados.

**Pooled assembly** — montagem que combina reads de múltiplos tecidos/condições/estágios numa única referência transcriptômica.

**Protein language model (PLM)** — modelo de aprendizado profundo (ex.: ESM-2) treinado sobre sequências de proteína, usado para prever função sem depender de homólogo próximo em banco.

**Quasi-mapping** — estratégia de pseudo-alinhamento (usada pelo Salmon) que localiza posições compatíveis de uma read no transcriptoma sem produzir um alinhamento base-a-base completo.

**Rcorrector** — corrige erros de sequenciamento com base num grafo de k-mers confiáveis; recomendado antes do trimming em montagens de novo.

**Read-through** — quando o fragmento sequenciado é mais curto que o comprimento da read, fazendo o sequenciador "ler através" do inserto até o adaptador da outra ponta.

**Redundância (de transcritos)** — presença de múltiplas sequências quase idênticas na montagem representando o mesmo transcrito/gene, por isoforma real ou artefato de montagem.

**Reference-free / reference-based** — respectivamente, análise que não depende de um genoma de referência vs. análise que usa um genoma de referência como base.

**rnaQUAST** — gera estatísticas descritivas de qualidade para montagens de transcriptoma, com ou sem genoma de referência.

**rnaSPAdes** — montador de transcriptoma de novo baseado em múltiplos valores de k-mer, derivado do motor SPAdes.

**RNA-Bloom2** — montador de novo reference-free para dados long-read (PacBio/ONT).

**RT-PCR (qRT-PCR)** — técnica de PCR quantitativa em tempo real usada para validar experimentalmente níveis de expressão gênica detectados por RNA-Seq.

**Salmon** — quantifica expressão por quasi-mapping, sem precisar de alinhamento completo contra a referência.

**Shrinkage (LFC)** — técnica estatística (ex.: `apeglm` no DESeq2) que ajusta estimativas de log2FoldChange para genes de baixa contagem, reduzindo ruído nas estimativas.

**Single-copy ortholog** — gene ortólogo presente em cópia única, base das métricas de completude do BUSCO.

**SortMeRNA** — filtra reads de rRNA contra bancos de referência (SILVA, Rfam).

**Splicing alternativo** — processo pelo qual um mesmo gene gera múltiplas isoformas de mRNA por combinação diferencial de éxons.

**Stranded (biblioteca)** — protocolo de preparo de biblioteca que preserva a informação de qual fita de DNA/RNA originou cada read.

**TERM2GENE / TERM2NAME** — tabelas customizadas usadas pelo `enricher()`/`GSEA()` do clusterProfiler para enriquecimento funcional sem depender de um OrgDb.

**TransDecoder** — identifica ORFs e regiões codificantes candidatas dentro de transcritos montados de novo.

**Transcriptoma** — conjunto completo de moléculas de RNA transcritas num organismo/tecido/condição específicos, em um dado momento.

**TransPi** — pipeline Nextflow que executa múltiplos montadores e consolida com EvidentialGene, reportando métricas padronizadas de qualidade.

**TransRate** — avalia montagens por score de contig, sem depender de referência.

**Trimming** — remoção/corte de bases de baixa qualidade, adaptadores ou artefatos das extremidades de uma read.

**Trinity** — montador de transcriptoma de novo baseado em grafo de De Bruijn, composto pelos módulos Inchworm, Chrysalis e Butterfly.

**Trinotate** — suíte clássica de anotação funcional integrada ao ecossistema Trinity.

**TSA (Transcriptome Shotgun Assembly)** — categoria de depósito no NCBI/GenBank específica para montagens de transcriptoma.

**tximport** — importa quantificações de ferramentas como o Salmon e agrega por gene via um mapa transcrito→gene (pacote R).

**Universo (enrichment universe)** — conjunto de genes/background contra o qual um teste de enriquecimento funcional compara o conjunto de interesse; definir incorretamente é o erro mais comum de enriquecimento em não-modelo.

**Unstranded (biblioteca)** — protocolo de biblioteca que não preserva informação de fita de origem.

**UTR (Untranslated Region)** — região do transcrito (5' ou 3') que não é traduzida em proteína, mas frequentemente contém elementos regulatórios.

**Volcano plot** — gráfico que relaciona magnitude (log2FoldChange) e significância estatística (−log10 padj) de genes num teste de expressão diferencial.

**Zenodo** — repositório de arquivamento científico que integra com GitHub para gerar DOIs permanentes de releases de código/dados.

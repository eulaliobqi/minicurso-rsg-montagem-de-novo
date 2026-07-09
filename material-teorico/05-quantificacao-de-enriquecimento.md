# Módulo 5 — Quantificação, Expressão Diferencial e Enriquecimento

**Duração prevista:** 55 min · **Pré-requisito:** Módulo 4 (montagem + anotação funcional já concluídas para *A. gemmatalis* e *M. spectabilis*)

## Conceitos-chave

Depois de montar um transcriptoma de novo (Módulo 3) e anotá-lo funcionalmente (Módulo 4), o objetivo deste módulo é responder à pergunta biológica: **quais transcritos mudam de expressão entre condições, e o que essas mudanças significam funcionalmente?**

Em organismo não-modelo essa pergunta se divide em três problemas técnicos encadeados, cada um com uma armadilha específica:

1. **Quantificar abundância de transcritos** sem um genoma de referência para alinhar contra, resolvido por pseudo-alinhamento/quasi-mapping (Salmon), que evita o custo computacional e os vieses de um alinhador completo (STAR/HISAT2) quando não há coordenadas genômicas confiáveis.
2. **Agregar transcritos em unidades biológicas ("genes")**: um problema que não existe em organismo modelo (onde o GTF já define genes), mas é central em montagem de novo. Trinity produz múltiplas isoformas por "gene" (`_i1`, `_i2`...) e às vezes contigs redundantes que só podem ser agrupados a posteriori, por padrão de leitura compartilhada (Corset) ou por hierarquia de isoforma+gene do próprio Trinity (via `tximport`).
3. **Testar diferença de expressão e interpretar biologicamente o resultado**: DESeq2/edgeR fazem a parte estatística de forma robusta mesmo sem anotação completa, mas a interpretação funcional (GO/KEGG) exige uma base de anotação de organismo (OrgDb) que **não existe** no Bioconductor para *A. gemmatalis* nem para *M. spectabilis*. Esse é o desafio central do módulo, detalhado adiante.

O fio condutor pedagógico é que cada etapa do fluxo padrão (Salmon → tximport/Corset → DESeq2 → clusterProfiler) tem uma versão "de livro-texto" pensada para organismo modelo e uma adaptação necessária para não-modelo. Entender onde e por que a adaptação é necessária evita o erro mais comum da área: tratar um pipeline de organismo modelo como plug-and-play para de novo.

## Ferramentas e fluxo

### Salmon — quantificação sem alinhamento

Salmon (Patro et al. 2017, *Nature Methods* 14:417–419, DOI [10.1038/nmeth.4197](https://doi.org/10.1038/nmeth.4197)) quantifica abundância de transcritos por **quasi-mapping**: em vez de produzir um alinhamento base-a-base, encontra o conjunto mínimo de posições no transcriptoma compatíveis com cada read (via k-mer/suffix array) e resolve a atribuição de reads multi-mapeadas por inferência estatística (EM/VBEM) num modelo que já embute vieses de sequenciamento.

Fluxo de comando:

```bash
mamba activate trinity   # ambiente com Salmon instalado

# 1. Construir índice a partir do transcriptoma de novo (Trinity.fasta)
salmon index \
  -t Trinity.fasta \
  -i salmon_index \
  -k 31 \
  -p 16

# 2. Quantificar cada amostra (paired-end)
salmon quant \
  -i salmon_index \
  -l A \
  -1 amostra_R1.fastq.gz -2 amostra_R2.fastq.gz \
  --validateMappings \
  --gcBias --seqBias \
  -p 16 \
  -o quant/amostra
```

Parâmetros que fazem diferença prática:
- `--validateMappings`: ativa um alinhamento local (SIMD) sobre os hits de quasi-mapping para eliminar falsos positivos espúrios, recomendado por padrão desde Salmon ≥0.14, especialmente importante em transcriptoma de novo, que tem mais redundância/isoformas parecidas do que um transcriptoma de referência curado.
- `--gcBias` / `--seqBias`: corrigem viés de conteúdo GC e viés de sequência específica de posição (hexâmeros de priming), respectivamente. Sem eles, transcritos com GC atípico (comum em regiões repetitivas ou UTRs longas de insetos) têm abundância sistematicamente sub ou superestimada. Salmon foi o primeiro quantificador a corrigir viés de GC no nível de transcrito, o que reduz falsos positivos na etapa de DE.
- `-l A`: detecção automática do tipo de biblioteca (orientação/strandedness); evitar fixar manualmente a menos que já se saiba o protocolo do kit.

### Corset e tximport — de transcritos a "genes"

Sem GTF, dois caminhos equivalentes resolvem o agrupamento:

**Corset** (Davidson & Oshlack 2014, *Genome Biology* 15:410, DOI [10.1186/s13059-014-0410-6](https://doi.org/10.1186/s13059-014-0410-6)) agrupa contigs hierarquicamente com base em reads compartilhadas entre eles e na correlação de expressão. É agnóstico ao montador, então funciona tanto para dúvidas sobre isoformas quanto para redundância de contigs que o próprio Trinity não colapsou.

**tximport** (Soneson, Love & Robinson 2015, *F1000Research* 4:1521, DOI [10.12688/f1000research.7563.2](https://doi.org/10.12688/f1000research.7563.2)) importa as quantificações do Salmon já agregando por gene, usando o `gene_trans_map` que o próprio Trinity gera (`Trinity.fasta.gene_trans_map`), com `countsFromAbundance = "lengthScaledTPM"`. Essa opção corrige o efeito de comprimento diferencial entre isoformas de um mesmo gene antes de entregar a matriz de contagem ao DESeq2, prática recomendada por Soneson et al. porque comparar contagens brutas de "genes" com número variável de isoformas por amostra infla falsos positivos.

Na prática do laboratório, para Trinity + Salmon, `tximport` é o caminho mais direto (o mapeamento gene→transcrito já existe); Corset é preferido quando não há esse mapeamento de montador (ex.: montagens redundantes de múltiplos assemblers) ou quando se quer detectar redundância de contigs não capturada pelo Trinity.

```r
library(tximport)
txi <- tximport(files = quant_files,
                 type = "salmon",
                 tx2gene = gene_trans_map,
                 countsFromAbundance = "lengthScaledTPM")
```

### DESeq2, edgeR e limma-voom — expressão diferencial

**DESeq2** (Love, Huber & Anders 2014, *Genome Biology* 15:550, DOI [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)) modela contagens com distribuição binomial negativa, estima dispersão por shrinkage empírico-Bayesiano (compartilhando informação entre genes com contagens médias parecidas) e aplica shrinkage de log2FoldChange (`lfcShrink`, método `apeglm` recomendado) para estabilizar estimativas de genes de baixa contagem, crítico em não-modelo, onde muitos transcritos têm expressão baixa/variável.

Fluxo padrão do laboratório:

```r
library(DESeq2)
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condicao)
dds <- dds[rowSums(counts(dds)) >= 10, ]     # filtragem de baixa contagem
dds <- DESeq(dds)
res <- lfcShrink(dds, coef = "condicao_exposto_vs_controle", type = "apeglm")
# padj <= 0.05 e |log2FoldChange| >= 1 (padrão do laboratório)
```

**edgeR** usa abordagem semelhante (binomial negativa + normalização TMM) mas com shrinkage de dispersão via *quantile-adjusted conditional maximum likelihood*, e tende a ser um pouco mais sensível (menos conservador) que DESeq2 quando o número de réplicas é muito baixo (n=2–3 por grupo), porque seu empirical Bayes compartilha informação entre genes de forma mais agressiva. Nenhum dos dois é estritamente "melhor"; a escolha usual do laboratório é DESeq2 por integração nativa com `tximport`, `lfcShrink` e `vst()` para PCA.

**limma-voom** é preferível quando: (a) o desenho experimental é complexo (múltiplos fatores/blocos, medidas repetidas) e se quer aproveitar o arcabouço de modelo linear do `limma`; (b) o número de amostras é grande o suficiente (tipicamente >10 por grupo) para a transformação `voom` (que converte contagens em log-CPM com peso de precisão modelado por média-variância) estimar bem a tendência média-variância. Não é a escolha default do laboratório para os experimentos atuais (poucas réplicas, desenho simples exposto-vs-controle).

### clusterProfiler / topGO / goatools — enriquecimento funcional

**clusterProfiler 4.0** (Wu et al. 2021, *The Innovation* 2:100141, DOI [10.1016/j.xinn.2021.100141](https://doi.org/10.1016/j.xinn.2021.100141)) é a ferramenta padrão do laboratório porque sua função genérica `enricher()` (ORA) e `GSEA()` **não exigem um OrgDb**: aceitam diretamente tabelas `TERM2GENE`/`TERM2NAME` construídas pelo usuário, exatamente o que se precisa para organismo sem anotação no Bioconductor (ver seção seguinte). `topGO` e `goatools` (Python) são alternativas válidas — `goatools` é frequentemente citado em pipelines totalmente Python, útil se o restante do fluxo de anotação já é Python, como eggNOG-mapper/InterProScan —, mas o laboratório padroniza em `clusterProfiler` por já ter o `enricher()` testado em produção no projeto *S. frugiperda*.

## O desafio do "sem OrgDb" em não-modelo

Nenhum dos dois organismos deste minicurso, *A. gemmatalis* e *M. spectabilis*, tem pacote `OrgDb` no Bioconductor (o equivalente de `org.Mm.eg.db` ou `org.At.tair.db`). Sem OrgDb, funções como `enrichGO()`/`enrichKEGG()` no modo "de livro-texto" simplesmente não funcionam, porque elas esperam consultar esse banco para saber quais genes pertencem a quais termos GO/KEGG.

A solução consolidada, já implementada e testada em produção no projeto irmão RNA-Seq-not-model (*S. frugiperda*, `C:\Users\eulal\.claude\RNA-Seq-not-model`), é construir esse mapeamento **na mão**, a partir da própria anotação funcional gerada no Módulo 4 (eggNOG-mapper). Passo a passo conceitual:

1. **Ponto de partida:** output do eggNOG-mapper (`*.emapper.annotations`), que para cada proteína/transcrito anotado já lista colunas com termos GO (coluna `GOs`) e ortólogos KEGG (coluna `KEGG_ko`).
2. **Resolver a granularidade transcrito → gene:** como a anotação roda sobre proteínas preditas (TransDecoder, ex. `TRINITY_DN1234_c0_g1_i1.p1`), é preciso mapear de volta para o nível de "gene" definido no passo de agrupamento (Corset/tximport), usando o `gene_trans_map` do Trinity. Sem esse passo, a mesma unidade biológica aparece fragmentada em múltiplos IDs de transcrito/proteína e o enriquecimento perde poder estatístico.
3. **Construir `TERM2GENE`:** uma tabela de duas colunas (termo GO ou KO; gene_id), obtida explodindo a lista de termos separados por vírgula da coluna `GOs`/`KEGG_ko` do eggNOG (um gene pode ter múltiplos termos).
4. **Construir `TERM2NAME`:** para GO, os nomes descritivos vêm do pacote `GO.db` (genérico, não específico de espécie, apenas define o vocabulário controlado GO, então pode ser usado para qualquer organismo); para KEGG, a própria coluna `Description` do eggNOG-mapper já fornece o nome do ortólogo.
5. **Definir o universo (`universe`) como o conjunto de genes efetivamente testados na análise de DE**, não o genoma inteiro nem uma espécie de referência (ver próxima seção: é o erro mais crítico do módulo).
6. **Chamar `enricher()`/`GSEA()` do clusterProfiler passando essas tabelas diretamente**, em vez de `enrichGO()`/`enrichKEGG()` com parâmetro `OrgDb=`.

Esse fluxo está implementado nos scripts `02_gene2go_build.R` (constrói os mapeamentos a partir do eggNOG + `gene_trans_map`) e `03_enrichment.R` (roda `enricher()`/`fgsea()` com esses mapeamentos) do projeto RNA-Seq-not-model, e é o template a ser replicado para *A. gemmatalis* e *M. spectabilis* neste minicurso.

## Prática usual vs. estado da arte (2024-2026)

| Aspecto | Prática usual/consolidada | Estado da arte 2024-2026 | Fonte (DOI) |
|---|---|---|---|
| Teste de DE em não-modelo, poucas réplicas | DESeq2 com shrinkage `apeglm`; edgeR como alternativa mais sensível quando n≤3/grupo | Estudo de replicabilidade em larga escala (18.000 subamostragens) confirma que resultados de DE/enriquecimento com coortes pequenas têm baixa replicabilidade, mas alta precisão mediana em vários datasets com >5 réplicas/grupo. Reforça, não substitui, DESeq2/edgeR como consenso; não emergiu um método substituto de adoção ampla | Degen & Medo 2025, *PLOS Comp. Biol.* 21:e1011630, DOI [10.1371/journal.pcbi.1011630](https://doi.org/10.1371/journal.pcbi.1011630) |
| Agregação transcrito→gene sem referência | Corset (agrupamento por leitura compartilhada) ou `tximport` com `gene_trans_map` do Trinity | Sem mudança de paradigma relatada 2024-2026; ambos seguem sendo os caminhos padrão citados em pipelines de novo recentes | Davidson & Oshlack 2014, DOI [10.1186/s13059-014-0410-6](https://doi.org/10.1186/s13059-014-0410-6); Soneson et al. 2015, DOI [10.12688/f1000research.7563.2](https://doi.org/10.12688/f1000research.7563.2) |
| Enriquecimento GO/KEGG sem OrgDb | Construir `TERM2GENE`/`TERM2NAME` manualmente a partir do eggNOG-mapper e usar `enricher()`/`GSEA()` genéricos do clusterProfiler; é de fato a prática usual reportada em fóruns e pipelines da comunidade (Bioconductor support, Biostars) | Surgiram pacotes dedicados que empacotam esse mesmo padrão: **FunFEA** (fungos, gera modelo de background diretamente do output do eggNOG-mapper) e uso mais sistemático de `AnnotationForge::makeOrgPackageFromNCBI()` para gerar um OrgDb customizado quando há orçamento de tempo para isso; para insetos sem genoma anotado no NCBI (caso de ambos organismos deste curso), a rota `enricher()` manual continua sendo a única viável | Charest et al. 2025, *BMC Bioinformatics* 26:138, DOI [10.1186/s12859-025-06164-7](https://doi.org/10.1186/s12859-025-06164-7); Wu et al. 2021, DOI [10.1016/j.xinn.2021.100141](https://doi.org/10.1016/j.xinn.2021.100141) |
| Quantificação sem alinhamento | Salmon quasi-mapping com `--validateMappings --gcBias --seqBias` | Sem mudança de paradigma relatada; Salmon segue sendo o quantificador de facto para RNA-Seq de novo em 2024-2026, sem substituto de adoção comparável | Patro et al. 2017, DOI [10.1038/nmeth.4197](https://doi.org/10.1038/nmeth.4197) |

**Leitura para o curso:** o consenso 2014-2016 (Salmon/DESeq2/tximport/Corset/clusterProfiler) permanece o estado da prática em 2026. O que mudou não foi a escolha de ferramenta, e sim (a) evidência mais robusta de que replicabilidade com poucas réplicas é limitada, o que reforça a necessidade de reportar resultados como "hipóteses a validar", não conclusões definitivas, e (b) o surgimento de pacotes que empacotam o padrão "eggNOG → enriquecimento sem OrgDb" que este laboratório já usa manualmente, validando a abordagem adotada.

## Erro comum: universo/background errado

O erro mais frequente, e mais silencioso, porque não gera mensagem de erro, apenas resultado inflado, é definir o `universe` do teste de enriquecimento incorretamente. Em `enricher()`/`enrichGO()`/`enrichKEGG()`, o `universe` é o conjunto de genes contra o qual a proporção observada nos genes DE é comparada estatisticamente (teste hipergeométrico).

**Erro típico:** usar como universo o genoma/transcriptoma inteiro de uma espécie modelo próxima (ex.: todos os genes anotados de *Drosophila melanogaster* como proxy para *A. gemmatalis*, ou todo o conjunto de genes com OrgDb disponível), ou mesmo todos os transcritos brutos da montagem Trinity, incluindo os que nunca passaram pelo filtro de expressão mínima do DESeq2.

**Por que isso é errado:** o teste hipergeométrico calcula a probabilidade de observar X genes de um termo GO entre os genes DE, dado quantos genes daquele termo existem no universo. Se o universo é artificialmente grande (genoma inteiro) ou vem de outra espécie, a proporção de "sucesso" (genes DE anotados por termo) é comparada contra uma base que **não reflete o que foi de fato testável no experimento**. Isso reduz sistematicamente o denominador de "genes do termo não-DE mas testados" incorretamente, ou infla o número total de genes candidatos que nunca tiveram chance real de aparecer como DE (porque nem estavam expressos na biblioteca ou nem foram anotados). O resultado prático: termos GO aparecem como "significativamente enriquecidos" com p-valores artificialmente baixos, que não sobreviveriam a um teste com o universo correto.

**Regra prática do laboratório (implementada em `03_enrichment.R`):** `universe` = o vetor de `gene_id` que sobrou após a filtragem de baixa contagem do DESeq2 e que tem pelo menos uma anotação (GO ou KEGG) via eggNOG, ou seja, exatamente o conjunto de genes que *poderiam*, em princípio, ter sido chamados como DE e anotados nesta análise. Nunca o genoma de referência de outra espécie, nunca a montagem Trinity bruta sem filtro.

## Exemplo aplicado (hipotético, claramente identificado como didático)

**Aviso:** o exemplo abaixo é inteiramente hipotético, construído apenas para ilustrar o fluxo de comandos e a interpretação de saída. Não representa um resultado experimental real deste laboratório.

Cenário didático: *A. gemmatalis* exposta a um inseticida (grupo `exposto`, n=3) vs. controle não exposto (grupo `controle`, n=3), com objetivo de identificar processos biológicos enriquecidos entre os genes diferencialmente expressos.

```bash
mamba activate trinity

# 1. Quantificação (Salmon) — repetir para as 6 amostras
salmon quant -i salmon_index -l A \
  -1 exposto_R1_1.fastq.gz -2 exposto_R1_2.fastq.gz \
  --validateMappings --gcBias --seqBias -p 16 -o quant/exposto_R1
# ... repetir para exposto_R2, exposto_R3, controle_R1, controle_R2, controle_R3
```

```r
# 2. tximport + DESeq2 (R, ambiente com DESeq2 instalado)
library(tximport); library(DESeq2)

txi <- tximport(quant_files, type = "salmon",
                 tx2gene = gene_trans_map, countsFromAbundance = "lengthScaledTPM")

coldata <- data.frame(condicao = factor(c("exposto","exposto","exposto",
                                           "controle","controle","controle"),
                                         levels = c("controle","exposto")))

dds <- DESeqDataSetFromTximport(txi, coldata, design = ~ condicao)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res <- lfcShrink(dds, coef = "condicao_exposto_vs_controle", type = "apeglm")
# Hipotético: 187 genes up, 94 down (padj <= 0.05, |log2FC| >= 1)
```

```r
# 3. Enriquecimento GO/KEGG sem OrgDb (mesmo padrão do RNA-Seq-not-model)
# 03_gene2go_build.R já rodado no Módulo 4 com o eggNOG-mapper de A. gemmatalis
Rscript scripts/03_enrichment.R \
  --results deseq2_results_all.tsv \
  --gene2go ag_gene2go.rds \
  --analysis go \
  --padj_cutoff 0.05 --lfc_cutoff 1.0 \
  --outdir enrichment_go/
```

**Interpretação didática do resultado hipotético:** suponha que o `go_bp_enrichment.tsv` resultante mostre enriquecimento significativo de termos como "resposta a estímulo xenobiótico" (GO:0009410) e "atividade de monooxigenase" (relacionada a citocromo P450) entre os genes up-regulados no grupo exposto. Esse padrão seria biologicamente coerente com uma resposta de detoxificação a inseticida. Mas, tratando-se de um exemplo hipotético com n=3 por grupo, o próximo passo real de qualquer análise assim seria (a) verificar se os termos sobrevivem a validação por qRT-PCR de genes-chave, e (b) reportar o resultado como hipótese, dado o que a literatura de replicabilidade (Degen & Medo 2025) mostra sobre poder estatístico limitado em coortes pequenas.

## Recomendação para o curso

1. **Demonstrar ao vivo apenas o Salmon `quant`** de uma amostra pequena (subsample de reads). O índice e a quantificação completa das 6+ amostras reais devem rodar em `screen` no servidor antes ou depois da aula, não durante.
2. **Focar a demonstração de R no par `tximport` → `DESeq2` → `lfcShrink`**, reforçando a leitura da tabela de resultados (`baseMean`, `log2FoldChange`, `padj`) antes de qualquer enriquecimento.
3. **Dedicar o tempo mais denso da aula ao bloco "sem OrgDb"**: abrir ao vivo o `02_gene2go_build.R` e `03_enrichment.R` do projeto RNA-Seq-not-model como estudo de caso real e funcional, mostrando a tabela `TERM2GENE` gerada a partir do eggNOG-mapper. Este é o ponto que mais distingue organismo não-modelo de organismo modelo, e onde a maioria dos tutoriais genéricos de RNA-Seq falha em cobrir.
4. **Fazer os alunos identificarem o erro de universo** num exercício rápido: dar dois `enricher()` rodados com universos diferentes (um correto, um com genoma inteiro) sobre o mesmo conjunto de genes DE hipotéticos, e pedir que expliquem por que os p-valores diferem.
5. **Deixar explícito, com o exemplo aplicado, que resultado de enriquecimento em experimento com poucas réplicas é hipótese, não conclusão**, citando Degen & Medo (2025), para calibrar a expectativa dos alunos antes que apliquem isso aos próprios dados de *A. gemmatalis*/*M. spectabilis*.

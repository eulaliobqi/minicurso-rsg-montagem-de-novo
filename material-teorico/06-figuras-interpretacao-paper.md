# Módulo 6 — Figuras, Interpretação Biológica, Reprodutibilidade e Escrita do Paper

> Duração: 60 min | Módulo final e mais abrangente do minicurso. Objetivo: sair de "tenho uma tabela de DEGs e um BUSCO score" para "tenho um manuscrito submetível, com dados depositados e pipeline reproduzível por terceiros".

## Figuras publicáveis

Cada figura de um paper de transcriptômica tem um papel narrativo específico. Escolher a figura errada para a pergunta errada é um dos erros mais comuns em primeira submissão.

| Figura | O que comunica | Quando usar | Pacote/ferramenta |
|---|---|---|---|
| **PCA / MDS de amostras** | Estrutura global da variação — as réplicas biológicas agrupam por condição/tecido? Há outlier técnico? | Sempre, logo após `vst()`/`rlog()` no DESeq2, antes de qualquer teste de expressão diferencial. Vai em Results como "controle de qualidade" | `plotPCA()` (DESeq2) ou `ggplot2` customizado |
| **Heatmap de expressão** | Padrões de co-expressão entre genes/amostras; agrupamento hierárquico revela clusters biológicos | Para os top-N genes mais variáveis, ou para um conjunto curado (ex.: família CYP450) | `ComplexHeatmap` (Gu, 2016) ou `pheatmap` (mais simples, menos customizável) |
| **Volcano plot** | Visão de conjunto de todos os genes testados: magnitude (log2FC, eixo x) vs. significância (−log10 padj, eixo y) num único gráfico | Painel-chave de Results de qualquer análise de expressão diferencial | `EnhancedVolcano` (Blighe, Rana & Lewis) |
| **MA-plot** | Relação entre expressão média (baseMean, eixo x) e log2FC (eixo y) — expõe viés dependente de expressão baixa | Diagnóstico complementar ao volcano; útil para justificar shrinkage de LFC (`lfcShrink`) | `plotMA()` (DESeq2) |
| **Barplot/dotplot de GO/KEGG** | Termos funcionais enriquecidos entre os DEGs, ranqueados por significância (dotplot também mostra gene ratio e contagem) | Para resumir a interpretação funcional de uma lista de dezenas/centenas de genes | `clusterProfiler` (Wu et al., 2021) |
| **Gráfico de completude BUSCO** | Proporção de ortólogos universais completos/fragmentados/ausentes — métrica-padrão de qualidade de montagem | Sempre no início de Results, como validação da montagem antes de qualquer análise downstream | script R/Python padrão do BUSCO (`generate_plot.py`) |
| **Curva ExN50** | Completude estrutural da montagem em função da profundidade de expressão — complementa BUSCO (que mede completude gênica, não estrutural) | Quando se quer justificar escolha de parâmetros de montagem/filtragem por expressão | script Trinity `TrinityStats.pl` + `ExN50` |
| **Mapa de via KEGG (pathview)** | Sobrepõe log2FC dos genes diferencialmente expressos diretamente no diagrama da via metabólica/sinalização | Quando uma via específica (ex.: metabolismo de xenobióticos via P450) aparece enriquecida e merece destaque visual | `pathview` (Luo & Brouwer, 2013) |

### Padrões de qualidade de figura

- **Resolução:** exportar em ≥300 dpi para raster (PNG/TIFF); jamais usar captura de tela.
- **Formato vetorial:** sempre que possível, exportar também em SVG ou PDF (`ggsave(..., device = "pdf")` ou `device = "svg"`). O vetorial permite edição fina no Illustrator/Inkscape sem perda de qualidade, e boa parte das revistas exige justamente esse formato para a arte final.
- **Paletas acessíveis para daltonismo:** usar `viridis`, `cividis` ou paletas colorblind-safe do `RColorBrewer` (`brewer.pal(n, "RdYlBu")` com `colorblindr::cvd_grid()` para checagem). Evitar o clássico vermelho-verde em heatmaps e volcano plots.
- **Fontes legíveis:** tamanho mínimo de texto em eixos/legendas equivalente a 6-7pt no tamanho final impresso; usar fontes sans-serif (Arial/Helvetica). A maioria das revistas exige isso explicitamente nas instruções aos autores.
- **Consistência entre figuras:** mesma paleta de cores para as mesmas condições/tecidos em todas as figuras do paper (ex.: sempre laranja = tratado, azul = controle).

## Interpretação biológica

Uma lista de DEGs ou termos GO/KEGG enriquecidos não é, por si só, uma conclusão biológica: é um dado intermediário. Interpretá-la é um processo em funil.

1. **Da lista de DEGs à anotação funcional.** Cruzar os IDs de transcritos/genes com a anotação (BLAST/eggNOG/InterProScan) para obter nomes de genes, domínios Pfam e termos GO associados.
2. **Do enriquecimento estatístico à via biológica.** O `clusterProfiler`/`pathview` aponta *quais* vias estão sobrerrepresentadas — mas o pesquisador precisa perguntar "essa via faz sentido para a biologia do meu organismo e da minha pergunta experimental?"
3. **Da via biológica à narrativa do organismo-praga.** Para insetos-praga como *A. gemmatalis* e *M. spectabilis*, a pergunta costuma ser: *o que essa expressão diferencial nos diz sobre defesa/detoxificação/resistência a inseticidas?* Genes-chave a monitorar nessa narrativa:
   - **Citocromo P450 (CYP450):** superfamília de monoxigenases envolvida na fase I de detoxificação de xenobióticos (inseticidas, aleloquímicos da planta hospedeira). Superexpressão de subfamílias específicas (ex.: CYP6, CYP9 em Lepidoptera) é uma assinatura clássica de resistência metabólica a inseticidas.
   - **Glutationa-S-transferase (GST):** fase II de detoxificação, conjuga glutationa a eletrófilos reativos (incluindo metabólitos de piretroides e organofosforados). Enriquecimento de GST em vias KEGG de "metabolismo de xenobióticos por citocromo P450" ou "glutationa" reforça a hipótese de resposta a estresse químico.
   - **UDP-glicosiltransferase (UGT):** fase II, glicosila compostos lipofílicos (incluindo aleloquímicos vegetais e inseticidas) para facilitar a excreção. Em pragas associadas a plantas hospedeiras específicas (soja para *A. gemmatalis*, gramíneas para *M. spectabilis*), UGTs costumam estar entre os genes mais expressos em tecidos de alimentação e detoxificação.
   - **Genes de imunidade/resposta a estresse:** peptídeos antimicrobianos, vias Toll/IMD, proteínas de choque térmico — relevantes quando o desenho experimental envolve infecção, endossimbiontes ou estresse abiótico.

**Regra de ouro:** a narrativa biológica só é defensável se (a) o termo GO/KEGG é estatisticamente robusto (padj ≤ 0.05 no teste de enriquecimento, não apenas nos DEGs individuais) e (b) há suporte de anotação funcional direta (domínio Pfam correspondente, não apenas homologia BLAST de baixa cobertura). Evitar generalizações do tipo "achamos genes de resistência a inseticidas" quando na verdade se encontrou uma CYP450 sem subfamília identificada. Especifique sempre que der.

## Reprodutibilidade computacional

Reprodutibilidade não é burocracia — é o que permite que outro laboratório (ou você mesmo, em 2 anos) reexecute a análise e obtenha o mesmo resultado. Camadas de reprodutibilidade, da mais básica à mais robusta:

1. **Registro de versões e parâmetros.** Toda ferramenta usada (Trinity, STAR, DESeq2, clusterProfiler etc.) deve ter sua versão exata registrada, junto com todos os parâmetros não-default e qualquer seed de aleatoriedade (ex.: `set.seed()` no R antes de qualquer amostragem/bootstrap). Isso vai na seção Methods do paper.
2. **Ambiente conda versionado.** Um arquivo `environment.yml` (gerado via `mamba env export --no-builds > environment.yml`) capturando todas as dependências e suas versões, versionado junto ao código.
3. **Containers (Docker/Singularity/Apptainer).** Um passo acima do conda: empacota não só as dependências Python/R, mas todo o sistema operacional subjacente, eliminando problemas de "funciona na minha máquina". Em ambiente acadêmico de HPC, **Singularity/Apptainer** é preferido a Docker por não exigir privilégios de root, o que importa em clusters institucionais como os da UFV/CENAPAD.
4. **Orquestração versionada (Nextflow/Snakemake).** Em vez de uma sequência de comandos soltos em um README, o pipeline inteiro é codificado como um workflow versionado, com cada processo isolado (frequentemente já contendo sua própria definição de container por etapa). Isso é o que hoje se aproxima de "padrão-ouro" em transcriptômica de novo (ver seção seguinte).
5. **Versionamento de código (GitHub) + arquivamento com DOI (Zenodo).** O repositório GitHub garante histórico e colaboração, mas não é uma referência citável permanente: um commit pode ser reescrito, uma branch deletada, uma conta suspensa. O passo final é arquivar a versão usada no paper no **Zenodo** (a integração nativa GitHub→Zenodo gera um DOI imutável) e citar esse DOI na seção de disponibilidade de dados/código do manuscrito.

## Prática usual vs. estado da arte (2024-2026)

| Aspecto | Prática usual/consolidada | Estado da arte 2024-2026 | Fonte (DOI/URL oficial) |
|---|---|---|---|
| Orquestração do pipeline de montagem de novo | Sequência de comandos bash documentada em README, ambiente conda `environment.yml` listado no Methods | Pipeline Nextflow versionado com containers por processo (Docker/Singularity), ex. `nf-core/denovotranscript` (2025) e seu antecessor `TransPi` (2022), ambos publicados especificamente para eucariotos não-modelo | TransPi: Rivera-Vicéns et al. 2022, *Mol Ecol Resour* 22:2070-2086, [doi:10.1111/1755-0998.13593](https://doi.org/10.1111/1755-0998.13593); nf-core/denovotranscript: Bhojwani et al. 2025, *Current Protocols*, [doi:10.1002/cpz1.70248](https://doi.org/10.1002/cpz1.70248) |
| Triagem de contaminação antes do depósito TSA | Inspeção manual/BLAST pontual contra bancos de contaminantes conhecidos | Uso obrigatório recomendado do **FCS-GX** (Foreign Contamination Screen) do NCBI antes de qualquer submissão, para evitar atraso/rejeição pós-submissão | NCBI TSA Submission Guide, atualizado out/2024, [ncbi.nlm.nih.gov/genbank/tsaguide](https://www.ncbi.nlm.nih.gov/genbank/tsaguide/) |
| Avaliação de qualidade do software/pipeline em si | Nenhuma métrica formal — "funciona" é o critério | Frameworks quantitativos de aderência aos princípios FAIR aplicados a *software* de pesquisa (não só dados), com indicadores automatizáveis via plataformas como OpenEBench | Serrano-Antón et al. (FAIRsoft), 2024, *Bioinformatics* 40(8):btae464, [doi:10.1093/bioinformatics/btae464](https://doi.org/10.1093/bioinformatics/btae464) |
| Disponibilidade de código citável | Link para repositório GitHub na seção de Methods (sem garantia de persistência) | Exigência crescente (não universal, mas cada vez mais comum em revistas como GigaScience) de arquivamento com DOI permanente (Zenodo) e licença open-source explícita no momento da submissão | GigaScience — política de dados FAIR e código open-source, [academic.oup.com/gigascience](https://academic.oup.com/gigascience); FAIR Guiding Principles, Wilkinson et al. 2016, *Sci Data* 3:160018, [doi:10.1038/sdata.2016.18](https://doi.org/10.1038/sdata.2016.18) |

**Nota sobre a pergunta "containers/Nextflow tornaram-se mandatórios para submissão?"** A busca não encontrou evidência de que revistas como *Genome Biology*, *GigaScience* ou *BMC Bioinformatics* exijam formalmente, como condição de submissão, que o pipeline seja entregue como container ou workflow Nextflow/Snakemake versionado. O que existe é uma tendência clara: (a) mandato de disponibilização de código sob licença aberta é praticamente universal nessas revistas; (b) o ecossistema `nf-core` virou o padrão de fato para pipelines de bioinformática publicados, inclusive com um pipeline oficial dedicado a montagem de novo de transcriptomas lançado em 2025; (c) revisores cada vez mais pedem, informalmente, que o pipeline seja reexecutável via container. **Recomendação para o curso:** trate `environment.yml` + registro de versões/parâmetros como o mínimo aceitável, e Nextflow/Snakemake + container como a prática recomendada e crescentemente esperada por revisores. Não dá para afirmar que seja uma exigência formal universal — isso seria impreciso.

## Submissão de dados (TSA/SRA) e escrita IMRaD

### Checklist de depósito de dados

1. **SRA primeiro, sempre.** As leituras brutas (FASTQ) devem ser depositadas no SRA *antes* da submissão da montagem — o TSA exige o número de acesso SRR, mais BioProject (PRJNA) e BioSample (SAMN) já existentes.
2. **Triagem de contaminação (FCS-GX)** rodada sobre a montagem final antes do envio.
3. **Regras de sequência para TSA:** comprimento mínimo de 200 pb; sem "N" nas extremidades da sequência; blocos internos de "N" com mais de 14 nucleotídeos exigem anotação de *assembly gap feature* (o portal de submissão orienta essa etapa).
4. **Formato de arquivo:** FASTA simples (extensão `.fsa`, com `moltype=transcribed_RNA`) se não houver anotação funcional incorporada; arquivo ASN.1/`.sqn` gerado via `table2asn` se houver anotação (nomes de proteína devem seguir as *International Protein Nomenclature Guidelines*).
5. **Submissão via TSA Submission Portal** do NCBI (não mais por e-mail/FTP direto).
6. **Código e pipeline no GitHub**, com release final arquivada no Zenodo (DOI citável no manuscrito).

### Escrita em formato IMRaD

- **Introduction:** contextualiza a lacuna (organismo não-modelo sem genoma de referência, relevância como praga agrícola) e formula a pergunta biológica.
- **Methods:** esta é a seção que precisa ser *literalmente reproduzível* — todo comando com seus parâmetros não-default, toda versão de ferramenta, ambiente conda/container referenciado, seeds de aleatoriedade, e os números de acesso SRA/BioProject/TSA assim que disponíveis. Frases-modelo úteis:
  - *"A montagem de novo do transcriptoma foi realizada com Trinity v2.15.1 (Grabherr et al., 2011), utilizando os parâmetros `--max_memory 30G --CPU 16 --seqType fq`, a partir de N bibliotecas de RNA-Seq pareadas (2×150 pb)."*
  - *"Transcritos redundantes foram agrupados com CD-HIT-EST v4.8.1 a 95% de identidade nucleotídica."*
  - *"A completude da montagem foi avaliada com BUSCO v5.x contra o conjunto de ortólogos `insecta_odb10`."*
  - *"A expressão diferencial foi estimada com DESeq2 v1.4x, considerando genes diferencialmente expressos aqueles com |log2FoldChange| ≥ 1 e padj ≤ 0,05 (teste de Wald, correção de Benjamini-Hochberg)."*
- **Results:** apresenta as figuras-chave (PCA, volcano, heatmap, GO/KEGG, BUSCO) com legendas autoexplicativas — uma legenda deve permitir que o leitor entenda a figura sem voltar ao texto principal (n de amostras, teste estatístico usado, significado de cores/símbolos).
- **Discussion:** conecta os achados de volta à biologia do organismo-praga e à literatura prévia — é aqui que a narrativa de detoxificação/resistência a inseticidas se desenvolve com apoio de citações.

## Checklist final "do dado ao paper"

- [ ] Leituras brutas depositadas no SRA, com BioProject e BioSample criados
- [ ] Montagem final passou por triagem de contaminação (FCS-GX)
- [ ] BUSCO, ExN50 e demais métricas de qualidade da montagem documentados e reportados em figura/tabela
- [ ] Todas as figuras exportadas em ≥300 dpi (raster) e/ou SVG/PDF (vetorial), com paleta colorblind-safe
- [ ] Lista de DEGs, termos GO/KEGG enriquecidos e vias KEGG anotadas revisadas quanto a plausibilidade biológica (não apenas significância estatística)
- [ ] `environment.yml` (ou container Docker/Singularity) versionado e disponível
- [ ] Pipeline completo versionado no GitHub, idealmente como workflow Nextflow/Snakemake
- [ ] Release do repositório arquivada no Zenodo com DOI
- [ ] Todas as versões de software e parâmetros não-default documentados na seção Methods
- [ ] Seeds de aleatoriedade registrados onde aplicável
- [ ] Montagem submetida ao TSA via TSA Submission Portal, seguindo regras de comprimento mínimo/N/anotação
- [ ] Declaração de disponibilidade de dados e código (Data/Code Availability statement) completa no manuscrito, com todos os DOIs e números de acesso

## Aplicação aos organismos-caso

**Anticarsia gemmatalis** (lagarta-da-soja, Lepidoptera): ao final do pipeline, uma narrativa plausível seria esta. A montagem de novo (validada por BUSCO contra `insecta_odb10` e ExN50) revelou, entre os genes diferencialmente expressos em tecido de intestino médio exposto a um inseticida ou a um genótipo de soja resistente, enriquecimento significativo de termos GO relacionados a "processo metabólico de xenobiótico" e da via KEGG "metabolismo de xenobióticos por citocromo P450". O heatmap dos genes dessa via mostraria subfamílias específicas de CYP450 (ex.: CYP6/CYP9) e GSTs superexpressas no grupo exposto, enquanto o mapa `pathview` da via correspondente tornaria visualmente evidente quais enzimas específicas da cascata de detoxificação (fase I → fase II → transporte/excreção) estão sendo mobilizadas. É esse tipo de detalhe que conecta a figura à pergunta prática de manejo de resistência a inseticidas em lavouras de soja.

**Mahanarva spectabilis** (cigarrinha-das-pastagens, Hemiptera): como o transcriptoma é de novo e sem genoma de referência publicado, a completude BUSCO e a anotação funcional (eggNOG/Pfam) carregam peso extra na credibilidade do estudo. Dado o histórico do inseto de se alimentar de seiva de gramíneas e abrigar endossimbiontes, uma narrativa esperada seria: os DEGs relacionados a glândula salivar mostram enriquecimento de GO ligados a resposta a estresse oxidativo e desintoxicação de aleloquímicos da planta hospedeira (UGTs, GSTs), enquanto uma fração dos transcritos anotados como de origem bacteriana (já discutida no módulo de anotação) reforça a hipótese de contribuição metabólica de endossimbiontes (padrão consistente com relatos de *Sulcia*/*Sodalis*-like em Auchenorrhyncha). Na Discussion, isso deve ser tratado como hipótese a validar com dados adicionais (ex.: FISH, metagenômica direcionada), nunca como conclusão definitiva a partir de RNA-Seq isolado.

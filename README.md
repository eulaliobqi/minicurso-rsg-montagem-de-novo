# Minicurso: Montagem de Transcriptoma de Eucariotos Não-Modelo via RNA-Seq

Minicurso de 6 horas sobre montagem *de novo* de transcriptomas, para quem trabalha com espécies **sem genoma de referência publicado** — o cenário mais comum na entomologia agrícola e em biodiversidade em geral, mas raramente coberto pelos tutoriais padrão de RNA-Seq (que assumem um genoma-modelo bem anotado).

O curso usa dois organismos-caso reais, ambos sem referência de qualidade disponível:

| Organismo | Nome popular | Grupo | Papel no curso |
|---|---|---|---|
| *Anticarsia gemmatalis* | Lagarta-da-soja | Lepidoptera | Herbívoro foliar; pragas de detoxificação (P450, esterases) |
| *Mahanarva spectabilis* | Cigarrinha-das-pastagens | Hemiptera | Sugador de floema; associado a endossimbiontes obrigatórios (Sulcia, Sodalis-like) |

## Para quem é este curso

Pós-graduandos e pesquisadores que já sabem o que é RNA-Seq, mas nunca precisaram montar um transcriptoma do zero porque seu organismo de interesse não tem genoma publicado. Não é preciso experiência prévia com Trinity, rnaSPAdes ou linha de comando avançada — os conceitos são construídos do zero em cada módulo.

## Estrutura do repositório

```
Minicurso-rsg/
├── material-teorico/     # conteúdo teórico dos 8 módulos (Markdown)
├── Script/                # scripts do pipeline, versão genérica e comentada
├── Artigos/               # artigos de referência sobre montagem de novo em não-modelo
└── Apresentação/          # slides finais em PowerPoint (.pptx)
```

> A fonte editável dos slides (deck Reveal.js usado para gerar o `.pptx`) não faz parte deste repositório.

## Módulos

O curso segue o fluxo padrão de um pipeline de RNA-Seq *de novo* — cada módulo resolve o problema deixado pelo anterior:

| # | Módulo | Duração | Conteúdo |
|---|---|---|---|
| [0](material-teorico/00-fundamentos.md) | Fundamentos e Desenho Experimental | — | O que é transcriptoma, organismo não-modelo, genome-guided vs. *de novo*, desenho experimental |
| [1](material-teorico/01-qualidade-preprocessamento.md) | Qualidade e Pré-processamento das Reads | — | FastQC/MultiQC, trimming, remoção de contaminantes e rRNA |
| [2](material-teorico/02-montagem-de-novo.md) | Montagem de Novo | 60 min | Grafo de De Bruijn, Trinity, rnaSPAdes, escolha de *k* |
| [3](material-teorico/03-avaliacao-montagem.md) | Redução de Redundância e Avaliação da Montagem | 50 min | CD-HIT-EST, EvidentialGene, métricas de qualidade (BUSCO, N50) |
| [4](material-teorico/04-anotacao-funcional.md) | Predição de ORFs e Anotação Funcional | 60 min | TransDecoder, busca por homologia, GO/KEGG/Pfam, filtragem de contaminantes |
| [5](material-teorico/05-quantificacao-de-enriquecimento.md) | Quantificação, Expressão Diferencial e Enriquecimento | 55 min | Salmon, agregação transcrito→gene (Corset/tximport), DESeq2, clusterProfiler |
| [6](material-teorico/06-figuras-interpretacao-paper.md) | Figuras, Interpretação Biológica, Reprodutibilidade e Escrita do Paper | 60 min | Figuras publicáveis (PCA, heatmap, volcano, MA-plot), reprodutibilidade, submissão |
| [6b](material-teorico/06b-avancado-ortologia-filogenia.md) | *(Avançado, extraclasse)* Ortologia e Filogenia | — | Ortólogos em espécie-modelo, expansão/contração de famílias gênicas, árvores filogenéticas |

Material de apoio, consultável a qualquer momento durante ou após o curso:

- **[Glossário](material-teorico/glossario.md)** — ~100 termos técnicos definidos, em ordem alfabética.
- **[FAQ](material-teorico/faq.md)** — ~50 perguntas frequentes, organizadas por módulo.
- **[Referências](material-teorico/referencias.md)** — bibliografia completa, com DOI conferido para cada citação, nos formatos ABNT e APA.

## Como usar este material

1. Siga os módulos na ordem (0 → 6), cada um assume que o anterior foi concluído.
2. O Módulo 6b é opcional e independente — só é necessário para quem quiser aprofundar em ortologia/filogenia após a anotação funcional (Módulo 4).
3. Consulte o glossário sempre que um termo técnico não estiver claro; consulte o FAQ para dúvidas comuns já respondidas.
4. Os slides finais (`Apresentação/apresentacao-minicurso-transcriptoma.pptx`) seguem a mesma estrutura de módulos e podem ser usados como roteiro de aula.

## Scripts e pipelines

**[`Script/scripts.md`](Script/scripts.md)** reúne, em um único documento comentado, a lógica central de cada etapa do pipeline usado no curso — versão generalizada e resumida, pensada para leitura em aula:

1. Controle de qualidade e trimagem de reads (FastQC, fastp, MultiQC)
2. Montagem *de novo* (Trinity, CD-HIT-EST, TransDecoder)
3. Avaliação de qualidade da montagem (TrinityStats, BUSCO)
4. Anotação funcional (DIAMOND, TaxonKit, eggNOG-mapper, Pfam/HMMER)
5. Visualização (figuras de anotação, QC de montagem, GO/KEGG)
6. Triagem de endossimbiontes a partir da anotação já gerada
7. Predição de secretoma (SignalP + TMHMM) — etapa pendente de execução

Cada seção aponta o script de produção completo correspondente (com tratamento de erro, retomada de execução e logging), além de duas lições reais do projeto que valem a leitura: um bug de parsing por índice de coluna hardcoded no eggNOG-mapper, e o critério de corte usado para secretoma clássico.

## Artigos de referência

A pasta **[`Artigos/`](Artigos/)** reúne quatro artigos usados como referência para o desenho do pipeline, todos aplicando montagem *de novo* de transcriptoma a organismos não-modelo:

| Arquivo | Referência |
|---|---|
| [`Assadi-2022.pdf`](Artigos/Assadi-2022.pdf) | Sheikh-Assadi, M. et al. Normalized Workflow to Optimize Hybrid De Novo Transcriptome Assembly for Non-Model Species: A Case Study in *Lilium ledebourii*. *Plants*, 2022. |
| [`Duarte-2021.pdf`](Artigos/Duarte-2021.pdf) | Duarte, G. T.; Volkova, P. Yu.; Geras'kin, S. A. A Pipeline for Non-model Organisms for de novo Transcriptome Assembly, Annotation, and Gene Ontology Analysis Using Open Tools: Case Study with Scots Pine. *Bio-protocol*, 11(03): e3912, 2021. |
| [`Joshi-2024.pdf`](Artigos/Joshi-2024.pdf) | Joshi, J.; Wang-Pruski, G. De novo transcriptome assembly and differential gene expression analysis in different developmental stages of *Agriotes sputator* (click beetle). *Scientific Reports*, 2024. |
| [`Santillan-2019.pdf`](Artigos/Santillan-2019.pdf) | Moreno-Santillán, D. D. et al. De Novo Transcriptome Assembly and Functional Annotation in Five Species of Bats. *Scientific Reports*, 2019. |

## Licença e uso

Material didático produzido para fins educacionais. Fotos dos organismos-caso usadas na apresentação têm licença Creative Commons com atribuição (Wikimedia Commons, iNaturalist).

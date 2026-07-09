# Módulo 2 — Montagem de Novo

*(60 min — módulo mais denso do minicurso)*

## Conceitos-chave

### Grafo de De Bruijn

Diferente da montagem de genomas guiada por sobreposição de reads (overlap-layout-consensus), os montadores de transcriptoma de novo usados neste curso (Trinity, rnaSPAdes) constroem um **grafo de De Bruijn**: cada read é quebrada em todos os seus k-mers (subsequências de comprimento fixo *k*), e o grafo conecta k-mers que se sobrepõem em *k*-1 bases. Caminhar por esse grafo reconstrói sequências mais longas sem nunca alinhar um read contra outro diretamente, o que barateia bastante a montagem de milhões de reads curtos.

A escolha de *k* é um trade-off clássico: k pequeno (ex. 25) conecta mais fragmentos mas confunde regiões repetitivas/parálogas; k grande (ex. 55) resolve melhor repetições mas fragmenta transcritos de baixa cobertura. Trinity usa k=25 fixo (não é parâmetro do usuário); rnaSPAdes usa uma abordagem **multi-k-mer**, testando vários valores de k e combinando os resultados. É por isso que costuma superar Trinity em transcritos de baixa cobertura (ver seção de comparação).

### Por que transcriptoma é mais difícil que genoma

Em RNA-Seq, a cobertura de um "nó" do grafo é proporcional ao nível de expressão do gene, não ao seu tamanho — genes altamente expressos geram milhões de k-mers idênticos, genes de baixa expressão geram poucos. Isso quebra a suposição de cobertura uniforme que sustenta montadores de genoma, e é a razão de existirem ferramentas dedicadas (Trinity, rnaSPAdes) em vez de reaproveitar montadores genômicos (SPAdes, Velvet) diretamente.

### Ambiguidade de isoformas e o que Trinity chama de "gene"

Um mesmo locus genético pode gerar múltiplas isoformas por splicing alternativo, e diferentes genes parálogos compartilham éxons quase idênticos. No grafo de De Bruijn, isso aparece como bifurcações (bubbles) dentro do mesmo componente conectado. Trinity resolve isso em duas etapas (Chrysalis + Butterfly, ver abaixo), mas vale entender a nomenclatura de saída antes de seguir:

```
>TRINITY_DN1000_c115_g5_i1
              │    │   │  │
              │    │   │  └─ isoform (i1, i2, ...)
              │    │   └──── "gene" Trinity = cluster de isoformas
              │    └──────── component (subgrafo pós-Chrysalis)
              └───────────── ID do De Bruijn graph (DN = "de novo")
```

**O "gene" (`g#`) do Trinity NÃO é um gene biológico confirmado**: é um agrupamento de isoformas que compartilham sequência/estrutura de grafo suficiente para serem consideradas do mesmo locus, segundo a heurística do Chrysalis. Na prática:
- Genes parálogos muito similares podem ser fundidos em um único "gene" Trinity (superestimando isoformas, subestimando genes).
- Genes verdadeiramente distintos com baixa expressão podem aparecer fragmentados em múltiplos componentes (subestimando o gene como um todo).

Por isso, a contagem "42.372 transcritos" (*A. gemmatalis*) ou "90.344 genes / 103.560 transcritos" (*M. spectabilis*) deve sempre ser lida como **unidades de montagem**, não como número definitivo de genes biológicos. A validação real vem da anotação funcional (Módulo 3) e, idealmente, de comparação com um genoma de referência, quando disponível.

---

## Ferramentas e fluxos

### Trinity — Inchworm → Chrysalis → Butterfly

Trinity (Grabherr et al. 2011, *Nat Biotechnol* 29:644–652, DOI [10.1038/nbt.1883](https://doi.org/10.1038/nbt.1883)) é composto por três módulos sequenciais:

1. **Inchworm** — monta os k-mers (k=25) em contigs lineares gulosos (greedy), gerando os transcritos mais prováveis sem considerar ainda isoformas alternativas. É a etapa mais rápida e menos exigente em memória.
2. **Chrysalis** — agrupa os contigs do Inchworm que compartilham k-mers em **componentes** (clusters), constrói um grafo de De Bruijn completo para cada componente e particiona os reads de volta entre eles. É a etapa mais exigente em RAM.
3. **Butterfly** — processa cada grafo de componente individualmente, usando os reads (incluindo informação de pareamento) para resolver caminhos no grafo que correspondem a isoformas e variantes alélicas/parálogas, reportando as sequências finais no formato `c#_g#_i#` acima.

**Parâmetros-chave:**

| Parâmetro | Função | Observação prática |
|---|---|---|
| `--max_memory 30G` | Teto de RAM para Chrysalis/jellyfish | Ver regra de dimensionamento abaixo |
| `--CPU 16` | Threads paralelas | Usar `--CPU 16` no servidor Debian (32 cores, deixar margem) |
| `--seqType fq` | Formato de entrada | FASTQ é o padrão em RNA-Seq Illumina |
| `--min_kmer_cov` | Cobertura mínima de k-mer para entrar no Inchworm | Default=1; aumentar (ex. 2) em datasets muito profundos para reduzir ruído/tempo, à custa de sensibilidade em transcritos raros |
| `--SS_lib_type` | Biblioteca stranded (`RF` ou `FR`) | Se o kit de biblioteca preservar strand (ex. dUTP), declarar evita quimeras antisense; **omitir** se a biblioteca não for stranded |

Comando de referência usado no projeto (`--max_memory 30G --CPU 16 --seqType fq`, conforme padrão do laboratório).

### rnaSPAdes — abordagem multi-k-mer

rnaSPAdes (Bushmanova et al. 2019, *GigaScience* 8(9):giz100, DOI [10.1093/gigascience/giz100](https://doi.org/10.1093/gigascience/giz100)) reaproveita o motor SPAdes, mas com um algoritmo de simplificação de grafo dedicado a transcriptomas (remoção de bolhas por diferença de cobertura, não só topologia). Diferente do k fixo do Trinity, ele constrói e combina automaticamente grafos com múltiplos valores de k. Na prática, isso melhora a montagem de transcritos de **baixa expressão** (menos k-mers, mais sensível a k pequeno), tende a produzir **menos transcritos totais e com N50 mais alto** que Trinity no mesmo dataset (à custa de fragmentar mais isoformas raras) e consome tipicamente **menos RAM**, o que o torna mais tratável em notebooks e servidores modestos.

Benchmarks do próprio artigo original, replicados por trabalhos posteriores (ex. transXpress, ver abaixo), mostram que **nenhum dos dois é estritamente superior**: rnaSPAdes tende a vencer em contiguidade/N50 e uso de recursos, Trinity tende a vencer em sensibilidade a isoformas raras e segue mais testado em artrópodes/insetos, por ser o padrão histórico da literatura entomológica.

### Pipelines multi-montador (estado da arte para completude)

A resposta consolidada na literatura recente para maximizar completude e reduzir o viés de um único montador é **não escolher um só montador**, e sim combinar vários e depois consolidar. As principais implementações hoje:

- **TransPi** (Rivera-Vicéns et al. 2022, *Mol Ecol Resour* 22(5):2070–2086, DOI [10.1111/1755-0998.13593](https://doi.org/10.1111/1755-0998.13593)) — pipeline Nextflow que roda automaticamente múltiplos montadores (Trinity, rnaSPAdes, entre outros com diferentes k), consolida com EvidentialGene e reporta métricas BUSCO/TransRate padronizadas. Reduz duplicação e aumenta completude BUSCO frente a montadores isolados, segundo o próprio artigo.
- **transXpress** (BMC Bioinformatics 24:133, 2023, DOI [10.1186/s12859-023-05254-8](https://doi.org/10.1186/s12859-023-05254-8)) — pipeline Snakemake que combina Trinity + rnaSPAdes e paraleliza em clusters heterogêneos, com foco em anotação integrada (SequenceServer) além da montagem.
- **Oyster River Protocol / ORP** (MacManes 2018, *PeerJ* 6:e5428, DOI [10.7717/peerj.5428](https://doi.org/10.7717/peerj.5428)) — combina Trinity, SPAdes e Shannon com múltiplos k, funde os transcritomas e extrai grupos de ortólogos; demonstrou escores de TransRate/DETONATE superiores a qualquer montador isolado no benchmark original.
- **HPC-T-Assembly** (Liberati et al. 2025, *BMC Bioinformatics* 26:113, DOI [10.1186/s12859-025-06121-4](https://doi.org/10.1186/s12859-025-06121-4)) — **citação do rascunho CONFIRMADA como real** (não é fabricada). Pipeline voltado a HPC/SLURM para montagem de novo de datasets multi-espécie em grande escala, focado em orquestração de recursos computacionais mais do que em algoritmo de montagem novo.

### Consolidação — EvidentialGene (tr2aacds)

Depois de gerar transcritos de múltiplos montadores/k, é preciso um "oráculo" que escolha o melhor conjunto não-redundante. O **EvidentialGene tr2aacds** (Gilbert, Indiana University), ferramenta de longa data na comunidade sem um artigo peer-reviewed único e canônico (citada primariamente via documentação/website do projeto), prediz ORFs, remove fragmentos e redundância por identidade de sequência e potencial codificante, e classifica transcritos em "primary" (representante do locus) ou "alternative" (isoformas adicionais). É o consolidador usado dentro do TransPi e a etapa padrão recomendada sempre que há mais de uma montagem a integrar. O Módulo 3 aprofunda esse pós-processamento e a avaliação de qualidade da montagem.

### Fronteira — long-read reference-free (PacBio Iso-Seq / Oxford Nanopore)

Para quem tem dados de long-read (não é o caso dos dois organismos-caso deste curso, que usaram apenas Illumina short-read), o campo já não depende só de montagem tipo De Bruijn:

- **RNA-Bloom2** (Nip et al. 2023, *Nat Commun*, DOI [10.1038/s41467-023-38553-y](https://doi.org/10.1038/s41467-023-38553-y)) — montagem reference-free de long-reads, usa 27,0–80,6% do pico de memória e 3,6–10,8% do tempo de execução de métodos concorrentes reference-free equivalentes, segundo o artigo original.
- **isONform** (Bioinformatics 39(Supplement_1):i222, 2023) — reconstrução de isoformas reference-free a partir de dados Oxford Nanopore, com recall substancialmente maior que RATTLE em cenários de muitas isoformas verdadeiras por locus.
- **RATTLE** — clustering/coleção de isoformas de long-read reference-free; comparações recentes (2024–2025) mostram que **perde precisão de recall conforme o número de isoformas verdadeiras aumenta**, sendo hoje considerado inferior a isONform/RNA-Bloom2 nesse quesito.
- **IsoQuant** — em benchmarks recentes desempenha-se bem tanto em modo reference-based quanto reference-free, ao lado de Bambu e StringTie2.
- **Iso-Seq3 (cluster/collapse)** — pipeline oficial da PacBio para clustering/colapso de reads full-length (FLNC) em isoformas; ainda é a via "padrão de fabricante" para quem tem dados PacBio puros.

O **LRGASP Consortium** (Nature Methods 21:1349–1363, julho de 2024, DOI [10.1038/s41592-024-02298-3](https://doi.org/10.1038/s41592-024-02298-3)) publicou o benchmark mais abrangente até hoje comparando esses métodos, usando mais de 427 milhões de long-reads (cDNA e RNA direta) em humano, camundongo e peixe-boi. **Confirmação da citação do rascunho:** o DOI e o ano (2024, Nature Methods) estão corretos. Achado central relevante para não-modelo: bibliotecas com reads mais longos e mais precisos geram transcritos mais corretos que aumentar apenas a profundidade; já profundidade maior melhora principalmente a quantificação, não a identificação de isoformas novas.

---

## Prática usual vs. estado da arte (2024-2026)

| Cenário | Ferramenta/prática usual | Estado da arte 2024-2026 | Fonte (DOI) |
|---|---|---|---|
| Montagem de transcriptoma de inseto/artrópode não-modelo (short-read) | Trinity isolado (é o que foi usado nos dois organismos-caso deste curso) | Pipeline multi-montador (TransPi/ORP/transXpress) + consolidação EvidentialGene; benchmarks mostram maior completude BUSCO e menor duplicação que montador único | TransPi: 10.1111/1755-0998.13593; ORP: 10.7717/peerj.5428 |
| Comparação de montadores em artrópodes não-modelo | Assume-se Trinity "padrão ouro" por tradição da literatura entomológica | Estudo de 2025 (bioRxiv) comparando Trinity, rnaSPAdes e IDBA-tran em inseto + crustáceo de água doce concluiu que **IDBA-tran é inadequado**, reforçando que a escolha entre Trinity/rnaSPAdes deve ser empírica (BUSCO/N50) por dataset, não assumida a priori | bioRxiv 10.1101/2025.08.01.668104 (preprint, não peer-reviewed até a data desta pesquisa) |
| Dimensionamento de RAM para Trinity | Regra de bolso "~1GB de RAM por 1 milhão de read pairs" | **Regra ainda citada e válida** na documentação oficial (Trinity Wiki — "Trinity Computing Requirements") e em múltiplas fontes de treinamento (Cornell BioHPC, Evomics) como estimativa de piso; complexidade do transcriptoma (heterozigosidade, tamanho do genoma, nº de genes) pode aumentar o requisito real bem acima da regra linear | Trinity Wiki: github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Computing-Requirements |
| Montagem long-read reference-free | Iso-Seq3 cluster/collapse (PacBio) isolado | RNA-Bloom2 e isONform superam ferramentas reference-free anteriores em memória/tempo/recall; LRGASP (2024) consolidou IsoQuant/Bambu/FLAIR como top-performers reference-based quando há genoma de qualidade | RNA-Bloom2: 10.1038/s41467-023-38553-y; LRGASP: 10.1038/s41592-024-02298-3 |
| Manutenção do Trinity | Assume-se desenvolvimento ativo contínuo | **Atenção:** financiamento dedicado do ITCR ao Trinity terminou em 2024 (financiado 2013–2024); o software continua funcional e de código aberto (última release confirmada: v2.15.2), mas o ritmo de manutenção/suporte oficial foi reduzido — relevante para decisões de adoção em pipelines de longo prazo | github.com/trinityrnaseq/trinityrnaseq (releases) |

---

## Tabela comparativa de montadores

| Montador | Contiguidade (N50) | Completude (BUSCO) | Uso de RAM | Facilidade de uso | Melhor cenário de aplicação |
|---|---|---|---|---|---|
| **Trinity** | Média-alta; sensível a isoformas | Alta em transcritos moderados/altamente expressos; pode fragmentar genes de baixa expressão | Alto (regra ~1GB/milhão de read pairs, pode escalar acima em transcriptomas complexos) | Alta — muito documentado, padrão histórico em não-modelo | Datasets com histórico de uso do grupo/comparabilidade com literatura prévia (ex. Lepidoptera, Hemiptera); quando isoformas alternativas são foco de interesse biológico |
| **rnaSPAdes** | Alta (multi-k-mer reduz fragmentação) | Comparável a Trinity, por vezes melhor em transcritos de baixa cobertura | Moderado-baixo — mais tratável em hardware modesto | Alta — linha de comando simples, poucos parâmetros obrigatórios | Datasets de RAM limitada, ou quando o objetivo é um conjunto de transcritos representativos (menos redundante) mais do que exaustivo em isoformas |
| **TransPi (multi-montador)** | Alta (herda o melhor de cada montador componente) | A mais alta entre as opções testadas nos benchmarks originais — reduz duplicação, aumenta BUSCO | Muito alto (roda vários montadores em paralelo/sequência) | Moderada — requer Nextflow + containers, mas automatiza a decisão de "qual montador usar" | Projetos de referência/publicação onde completude e reprodutibilidade justificam o custo computacional extra |
| **Oyster River Protocol (ORP)** | Alta | Alta (métricas TransRate/DETONATE superiores no artigo original) | Muito alto | Moderada-baixa — setup mais manual que TransPi | Grupos já habituados ao ecossistema MacManes-lab / quando se quer controle fino sobre cada montador componente |
| **transXpress** | Alta | Boa (Trinity+rnaSPAdes combinados) | Alto | Alta — Snakemake, cluster heterogêneo, inclui anotação | Ambientes de cluster heterogêneo (mistura de nós), quando anotação integrada (SequenceServer) é desejada junto da montagem |
| **RNA-Bloom2 (long-read)** | Muito alta (isoformas full-length nativas) | Alta, competitiva com métodos reference-based no artigo original | Baixo relativo a concorrentes long-read reference-free (27–80% do pico de memória de alternativas) | Moderada | Dados PacBio/ONT sem genoma de referência disponível |
| **IsoQuant (long-read)** | Alta (quando há genoma) | Muito alta em modo reference-based (top do benchmark LRGASP) | Moderado | Alta | Quando existe genoma de referência de qualidade, mesmo que de espécie próxima |

---

## Recomendação (TOP escolha)

**Para short-read em inseto/artrópode não-modelo sem genoma de referência (cenário dos dois organismos-caso deste curso): a recomendação de estado da arte é um pipeline multi-montador com consolidação EvidentialGene — preferencialmente TransPi (Nextflow, mais próximo do ecossistema já usado no laboratório) — em vez de Trinity isolado.**

Justificativa: os benchmarks publicados (TransPi, ORP) mostram de forma consistente que combinar Trinity + rnaSPAdes (± outros k) e consolidar com EvidentialGene produz maior completude BUSCO e menor duplicação do que qualquer montador único, sem custo adicional relevante de tempo de analista (a automação via Nextflow/Snakemake absorve a complexidade). O único trade-off real é computacional: rodar múltiplos montadores custa mais RAM e tempo de CPU que rodar um único Trinity.

**Ressalva honesta sobre este curso:** as montagens reais usadas como estudo de caso (*A. gemmatalis*, 42.372 transcritos; *M. spectabilis*, 90.344 genes/103.560 transcritos, N50=738) foram feitas com **Trinity puro**, não com um pipeline multi-montador. Isso reflete o que estava disponível/decidido nos projetos no momento da execução — não uma afirmação de que Trinity isolado seja a melhor prática atual. Um reprocessamento futuro desses mesmos dados brutos por TransPi ou ORP provavelmente aumentaria a completude BUSCO e reduziria a redundância de isoformas, ao custo de mais tempo de máquina.

Para quem já tem dados long-read (PacBio Iso-Seq ou ONT) e nenhum genoma de referência: **RNA-Bloom2** é hoje a escolha mais bem documentada e validada (Nip et al. 2023; consistente com os achados do LRGASP 2024 sobre a importância de reads longos/precisos).

---

## Otimização de recursos

- **RAM:** regra de bolso do Trinity ainda válida e citada oficialmente: **~1 GB de RAM por 1 milhão de pares de reads** para as etapas Inchworm/Chrysalis (as mais pesadas). É um piso, não um teto — transcriptomas de organismos com alta heterozigosidade, muitos parálogos ou tamanho de genoma grande podem exigir bem mais. Reserve sempre uma margem (`--max_memory` abaixo do total físico disponível) para não travar o sistema operacional do servidor.
- **CPU:** `--CPU 16` no servidor (32 cores) deixa margem para outros processos/screens simultâneos (ex. GROMACS rodando em paralelo). Butterfly paraleliza bem por componente; Chrysalis é o gargalo de RAM mais que de CPU.
- **Disco:** montagens de novo geram arquivos intermediários grandes (`.fasta`, k-mer counts, arquivos de checkpoint do Chrysalis/Butterfly). Verifique o espaço livre antes de iniciar (`df -h /home`); é uma regra do laboratório desde que `af3_databases` chegou a 96% de ocupação.
- **HPC/SLURM:** para datasets muito grandes (múltiplas amostras, multi-espécie), pipelines como HPC-T-Assembly foram desenhados especificamente para orquestrar filas SLURM e evitar que um único job monopolize o cluster. No servidor único do laboratório (32 cores, sem fila SLURM), o equivalente prático é sempre rodar dentro de `screen -S <nome>`, para não perder o processo por SIGTTOU e poder monitorar/desconectar com segurança.
- **Tempo:** regra de bolso complementar, útil para planejamento de aula/execução: ~30–60 minutos por milhão de reads nas etapas Inchworm+Chrysalis, altamente dependente da complexidade do transcriptoma. Trate como estimativa aproximada, não garantia.

---

## Aplicação aos organismos-caso

### *Anticarsia gemmatalis* (Lepidoptera — lagarta-da-soja)

Montagem Trinity real: **42.372 transcritos**. Trinity puro, sem consolidação multi-montador. Este dataset é o exemplo usado ao longo do curso para as etapas subsequentes (Módulo 3: CD-HIT, TransDecoder, anotação eggNOG/Pfam).

### *Mahanarva spectabilis* (Hemiptera — cigarrinha-das-pastagens, "amarelão")

Montagem Trinity real: **90.344 genes / 103.560 transcritos, N50 = 738**. Projeto focado em glândula salivar, com interesse adicional em toxinas e endossimbiontes bacterianos (Sulcia/Sodalis-like). O N50 relativamente baixo (738 pb) é consistente com um transcriptoma de novo que ainda não passou por consolidação/redução de redundância (EvidentialGene) nem por filtragem agressiva de fragmentos: é um candidato natural para reprocessamento com TransPi/ORP, caso o objetivo do projeto evolua para publicação com métricas de completude mais rigorosas.

**Nota metodológica honesta para os dois casos:** em nenhum dos dois projetos foi empregado rnaSPAdes, TransPi, ORP ou qualquer outro pipeline multi-montador. A escolha de Trinity refletiu o que já estava disponível e testado nos ambientes (`trinity` no servidor) no momento da execução desses projetos, não uma comparação formal de benchmarks feita pelo laboratório. Isso não invalida os resultados obtidos (ambos passaram por validação funcional a jusante: BUSCO, eggNOG, Pfam), mas é uma limitação metodológica que vale declarar em qualquer manuscrito decorrente.

---

## Fontes verificadas nesta pesquisa

| Citação | Status | DOI |
|---|---|---|
| Trinity (Grabherr et al. 2011) | Confirmado | 10.1038/nbt.1883 |
| Trinity — versão atual / status de manutenção | Confirmado (v2.15.2; financiamento ITCR encerrado em 2024) | github.com/trinityrnaseq/trinityrnaseq/releases |
| rnaSPAdes (Bushmanova et al. 2019) | Confirmado | 10.1093/gigascience/giz100 |
| TransPi (Rivera-Vicéns et al. 2022) | Confirmado | 10.1111/1755-0998.13593 |
| transXpress (2023) | Confirmado — draft do curso citava "2023 BMC Bioinformatics", correto | 10.1186/s12859-023-05254-8 |
| Oyster River Protocol (MacManes 2018) | Confirmado | 10.7717/peerj.5428 |
| RNA-Bloom2 (Nip et al. 2023) | Confirmado | 10.1038/s41467-023-38553-y |
| HPC-T-Assembly (Liberati et al. 2025) | **Confirmado como real** — a suspeita do rascunho de que fosse citação inventada NÃO se confirmou; o artigo existe, na BMC Bioinformatics, com esse DOI exato | 10.1186/s12859-025-06121-4 |
| LRGASP Consortium (2024, Nature Methods) | Confirmado | 10.1038/s41592-024-02298-3 |
| EvidentialGene tr2aacds | **NÃO HÁ DOI CANÔNICO ÚNICO** — ferramenta amplamente citada via documentação/website do autor (Don Gilbert, Indiana University), sem um artigo peer-reviewed central equivalente aos demais; citações na literatura variam (2013 em diante) | não aplicável |
| isONform | Confirmado (Bioinformatics 39(Suppl_1):i222, 2023) | 10.1093/bioinformatics/btad264 (não solicitado explicitamente, verificado por completude) |

Nenhuma citação do rascunho precisou ser removida ou marcada como incorreta: a única sinalizada como "verificar suspeita de fabricação" (HPC-T-Assembly) foi confirmada como real após busca direta na BMC Bioinformatics/PMC. A única lacuna genuína é o DOI do EvidentialGene, que não existe como citação única canônica na literatura, e isso deve ser declarado ao usar a ferramenta em qualquer manuscrito.

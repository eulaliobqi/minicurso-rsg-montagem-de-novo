# Módulo 3 — Redução de Redundância e Avaliação da Montagem

**Duração:** 50 min | **Pré-requisito:** Módulo 2 (montagem de novo concluída, ex.: `Trinity.fasta`)

## Conceitos-chave

Um assembler de novo como o Trinity monta transcritos a partir de grafos de de Bruijn construídos separadamente para cada "componente" de reads que compartilham k-mers. O resultado quase sempre contém uma fração grande de sequências redundantes: isoformas de splicing genuínas, sim, mas também fragmentos alélicos, artefatos de sobreposição entre bolhas do grafo e transcritos truncados que representam o mesmo gene várias vezes. Em *M. spectabilis*, por exemplo, o Trinity devolveu 103.560 "transcritos" agrupados em 90.344 "genes": a proporção 1,15:1 já sinaliza que a maior parte dos genes tem uma ou poucas isoformas montadas, mas a cauda de genes com dezenas de variantes quase sempre é ruído de montagem, não splicing biológico real.

Reduzir essa redundância antes de seguir para anotação e quantificação importa por dois motivos práticos. Primeiro, ferramentas de anotação (TransDecoder, eggNOG-mapper, InterProScan) cobram por sequência processada: rodar contra 103 mil transcritos quando o conjunto biológico real tem um décimo disso desperdiça tempo de CPU e, em pipelines com licença ou custo por chamada de API, dinheiro. Segundo, e mais sério: quantificadores como o Salmon dividem reads multi-mapeadas entre todas as isoformas compatíveis (modelo de verossimilhança da EM), então um inflacionamento artificial de "isoformas" quase idênticas dilui a contagem de reads entre cópias redundantes e distorce a expressão diferencial rio abaixo.

As duas estratégias consolidadas para isso são complementares, não concorrentes:

- **CD-HIT-EST** faz clustering por identidade de sequência pura (tipicamente 90–95% para transcriptomas), sem olhar cobertura de reads ou estrutura do grafo. É rápido, determinístico e fácil de justificar em métodos ("removemos redundância com CD-HIT-EST a 95% de identidade"), mas é cego a contexto biológico: pode colapsar isoformas reais que diferem por um éxon curto se a identidade global ficar acima do corte, ou deixar passar fragmentos parciais que têm baixa identidade só porque cobrem regiões diferentes do mesmo transcrito.
- **EvidentialGene** (pipeline `tr2aacds`, de Don Gilbert, Indiana University) faz algo mais sofisticado: traduz os transcritos, agrupa por similaridade de proteína e de estrutura de CDS, e classifica cada transcrito como "principal" (`okay`), "alternativo válido" (`okalt`) ou "redundante/drop" com base em sobreposição de ORF e cobertura, não só identidade bruta de nucleotídeo. O resultado costuma ser um conjunto primário mais biologicamente defensável, mas o pipeline é mais pesado computacionalmente e menos transparente: a decisão de manter ou descartar cada transcrito passa por várias heurísticas internas que a maioria dos usuários nunca lê no código-fonte.

Nenhuma das duas ferramentas tem "resposta certa" objetiva, porque não existe genoma de referência para a maioria dos organismos não-modelo trabalhados neste curso. Sem genoma, a única forma de saber se um cluster de sequências parecidas é "a mesma isoforma duas vezes" ou "duas isoformas reais" é inferência estatística sobre padrões de cobertura, não fato observável.

## Ferramentas de avaliação

Depois de reduzir redundância, a pergunta central muda de "removi o suficiente?" para "essa montagem representa bem o transcriptoma real do organismo?". Nenhuma métrica isolada responde isso, e a armadilha mais comum em quem está começando é olhar só para N50 e declarar vitória.

**Por que N50 sozinho engana:** N50 é o comprimento tal que 50% das bases montadas estão em contigs desse tamanho ou maior. Um assembler "preguiçoso" pode inflar o N50 artificialmente concatenando fragmentos não relacionados ou falhando em separar parálogos próximos em transcritos distintos: produz números bonitos e biologia errada. N50 mede tamanho, não correção nem completude. É por isso que este módulo trata N50 como ponto de partida, nunca como critério de aprovação.

**BUSCO (Benchmarking Universal Single-Copy Orthologs).** Testa se um conjunto conservado de genes ortólogos de cópia única — esperados em praticamente todo genoma de um clado — está presente e completo na montagem. A saída central é a notação C/S/D/F/M:
- **C (Complete)** — o ortólogo foi encontrado com comprimento dentro da faixa esperada.
  - **S (Single-copy)** — completo e presente uma única vez (o esperado biologicamente).
  - **D (Duplicated)** — completo, mas presente em múltiplas cópias: pode ser duplicação gênica real (alguns clados têm isso) ou sinal de que a redução de redundância (CD-HIT/EvidentialGene) não foi suficiente.
- **F (Fragmented)** — encontrado, mas com comprimento menor que o esperado: sugere montagem incompleta daquele transcrito.
- **M (Missing)** — não encontrado. Pode ser gene genuinamente ausente ou divergente no organismo, ou simplesmente não capturado pela profundidade de sequenciamento.

Um D alto (>5–8%) junto com C alto é sinal clássico de sub-redução de redundância: vale revisitar o CD-HIT-EST/EvidentialGene antes de seguir. Um M alto historicamente é interpretado como montagem incompleta, mas também pode refletir o tecido amostrado (uma biblioteca de glândula salivar, como em *M. spectabilis*, não vai expressar genes específicos de músculo, então "ausente" na montagem não significa "ausente no genoma").

**rnaQUAST.** Gera estatísticas descritivas clássicas de montagem (N50, comprimento médio, número de transcritos, GC%) e, quando há genoma de referência disponível, mapeia os transcritos de volta a ele para calcular métricas de completude e correção posicional. Para os dois organismos-caso deste curso (sem genoma de referência publicado), o rnaQUAST roda em modo de novo, usando BUSCO e GeneMarkS-T internamente; nesse cenário ele funciona mais como um agregador de relatório do que como uma métrica nova.

**TransRate.** Calcula um score por contig e um score de montagem agregado, ambos sem depender de referência: usa apenas as reads originais remapeadas contra a própria montagem para detectar quimeras, erros estruturais e bases mal suportadas. Ainda aparece em publicações de 2025, mas o desenvolvimento ativo da ferramenta parou há vários anos (o repositório oficial está essencialmente congelado); trate-a como uma métrica complementar histórica, não como padrão-ouro atual (mais detalhes na seção de estado da arte abaixo).

**ExN50 (utilitário do Trinity, `contig_ExN50_statistic.pl`).** Em vez de calcular N50 sobre todos os transcritos igualmente, o ExN50 pondera cada transcrito pela sua contribuição real à expressão (via contagem de reads do Salmon/Kallisto) e calcula o N50 considerando só os transcritos mais expressos, acumulando até uma fração `x` da expressão total (Ex). Um gráfico de Ex50 (eixo x = % da expressão acumulada, eixo y = N50 correspondente) que sobe e depois estabiliza em um platô claro indica que os transcritos biologicamente relevantes (os mais expressos) estão bem montados, mesmo que o N50 bruto, dominado por transcritos de baixíssima expressão e provável ruído, pareça mediano. É considerado superior ao N50 tradicional justamente porque penaliza menos o "lixo" de baixa cobertura que quase toda montagem de novo acumula, e foca no que realmente importa para a análise de expressão subsequente.

**Read representation / back-mapping (Bowtie2 ou Salmon).** Remapear as reads originais de volta contra a montagem e medir a porcentagem que se alinha é o proxy mais direto de completude: se 15% das reads não encontram lugar nenhum na montagem, ou esse material é contaminação (rRNA residual, DNA genômico, microbioma associado ao tecido), ou a montagem realmente perdeu uma fração relevante do transcriptoma. Valores acima de 80–90% de mapeamento (Bowtie2 com `--local` ou Salmon em modo `--validateMappings`) são geralmente considerados adequados para transcriptomas de novo; abaixo disso, vale investigar antes de prosseguir.

## Prática usual vs. estado da arte (2024-2026)

| Aspecto | Prática usual/consolidada | Estado da arte 2024-2026 | Fonte (DOI) |
|---|---|---|---|
| Redução de redundância | CD-HIT-EST a 90–95% identidade | EvidentialGene `tr2aacds` (classificação por estrutura de ORF/CDS, não só identidade) — cada vez mais citado em pipelines de novo por dar conjunto primário mais defensável biologicamente | CD-HIT: Fu et al. 2012, *Bioinformatics* 28(23):3150–3152, doi.org/10.1093/bioinformatics/bts565. EvidentialGene: Gilbert D.G. 2013, comunicação técnica (7th Annual Arthropod Genomics Symposium) — sem DOI de periódico revisado por pares até o momento |
| Software BUSCO | BUSCO v5 (Manni et al. 2021), linhagens OrthoDB v10 (`_odb10`) | **BUSCO v6.1.0** é a versão estável atual, com suporte a linhagens **OrthoDB v12.2** (`_odb12`, incluindo `insecta_odb12`, `hemiptera_odb12`, `arthropoda_odb12`, atualizadas em 03/2025). odb10 continua funcional para continuidade com literatura publicada, mas não recebe mais expansão de cobertura de genomas | Manni et al. 2021, *Mol Biol Evol* 38(10):4647–4654, doi.org/10.1093/molbev/msab199. Tegenfeldt et al. 2025, *Nucleic Acids Research* 53(D1):D516–D522, doi.org/10.1093/nar/gkae987 |
| Avaliação reference-free por contig | TransRate (score por contig + score agregado) | TransRate ainda é citado em papers de 2025, mas o repositório está sem desenvolvimento ativo há anos. **CATS** (Comprehensive Assessment of Transcript Sequences), publicado em 2026, propõe um framework reference-free (CATS-rf) e reference-based (CATS-rb) mais interpretável, com quatro componentes de score decompostos | TransRate: Smith-Unna et al. 2016, *Genome Research* 26(8):1134–1144, doi.org/10.1101/gr.196469.115. CATS: Bodulić & Vlahoviček 2026, *Nature Communications*, doi.org/10.1038/s41467-026-72171-8 |
| Estatísticas descritivas de montagem | rnaQUAST (com ou sem referência) | rnaQUAST segue mantido e ativo — release v2.3.1 em 06/06/2025 — permanece prática consolidada, não superada | Bushmanova et al. 2016, *Bioinformatics* 32(14):2210–2212, doi.org/10.1093/bioinformatics/btw218 |
| Completude ponderada por expressão | ExN50 (`contig_ExN50_statistic.pl`, Trinity) | Segue sendo a prática padrão recomendada pelos próprios autores do Trinity; nenhuma alternativa consolidada substituiu o conceito até 2026 | Haas et al. 2013, *Nature Protocols* 8:1494–1512, doi.org/10.1038/nprot.2013.084 |
| Back-mapping de reads | Bowtie2 `--local` (par de reads no par correto) | Salmon com `--validateMappings` vem substituindo Bowtie2 puro como back-mapping padrão porque já gera as contagens usadas depois no ExN50 e na quantificação — um mapeamento serve para dois propósitos | Salmon: Patro et al. 2017, *Nat Methods* 14:417–419, doi.org/10.1038/nmeth.4197. Bowtie2: Langmead & Salzberg 2012, *Nat Methods* 9:357–359, doi.org/10.1038/nmeth.1923 |

## Rubrica de decisão "aprovar / refazer"

Nenhum critério isolado decide sozinho. A tabela abaixo é o que costumo aplicar em revisão de montagem antes de liberar para anotação. Pense nela como um checklist de "sinal vermelho", não uma nota única.

| Critério | Aprovar (seguir para Módulo 4) | Sinal de alerta — investigar | Refazer/reamostrar |
|---|---|---|---|
| BUSCO C (Complete) | ≥ 90% para linhagem apropriada (ex. `insecta_odb12`) | 75–90% | < 75% |
| BUSCO D (Duplicated) dentro do C | < 5% | 5–10% (revisar redução de redundância) | > 10% (redundância não tratada) |
| ExN50 | Platô claro (Ex80–Ex90 estabiliza) | Platô tardio (>Ex95) ou instável | Sem platô — N50 cresce até Ex100 |
| Back-mapping (% reads remapeadas) | ≥ 85% | 70–85% | < 70% |
| N50 bruto | Usado só como contexto — nunca critério isolado | — | — |
| Contaminação (BUSCO de outro reino no top hits do DIAMOND/BLAST) | Ausente ou traço (<1% dos transcritos) | 1–5% de hits não-metazoários consistentes | > 5% — investigar biblioteca/preparo |
| Fragmentação (BUSCO F) | < 5% | 5–10% | > 10% — sugere profundidade insuficiente ou k-mer inadequado |

Quando os sinais são mistos (por exemplo, BUSCO C alto mas back-mapping baixo), a decisão não é binária: geralmente vale a pena investigar a causa específica (contaminação explica back-mapping baixo com BUSCO bom, porque os BUSCOs são de cópia única conservada e pouco afetados por contaminação de baixo nível) antes de decidir entre aprovar com ressalva documentada ou remontar.

## Aplicação aos organismos-caso

**Mahanarva spectabilis**: BUSCO C=96,2% na linhagem `insecta_odb10`. Pela rubrica acima, isso está confortavelmente na faixa de aprovação (≥90%), e é um número que resiste bem mesmo pensando no valor comparável em `insecta_odb12`. Dados de cobertura mais ampla do OrthoDB v12 tendem a puxar o C observado levemente para baixo (mais ortólogos de referência para checar), mas dificilmente o suficiente para tirar a montagem da faixa "boa", já que o C original era alto. O ponto de atenção real aqui não é o C, é a proporção de 103.560 transcritos para 90.344 genes: a razão de ~1,15 isoforma por gene é baixa o bastante para sugerir que a redução de redundância (CD-HIT-EST/EvidentialGene, conforme aplicado no pipeline do projeto) fez o trabalho esperado, mas vale conferir o D (Duplicated) do relatório BUSCO original antes de assumir que não sobrou redundância residual: um C=96,2% com D dentro de S não é a mesma coisa que um C=96,2% inflado por duplicações.

**Anticarsia gemmatalis**: a montagem de 42.372 transcritos não teve, até o momento de escrita deste módulo, um relatório BUSCO documentado na mesma extensão que *M. spectabilis*. Isso por si só já é um sinal para a aula: **42.372 transcritos sem métrica de completude anexada é exatamente o tipo de número que parece "razoável" a olho nu, nem absurdamente alto nem baixo para um transcriptoma de Lepidoptera, mas que não permite decisão nenhuma sobre aprovar ou refazer sem rodar BUSCO** (idealmente contra `lepidoptera_odb12`, a linhagem mais específica disponível para esse clado, seguindo a recomendação do próprio BUSCO de sempre usar o dataset mais específico possível em vez do genérico `insecta` ou `arthropoda`). Isso ilustra bem a diferença entre "ter um FASTA de saída do Trinity" e "ter uma montagem avaliada": são etapas discretas do pipeline, e é comum em projetos reais (inclusive nos do próprio laboratório) a segunda ficar pendente por mais tempo que a primeira.

---

**Leitura recomendada para o instrutor:** vale abrir ao vivo o relatório BUSCO completo (`short_summary.txt`) de *M. spectabilis* durante a aula e mostrar a decomposição C/S/D/F/M linha por linha: o número agregado "96,2%" esconde informação que só aparece quando se olha a distribuição completa.

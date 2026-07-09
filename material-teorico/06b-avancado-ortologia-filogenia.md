# Módulo 6b (Avançado) — Ortologia e Filogenia

> Módulo extraclasse, fora do cronograma cronometrado de 6h do minicurso. Destinado a quem já concluiu a montagem e anotação de um transcriptoma de novo (Módulo 4) e quer dar o próximo passo natural: comparar os genes encontrados com os de outras espécies em um contexto evolutivo explícito. Não é pré-requisito para os módulos anteriores nem exigido para entregar um projeto de RNA-Seq padrão — é aprofundamento para quem for, por exemplo, escrever um artigo com componente de evolução molecular ou investigar famílias gênicas expandidas em uma praga.

## Por que isso importa após montar um transcriptoma de novo

Depois do Módulo 4 você tem um conjunto de proteínas preditas (via TransDecoder) para *Anticarsia gemmatalis* ou *Mahanarva spectabilis*, anotadas funcionalmente por similaridade (BLAST/DIAMOND, eggNOG, Pfam). Essa anotação por similaridade responde "a que essa proteína se parece", mas não responde de forma rigorosa a três perguntas que aparecem com frequência em projetos de praga agrícola sem genoma de referência:

1. **Qual é o gene ortólogo correspondente em uma espécie modelo bem anotada** (ex.: *Bombyx mori* ou *Drosophila melanogaster* para Lepidoptera; *Nilaparvata lugens* ou *Rhodnius prolixus* para Hemiptera), para transferir com confiança uma anotação funcional detalhada (não apenas "hit de BLAST", mas relação evolutiva de ortologia 1:1)?
2. **Uma família gênica de interesse está expandida ou contraída** na espécie-praga em relação a parentes próximos? Isso é especialmente relevante para famílias de detoxificação (citocromos P450, esterases, glutationa S-transferases) e para receptores de inseticida, onde expansão gênica está associada a mecanismos de resistência.
3. **Onde a espécie se posiciona filogeneticamente** em relação a outras espécies já sequenciadas do mesmo grupo — útil tanto para contextualizar resultados quanto para escolher espécies de referência mais próximas em análises futuras (ex.: mapeamento genome-guided cross-species, Módulo 0).

Ortologia (inferir quais genes descendem de um único gene ancestral por especiação, e não por duplicação) e filogenia (reconstruir a árvore de relações evolutivas entre sequências ou espécies) são as ferramentas metodológicas para responder a essas três perguntas de forma rastreável e reprodutível, em vez de por inspeção visual de hits de BLAST.

## Conceitos-chave

**Ortogrupo.** Conjunto de genes de duas ou mais espécies que descendem de um único gene no ancestral comum de todas essas espécies. É a unidade básica de comparação em análises multi-espécies — mais robusta que "par de ortólogos" porque acomoda naturalmente duplicações e perdas gênicas dentro do grupo.

**Ortólogos vs. parálogos.** Ortólogos são genes que divergiram por um evento de **especiação** (o gene ancestral estava em uma única espécie que depois se dividiu em duas linhagens). Parálogos são genes que divergiram por um evento de **duplicação gênica** dentro da mesma linhagem. A distinção importa porque ortólogos tendem a reter função equivalente entre espécies (base da transferência de anotação), enquanto parálogos frequentemente divergem funcionalmente (neofuncionalização/subfuncionalização). É exatamente esse mecanismo de duplicação-e-divergência que gera expansões de famílias como P450 em insetos-praga.

**Árvore de genes vs. árvore de espécies.** Uma árvore de genes é construída a partir do alinhamento de sequências de um gene (ou ortogrupo) específico e reflete a história evolutiva *daquele gene*. Uma árvore de espécies reflete a história evolutiva das *linhagens*. As duas nem sempre coincidem: eventos de duplicação/perda gênica, transferência horizontal (rara em eucariotos, mas não nula) e sorting incompleto de linhagem (ILS: quando um polimorfismo ancestral persiste através de múltiplos eventos de especiação) podem fazer com que a árvore de um gene individual difira da árvore real de espécies. Por isso a prática recomendada é inferir a árvore de espécies a partir de muitos genes de cópia única (concatenação ou métodos de coalescência), não de um gene isolado.

## Ferramentas e fluxo

O fluxo padrão, partindo dos proteomas preditos de múltiplas espécies (incluindo o TransDecoder output do Módulo 4), é:

```
proteomas (FASTA, 1 isoforma representativa/gene) 
        │
        ▼
   OrthoFinder  →  ortogrupos + ortólogos par-a-par + genes de cópia única
        │
        ▼
   MAFFT (por ortogrupo de cópia única)  →  alinhamentos múltiplos
        │
        ▼
   IQ-TREE (ou FastTree)  →  árvores de genes / árvore de espécies
```

**OrthoFinder** (Emms & Kelly, 2019, *Genome Biology* 20:238, DOI: [10.1186/s13059-019-1832-y](https://doi.org/10.1186/s13059-019-1832-y)) infere ortogrupos a partir de comparações de similaridade par-a-par (DIAMOND por padrão) entre todos os proteomas de entrada, normaliza por comprimento e conteúdo de sequência para reduzir vieses entre espécies com genomas/anotações de qualidade desigual, e depois usa árvores filogenéticas (não apenas clustering por similaridade) para resolver ortólogos e parálogos dentro de cada ortogrupo, inferir eventos de duplicação e enraizar a árvore de espécies. Uso típico:

```bash
mamba activate kerson-paper   # ou ambiente equivalente com orthofinder instalado
orthofinder -f pasta_com_proteomas/ -t 32 -a 8
```

- Input: um FASTA de proteínas por espécie na mesma pasta, contendo **uma isoforma/proteína representativa por gene** — não o output bruto do TransDecoder com todas as isoformas, que infla artificialmente o número de "genes" e distorce a estimativa de ortogrupos (ver seção de limitações abaixo).
- `-t`: threads para as buscas de similaridade (BLAST/DIAMOND); `-a`: threads para a etapa de análise de árvores — em geral `-a` bem menor que `-t` é suficiente, pois essa etapa é mais leve por ortogrupo.
- Output-chave: `Orthogroups/Orthogroups.tsv` (matriz espécie × ortogrupo), `Single_Copy_Orthologue_Sequences/` (os genes de cópia única presentes em todas as espécies, ponto de partida ideal para a árvore de espécies), e `Species_Tree/SpeciesTree_rooted.txt` (árvore de espécies já inferida automaticamente pelo próprio OrthoFinder).
- Em 2025/2026, uma versão substancialmente revisada foi publicada (bioRxiv, DOI: [10.1101/2025.07.15.664860](https://doi.org/10.1101/2025.07.15.664860); versão final em *Nature Methods*, DOI: [10.1038/s41592-026-03126-6](https://doi.org/10.1038/s41592-026-03126-6)), com ganho relatado de acurácia na delineação filogenética de ortogrupos e escalabilidade para milhares de genomas com menor uso de RAM/tempo de execução. Vale considerar se o projeto crescer para incluir dezenas de proteomas de referência de insetos.

**MAFFT** (Katoh & Standley, 2013, *Molecular Biology and Evolution* 30(4):772-780, DOI: [10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)) alinha as sequências de cada ortogrupo de cópia única. Para conjuntos de até algumas centenas de sequências, o modo mais preciso (`--auto` ou explicitamente `--linsi` para conjuntos pequenos) é viável computacionalmente:

```bash
mamba activate kerson-paper   # tem mafft
mafft --auto ortogrupo_OG0001234.fasta > ortogrupo_OG0001234.aln.fasta
```

**IQ-TREE** (versão atual: IQ-TREE 2; Minh et al., 2020, *Molecular Biology and Evolution* 37(5):1530-1534, DOI: [10.1093/molbev/msaa015](https://doi.org/10.1093/molbev/msaa015)) infere a árvore por máxima verossimilhança, com três componentes centrais: **ModelFinder** (seleção automática do melhor modelo de substituição por critério de informação, evita ter que "adivinhar" o modelo), busca de árvore eficiente, e **ultrafast bootstrap (UFBoot2)** para suporte de ramo, muito mais rápido que bootstrap tradicional sem inflar artificialmente os valores de suporte (correção documentada em Hoang et al., 2018, *Molecular Biology and Evolution*, DOI: [10.1093/molbev/msx281](https://doi.org/10.1093/molbev/msx281)). Uso típico, seguindo a mesma convenção já usada no Módulo Eugenio SPL deste laboratório:

```bash
mamba activate iqtree
iqtree2 -s alinhamento_concatenado.fasta -m MFP -bb 1000 -T AUTO -o especie_outgroup
```

- `-m MFP`: ModelFinder Plus, seleciona o melhor modelo antes de construir a árvore.
- `-bb 1000`: 1000 réplicas de ultrafast bootstrap.
- `-T AUTO`: detecção automática do número ótimo de threads (regra não-negociável do laboratório).
- `-o`: define outgroup explícito, essencial para enraizar a árvore de forma biologicamente interpretável (ex.: uma espécie de ordem mais basal fora do clado de interesse).
- IQ-TREE é a melhor escolha quando o número de sequências/táxons é gerenciável (dezenas a poucas centenas) e a acurácia do suporte de ramo importa. É o caso típico deste módulo: situar 1-2 espécies-praga entre parentes já sequenciados.

**FastTree** (FastTree 2; Price, Dehal & Arkin, 2010, *PLoS ONE* 5(3):e9490, DOI: [10.1371/journal.pone.0009490](https://doi.org/10.1371/journal.pone.0009490)) é uma alternativa de máxima verossimilhança aproximada, ordens de magnitude mais rápida que IQ-TREE/RAxML, ao custo de precisão de topologia e de suporte de ramo (usa um teste local tipo Shimodaira-Hasegawa em vez de bootstrap completo). Faz sentido quando:

- O número de táxons é muito grande (milhares) e um IQ-TREE completo com bootstrap se tornaria computacionalmente inviável no tempo disponível.
- A árvore é exploratória — por exemplo, uma primeira inspeção rápida de um ortogrupo grande e ruidoso antes de decidir se vale a pena refinar com IQ-TREE, ou uma árvore de triagem para depurar um alinhamento antes da análise "de verdade".
- **Não** é o passo final de uma análise que vai para publicação com interpretação evolutiva fina — para isso, IQ-TREE com ModelFinder + UFBoot é o padrão esperado por revisores.

Para o cenário intermediário — muitos táxons, mas ainda perto de máxima verossimilhança completa — existe o **VeryFastTree** (Piñeiro, Abuín & Pichel, 2020, *Bioinformatics* 36(17):4658-4659, DOI: [10.1093/bioinformatics/btaa582](https://doi.org/10.1093/bioinformatics/btaa582); versão 4.0 com melhorias adicionais de paralelização e computação em disco descrita em 2024, *GigaScience*, DOI: [10.1093/gigascience/giae055](https://doi.org/10.1093/gigascience/giae055)), reimplementação do FastTree 2 com os mesmos métodos e heurísticas (mesma interface de linha de comando, resultados equivalentes), mas paralelizada e vetorizada para rodar 3-8× mais rápido em datasets muito grandes (testado com até ~1 milhão de sequências). Não é uma ferramenta diferente em termos de método — é a mesma álgebra do FastTree 2, otimizada para hardware moderno. Não está instalada em nenhum ambiente do laboratório no momento; se um projeto futuro precisar processar milhares de táxons, é candidata natural via `mamba install -c bioconda veryfasttree`.

## Limitações de usar transcriptoma (não genoma) para ortologia

Esta é a ressalva mais importante do módulo, e é honesta: ferramentas como OrthoFinder foram desenhadas e validadas majoritariamente com **proteomas derivados de genomas anotados** (um conjunto limpo de um gene = uma proteína representativa). Usar como input um proteoma derivado de uma montagem de transcriptoma de novo (via TransDecoder) introduz problemas reais:

1. **Isoformas redundantes inflam artificialmente o proteoma.** Um assemblador de novo como Trinity tipicamente produz muito mais "genes" do que o real: no benchmark do UnigeneFinder (Xue et al., 2025, *Plant Direct*, DOI: [10.1002/pld3.70056](https://doi.org/10.1002/pld3.70056)), uma montagem Trinity de *Arabidopsis thaliana* gerou 117.928 transcritos contra apenas 19.366 genes reais no genoma, uma inflação de quase 6×; padrão semelhante apareceu nas outras três espécies testadas (172 mil a 400 mil transcritos brutos contra 22-25 mil genes esperados). Jogar todas essas isoformas no OrthoFinder sem filtragem produz ortogrupos artificialmente inflados, contagens de "duplicação gênica" espúrias e comparações de tamanho de família gênica que não são confiáveis — justo a análise que mais interessaria para investigar expansão de famílias de detoxificação. **A mitigação usada neste curso é a etapa de CD-HIT (Módulo 3, 95% identidade) seguida da seleção da isoforma mais longa/representativa por gene do TransDecoder antes de qualquer análise de ortologia.** Não pular essa etapa é o ponto mais crítico desta seção.
2. **Transcritos parciais/fragmentados.** Genes com baixa expressão no(s) tecido(s) sequenciado(s) frequentemente são reconstruídos de forma incompleta (5' ou 3' truncados, ou até só um éxon interno). Uma proteína predita parcial alinha mal, pode ser erroneamente classificada como não-ortóloga (falso negativo de ortologia) ou, pior, ancorar um alinhamento múltiplo de forma a distorcer a topologia da árvore resultante nesse ortogrupo.
3. **Cobertura de genes dependente do desenho experimental.** Um transcriptoma só contém os genes expressos nos tecidos/condições amostrados — diferente de um genoma, que captura (em teoria) 100% do conteúdo gênico independentemente de expressão. Genes de expressão muito restrita (certo estágio de desenvolvimento, resposta a estresse específico) podem estar simplesmente ausentes do transcriptoma de novo, gerando "perdas gênicas" espúrias em comparações de tamanho de família entre espécies com transcriptoma vs. espécies com genoma completo.
4. **Impacto documentado na acurácia filogenômica.** Este não é um problema novo: Yang & Smith (2014, *Molecular Biology and Evolution* 31(11):3081-3092, DOI: [10.1093/molbev/msu245](https://doi.org/10.1093/molbev/msu245)) já demonstravam que a inferência de ortologia a partir de transcriptomas (sem genoma) exige etapas adicionais de refinamento filogenético (poda de ramos espúrios em árvores de homólogos) para atingir acurácia e ocupação de matriz comparáveis às obtidas com genomas. Esse método está hoje incorporado em pipelines como PhyloTreePruner e inspirou a própria abordagem filogenética do OrthoFinder.

**Recomendação prática para os organismos-caso:** usar como proteoma de entrada no OrthoFinder o output pós-CD-HIT + TransDecoder + seleção de isoforma representativa (não o Trinity.fasta bruto), tratar contagens de expansão/contração de família gênica com cautela adicional (idealmente confirmando manualmente os genes de interesse via inspeção do alinhamento e da árvore de gene individual, não só da contagem no ortogrupo), e, sempre que possível, complementar a interpretação com evidência independente (ex.: presença/ausência de domínios Pfam esperados, RT-PCR confirmatório) antes de reportar uma expansão de família como resultado biológico.

## Prática usual vs. estado da arte (2024-2026)

| Aspecto | Prática usual/consolidada | Estado da arte 2024-2026 | Fonte (DOI) |
|---|---|---|---|
| Inferência de ortogrupos | OrthoFinder (2019), DIAMOND + árvores de gene, dezenas a centenas de proteomas | OrthoFinder revisado com delineação filogenética aprimorada de ortogrupos (+7% de acurácia relativa relatada) e escalabilidade para milhares de genomas com menor RAM/tempo | Nature Methods, 10.1038/s41592-026-03126-6 (preprint: bioRxiv 10.1101/2025.07.15.664860) |
| Redução de redundância pré-ortologia em transcriptoma sem genoma | CD-HIT (95% identidade) + seleção manual/heurística da isoforma mais longa | Pipelines automatizados dedicados que combinam múltiplos métodos de clustering para gerar transcrito primário, CDS e proteína de forma unificada (equivalente a um "genoma sintético" para fins de ortologia) | Plant Direct, UnigeneFinder, 10.1002/pld3.70056 |
| Filogenia em conjuntos pequenos/médios (dezenas a centenas de táxons) | IQ-TREE 2 com ModelFinder + UFBoot2 | Sem mudança de paradigma — IQ-TREE 2 continua padrão-ouro nessa faixa; ganhos recentes são de eficiência/paralelização, não de método novo | MBE, 10.1093/molbev/msaa015 |
| Filogenia em conjuntos muito grandes (milhares a ~1 milhão de táxons) | FastTree 2 (rápido, menos preciso) | VeryFastTree (reimplementação paralelizada/vetorizada do mesmo método do FastTree 2, 3-8× mais rápida; v4.0 testada com ~1 milhão de sequências em ~36h em um único servidor) | Bioinformatics 10.1093/bioinformatics/btaa582; GigaScience (v4.0) 10.1093/gigascience/giae055 |
| Ortologia a partir de transcriptoma sem genoma | Uso direto do proteoma TransDecoder, aceitando ruído residual | Reconhecimento formal, na literatura, de que refinamento filogenético (poda de homólogos espúrios) é necessário para acurácia comparável a genomas — base metodológica ainda válida e referenciada, não superada por método novo dedicado a "OrthoFinder para transcriptoma" | MBE, Yang & Smith 2014, 10.1093/molbev/msu245 |

## Recomendação para quem quer aprofundar

Para o escopo típico deste curso (situar 1-2 espécies-praga sem genoma entre um punhado de espécies de referência já sequenciadas, dezenas de proteomas no máximo), a combinação **OrthoFinder (2019 ou a versão 2026 revisada, se disponível no ambiente) → MAFFT → IQ-TREE 2 com ModelFinder + UFBoot2** continua sendo a escolha correta e é a que reflete o estado da arte atual, sem exigir infraestrutura especial. FastTree/VeryFastTree só entram em cena se o projeto crescer para centenas ou milhares de proteomas de referência (cenário improvável para os organismos-caso deste curso, mas relevante, por exemplo, para uma futura análise pan-Lepidoptera ou pan-Hemiptera de uma família gênica específica entre dezenas de espécies publicadas). O ponto que mais vale revisitar antes de rodar qualquer análise de ortologia neste laboratório é a etapa de preparação do proteoma de entrada (isoforma única por gene, pós-CD-HIT), não a escolha entre as ferramentas de árvore em si — é aí que mora o maior risco de viés silencioso nos resultados.

## Aplicação aos organismos-caso

**Anticarsia gemmatalis** (Lepidoptera, lagarta-da-soja) pode ser posicionada filogeneticamente ao lado de proteomas de referência já publicados e genomicamente anotados de outros lepidópteros — candidatos naturais incluem *Bombyx mori* (bicho-da-seda, genoma de referência histórico do grupo), *Spodoptera frugiperda* (já presente neste próprio laboratório via o pipeline RNA-Seq-not-model, o que facilita reuso do proteoma) e *Helicoverpa armigera* (praga agrícola global bem estudada). Rodar OrthoFinder com esses proteomas + o de *A. gemmatalis* permite (a) confirmar a topologia esperada (posição relativa dentro de Noctuoidea), como controle de qualidade da própria análise, e (b) identificar, dentro dos ortogrupos de citocromos P450 e esterases, se *A. gemmatalis* apresenta contagem de genes destoante das demais espécies — um primeiro sinal quantitativo de expansão de família associável a resistência a inseticida, para investigação mais aprofundada.

**Mahanarva spectabilis** (Hemiptera, cigarrinha-das-pastagens) segue a mesma lógica, mas com referências de Hemiptera — candidatos incluem *Nilaparvata lugens* (praga do arroz, genoma publicado, mesma ordem), *Rhodnius prolixus* (vetor de Chagas, bem anotado, ainda que filogeneticamente mais distante dentro de Hemiptera) e, se disponível, *Nezara viridula* ou outro pentatomídeo/cercopídeo próximo. Como o transcriptoma de *M. spectabilis* já tem um componente de investigação de endossimbiontes e toxinas de glândula salivar (ver projeto trinity_maharnava), a análise de ortologia aqui tem valor adicional: genes de secreção salivar sem ortólogo claro em nenhuma espécie de referência (ortogrupos espécie-específicos ou "genes órfãos") são candidatos interessantes a toxinas ou efetores específicos da linhagem, potencialmente relacionados à indução de sintomas na planta hospedeira (o "amarelão") — uma leitura direta do resultado de ortologia que vai além de apenas "onde essa espécie se encaixa na árvore".

Em ambos os casos, o output mais imediatamente acionável do fluxo OrthoFinder → IQ-TREE não é a árvore de espécies em si (que tende a confirmar relações já conhecidas da taxonomia), mas a **tabela de contagem de genes por ortogrupo** (`Orthogroups.GeneCount.tsv`) cruzada com a anotação funcional do Módulo 4 — é ali que aparecem, de forma quantitativa e reprodutível, os candidatos a expansão de família ou a genes órfãos que merecem investigação funcional mais aprofundada.

# Módulo 4 — Predição de ORFs e Anotação Funcional

*(60 min — parte teórica do minicurso "Montagem de Transcriptoma de Eucariotos Não-Modelo via RNA-Seq")*

## Conceitos-chave

Depois que a montagem de novo produz um conjunto de transcritos (ex.: 42.372 transcritos Trinity em *Anticarsia gemmatalis*), o transcrito bruto ainda não diz nada sobre função biológica: é só uma sequência de nucleotídeos, muitas vezes carregando UTRs, isoformas redundantes e artefatos de montagem. O Módulo 4 cobre a ponte entre "transcrito montado" e "gene com função conhecida", que passa por três etapas conceituais distintas:

1. **Predição de ORFs (Open Reading Frames):** identificar, dentro de cada transcrito, qual é a região codificante mais provável e qual quadro de leitura ela usa. Isso é necessário porque um transcrito Trinity não vem anotado com CDS — ele é a sequência montada, e pode conter a proteína completa, um fragmento 5' ou 3' incompleto, ou ser puramente não-codificante (lncRNA, contaminação, UTR isolado por erro de montagem).
2. **Busca por homologia (annotation by similarity):** comparar as proteínas preditas contra bancos de referência (Swiss-Prot, TrEMBL, NCBI nr) para inferir função por similaridade de sequência. É a estratégia mais indicada para organismos não-modelo, porque dispensa um genoma de referência anotado do próprio organismo.
3. **Anotação funcional estruturada:** converter os hits de homologia (e buscas complementares por domínios/família) em termos padronizados e comparáveis entre estudos — GO (Gene Ontology), KEGG (KO/pathways), Pfam/InterPro, COG —, que alimentam as análises de enriquecimento do módulo seguinte (clusterProfiler).

Para organismos não-modelo, há um ponto que este módulo repete de propósito: anotar não é rodar BLAST e ficar com o melhor hit. Amostras de insetos costumam carregar material genético de simbiontes, patógenos e contaminantes ambientais junto com o RNA do hospedeiro, e sem uma etapa explícita de filtragem, uma fração nada trivial do transcriptoma "anotado" pode pertencer a bactérias, fungos ou até ao operador que preparou a biblioteca (contaminação humana). Isso vale sobretudo para os dois organismos-caso do curso: *A. gemmatalis* (Lepidoptera, herbívoro foliar, pode carregar vírus e microbiota intestinal) e *Mahanarva spectabilis* (Hemiptera sugador de floema, grupo taxonômico classicamente associado a endossimbiontes obrigatórios como Sulcia e Sodalis-like).

## Ferramentas e fluxo

### 1. TransDecoder — predição de ORFs

TransDecoder (Haas et al., parte do ecossistema Trinity/CO-genome) identifica candidatos a região codificante dentro de transcritos montados de novo. A versão estável usada em produção é a **5.7.1** (bioconda, build de julho/2023); já existe uma **6.0.0** mais recente no bioconda, e o desenvolvedor original avisou que o projeto histórico não recebe mais suporte ativo desde março/2024, com um sucessor em desenvolvimento chamado **TransDecoder2 (TD2)**. Para o minicurso, a recomendação é fixar a 5.7.1 (mamba: `transdecoder=5.7.1`), a mesma versão usada nos pipelines já validados do laboratório (`RNA-Seq-not-model`, `trinity_maharnava`).

Fluxo em duas etapas:

```bash
# Etapa 1 — extrai todas as ORFs candidatas longas (default: >=100 aa)
TransDecoder.LongOrfs -t transcritos.fasta

# Etapa 2 (opcional, mas recomendada) — evidência de homologia para reduzir falsos negativos
# ORFs curtas e reais (ex.: peptídeos secretados pequenos) são frequentemente descartadas
# pelo critério puramente estatístico do TransDecoder; buscas de homologia "resgatam" essas ORFs.
diamond blastp -d swissprot.dmnd -q transdecoder_dir/longest_orfs.pep \
  --max-target-seqs 1 --outfmt 6 -e 1e-5 > blastp.outfmt6

hmmsearch --cpu 16 --domtblout pfam.domtblout Pfam-A.hmm \
  transdecoder_dir/longest_orfs.pep

# Etapa 3 — predição final, informada pelas evidências acima
TransDecoder.Predict -t transcritos.fasta \
  --retain_pfam_hits pfam.domtblout \
  --retain_blastp_hits blastp.outfmt6 \
  --single_best_only
```

`--single_best_only` mantém apenas uma ORF por transcrito (a de maior score) — prática usual quando o objetivo é um catálogo de genes/proteínas não-redundante, como no pipeline de *M. spectabilis*, que chegou a 12.445 proteínas preditas a partir da montagem Trinity.

### 2. DIAMOND — busca de homologia rápida

DIAMOND (Buchfink, Xie & Huson, *Nature Methods* 12, 59–60, 2015, DOI: `10.1038/nmeth.3176`) é hoje o substituto de facto do BLAST clássico em escala de transcriptoma. O paper original relata ganhos de **até 20.000× sobre BLASTX** em leituras curtas, com sensibilidade comparável — essa é a cifra correta da literatura primária (não "~1000×"; o ganho de ~1000× aparece em benchmarks de blastp proteína-proteína completo, enquanto os 20.000× são especificamente para blastx de reads curtas contra bancos grandes). Em 2021 saiu o DIAMOND v2 (Buchfink, Reuter & Drost, *Nature Methods* 18, 366–368, DOI: `10.1038/s41592-021-01101-x`), que trouxe os modos de sensibilidade escalonados — `--fast`, `--mid-sensitive`, `--sensitive`, `--more-sensitive`, `--very-sensitive`, `--ultra-sensitive` — permitindo alinhamentos proteína-proteína em escala "árvore da vida" sem abrir mão de sensibilidade próxima à do BLASTP tradicional.

Parâmetros-chave para anotação de transcriptoma não-modelo:

```bash
diamond blastp \
  -q proteinas.transdecoder.pep \
  -d uniprot_sprot.dmnd \
  --more-sensitive \
  --evalue 1e-5 \
  --query-cover 50 \
  --subject-cover 50 \
  --max-target-seqs 1 \
  --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
  --threads 16 \
  -o diamond_sprot.outfmt6
```

- `--evalue 1e-5`: limiar padrão para reduzir hits espúrios em transcriptomas de novo, onde erros de montagem geram ORFs artificiais.
- `--query-cover` / `--subject-cover` ≥ 50%: evita anotar por hits parciais (domínio compartilhado, mas proteína diferente) — importa em famílias multigênicas como LRR-RLPs (ver projeto Kerson-paper) e serino-proteases (tripsinas).
- `--sensitive`/`--ultra-sensitive`: vale usar quando o banco alvo é filogeneticamente distante do organismo (ex.: comparar Hemiptera contra Swiss-Prot, dominado por vertebrados/organismos modelo); custo computacional maior, mas ainda ordens de magnitude abaixo do BLAST clássico.
- Bancos recomendados, em ordem de uso: **Swiss-Prot** (curado, alta confiança, primeira passada) → **TrEMBL/UniProt** ou **NCBI nr/RefSeq** (cobertura ampla, segunda passada para o que não anotou contra Swiss-Prot).

### 3. EnTAP — pipeline de anotação dedicado a eucariotos não-modelo

EnTAP (Eukaryotic Non-Model Transcriptome Annotation Pipeline; Hart, Wagner & Wegrzyn, *Molecular Ecology Resources* 20:591–604, 2020, DOI: `10.1111/1755-0998.13106`) é hoje a referência mais citada para exatamente o problema deste módulo: acelerar e tornar mais confiável a anotação de transcriptomas de eucariotos sem genoma de referência. A versão certa é a **2.2.0** — a documentação oficial em `entap.readthedocs.io` confirma; vale conferir o changelog antes da aula, já que o desenvolvimento ativo migrou do GitHub (`harta55/EnTAP`) para o GitLab (`gitlab.com/PlantGenomicsLab/EnTAP`). A release já empacota versões testadas de dependências (RSEM 1.3.3, TransDecoder 5.7.1, DIAMOND 2.1.8), o que ajuda na reprodutibilidade.

O que o EnTAP faz além de "rodar DIAMOND": (1) filtragem por expressão (remove transcritos com suporte de expressão insuficiente, reduzindo transcritos fragmentados/artefatos de montagem); (2) seleção de frame consistente com TransDecoder; (3) busca de similaridade contra múltiplos bancos simultaneamente, com deduplicação e escolha do "melhor hit informativo" (evita que um hit "hypothetical protein" de bitscore alto vença um hit funcionalmente informativo de bitscore ligeiramente menor); (4) filtro filogenético/contaminante configurável por lineage — o ponto que importa no cenário de inseto com simbiontes; (5) atribuição de domínios (Pfam/InterPro via EggNOG), GO e KEGG.

### 4. Trinotate — alternativa clássica integrada

Trinotate é a suíte de anotação historicamente acoplada ao ecossistema Trinity: BLAST+/Swiss-Prot para homologia, HMMER/Pfam para domínios, SignalP/TMHMM para peptídeo sinal e domínios transmembrana, e integração com bancos eggNOG/GO/KEGG, tudo compilado em um banco SQLite consultável e exportável como "annotation report" tabular. O repositório oficial confirma que Trinotate não está mais em desenvolvimento ativo desde março/2024 — os mantenedores recomendam continuar usando enquanto atender às necessidades do projeto, com fork livre para quem precisar estender. A versão mais recente amplamente disponível é a **4.0.2**. Para o minicurso, Trinotate continua sendo útil por integrar várias análises num único relatório, mas não deve ser apresentado como "estado da arte" — é a opção clássica/legada, enquanto EnTAP é a opção ativa recomendada para não-modelo.

### 5. eggNOG-mapper — ortologia e GO/KEGG/COG

eggNOG-mapper atribui GO terms, KEGG Orthology (KOs) e categorias COG por transferência de ortologia (não por melhor-hit simples), o que tende a ser mais preciso que anotação por similaridade pura quando o objetivo é uma classificação funcional comparável entre espécies. A série 2.1.x segue ativa: a **2.1.13** foi lançada em 18/06/2025 (correção de parsing de GFF do Prokka) e a **2.1.14** em 04/05/2025, corrigindo a URL do script de download do banco — ambas disponíveis via bioconda. Não há, até a data da pesquisa, uma versão "2.2" ou major nova publicada; a citação "v2" do artigo de referência (Cantalapiedra et al. 2021, *Molecular Biology and Evolution*, "eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale") continua sendo a descrição correta da arquitetura atual.

```bash
emapper.py -i proteinas.transdecoder.pep \
  --itype proteins \
  -m diamond \
  --cpu 16 \
  --output eggnog_annotation \
  --output_dir resultados_eggnog/
```

### 6. InterProScan — domínios, famílias e GO via evidência combinada

InterProScan integra cerca de 15 bancos de assinaturas (Pfam, PROSITE, PRINTS, SMART, CDD, PANTHER etc.) num único relatório com mapeamento para GO e pathways. Na série 5, a versão mais recente publicada é a **5.78-109.0** (o segundo número indica a versão do banco InterPro sincronizada, não um "minor" convencional; a 5.77-108.0 foi lançada em 29/01/2026). O artigo de referência histórico é Jones et al., *Bioinformatics* 30(9):1236–1240, 2014, DOI: `10.1093/bioinformatics/btu031`. Vale registrar uma novidade que não constava no material original: em 2026 a EBI publicou "InterProScan 6: a modern large-scale protein function annotation pipeline" (Blum, Hobbs, Florentino & Bateman, *Bioinformatics Advances* 6(1):vbag141, DOI: `10.1093/bioadv/vbag141`), anunciando uma reimplementação completa como pipeline Nextflow, com foco em escalabilidade, portabilidade (local/HPC/cloud) e reprodutibilidade — o que é relevante para este público, que já usa Nextflow no pipeline MD-GROMACS. Vale citar como novidade a acompanhar, mas ainda não é a opção default recomendada em produção enquanto a comunidade migra e valida (a série 5 é o padrão estabelecido no momento do curso).

Paralelização prática da série 5 (relevante no servidor de 32 cores):

```bash
interproscan.sh -i proteinas.transdecoder.pep \
  -f tsv,gff3 \
  -goterms -pa \
  --cpu 32 \
  -T /tmp/ips_tmp
```

`-pa` habilita mapeamento para pathways (KEGG/Reactome/MetaCyc via associação com GO/InterPro); `-goterms` adiciona GO. Em transcriptomas com mais de 10 mil proteínas, rodar em lotes (ex.: 2.000 sequências por chamada) evita estouro de memória e permite checkpoint/retomada.

## Prática usual vs. estado da arte (2024-2026)

| Aspecto | Prática usual/consolidada | Estado da arte 2024-2026 | Fonte (DOI) |
|---|---|---|---|
| Predição de ORFs | TransDecoder 5.x (LongOrfs+Predict, com evidência Pfam/BLAST) | TransDecoder2 (TD2) anunciado como sucessor; TransDecoder 6.0.0 já disponível, mas sem paper de validação amplo ainda | github.com/TransDecoder (sem DOI formal para TD2 até a pesquisa) |
| Busca de homologia | DIAMOND blastp/blastx contra Swiss-Prot/nr (v2.1.x) | DIAMOND v2 com modos de sensibilidade escalonados para buscas "árvore da vida" | 10.1038/s41592-021-01101-x |
| Pipeline de anotação não-modelo | EnTAP 2.2.0 (filtro de expressão+frame+contaminante) | Ainda EnTAP é a referência consolidada citada como melhor prática para filtragem de contaminantes em não-modelo; não foi identificado substituto amplamente adotado equivalente no período 2024-2026 | 10.1111/1755-0998.13106 |
| Domínios/famílias | InterProScan série 5 (5.78-109.0) | InterProScan 6 — reimplementação Nextflow, escalável, publicada em 2026; ainda em adoção inicial pela comunidade | 10.1093/bioadv/vbag141 |
| Ortologia/GO/KEGG | eggNOG-mapper 2.1.x (versão atual 2.1.13/2.1.14, 2025) | Mesma arquitetura v2 (Cantalapiedra et al. 2021); sem major release nova identificada até a pesquisa | eggNOG-mapper v2, *Mol Biol Evol* 2021 |
| Anotação funcional "por similaridade" | BLAST/DIAMOND melhor-hit + domínios | Modelos de linguagem de proteínas (protein language models, PLMs) — ESM-2 como encoder, ferramentas como ProtNLM, ProteInfer, ProtNote, NetGO 3.0 — usados para prever função sem depender de homólogo próximo no banco; ganhando tração acadêmica mas ainda pouco adotados em pipelines de rotina de transcriptômica não-modelo | ProtNote: 10.1093/bioinformatics/btaf170; NetGO 3.0: bioRxiv 10.1101/2022.12.05.519073 |

**Leitura do quadro para o curso:** a pesquisa não encontrou evidência de que EnTAP tenha sido superado como prática recomendada para filtragem de contaminantes/anotação integrada em transcriptomas de novo de insetos — ele continua sendo a citação padrão em publicações recentes da área. O avanço mais concreto do período 2024-2026 é o anúncio do InterProScan 6, que merece menção como tendência a observar, mas não substitui a série 5 no fluxo prático deste minicurso. Já os modelos de linguagem de proteínas (PLMs) representam a fronteira de pesquisa em predição de função quando não há homólogo caracterizado — situação comum em Hemiptera/Lepidoptera não-modelo —, mas ainda são ferramentas de pesquisa/benchmark (CAFA), não pipelines turnkey para o tipo de anotação em lote que EnTAP/Trinotate oferecem. Vale mencionar como direção futura, sem recomendar substituição do pipeline consolidado agora.

## Filtragem de contaminantes e endossimbiontes

Transcriptomas de novo de insetos capturam RNA de tudo que estava biologicamente ativo na amostra, não só do hospedeiro. As três categorias de contaminação mais comuns:

1. **Contaminação técnica:** DNA/RNA humano (manuseio de laboratório), leveduras/fungos de bancada, vetores de clonagem residuais.
2. **Micróbios ambientais/patógenos oportunistas:** bactérias de superfície corporal ou intestino, especialmente em amostras de campo (não axênicas) como as coletadas para *M. spectabilis*.
3. **Endossimbiontes obrigatórios/facultativos — o caso mais importante para a aula:** insetos sugadores de seiva (Hemiptera, como cigarrinhas-das-pastagens do gênero *Mahanarva*) dependem quase universalmente de simbiontes bacterianos para suplementar uma dieta pobre em nutrientes (aminoácidos essenciais, vitaminas B). Os gêneros clássicos são:
   - ***Buchnera aphidicola*** — simbionte primário obrigatório de afídeos (Hemiptera: Aphididae), o exemplo mais estudado de coevolução hospedeiro-simbionte com genoma extremamente reduzido.
   - ***Sulcia muelleri*** e simbiontes **Sodalis-like** — comuns em Auchenorrhyncha (cigarrinhas, cigarras), o clado ao qual pertence *Mahanarva*; no projeto real de *M. spectabilis*, a anotação de fato confirmou hits contra proteínas Sulcia/Sodalis-like (58 hits), mostrando que o transcriptoma de glândula salivar capturou não só o hospedeiro, mas também sua microbiota simbiótica.
   - ***Wolbachia*** — endossimbionte facultativo intracelular extremamente disseminado em artrópodes (estima-se presença em >40% das espécies de insetos), manipulador reprodutivo, frequentemente encontrado "de graça" em bibliotecas de RNA-Seq não direcionadas a ele.

Isso não é necessariamente um problema — pode até ser um achado biológico relevante, como no caso de *M. spectabilis* — mas precisa ser identificado e tratado de forma explícita, nunca deixado misturado ao catálogo de genes do hospedeiro. Estratégia recomendada:

- No EnTAP, configurar o filtro de contaminantes por lineage/táxon (ex.: excluir/marcar hits cujo melhor lineage LCA caia em Bacteria quando o organismo-alvo é Insecta), mantendo esses transcritos em uma tabela separada em vez de descartá-los — eles podem virar uma seção própria de "achados de simbiose" no artigo.
- Cruzar independentemente contra bancos de genomas de simbiontes conhecidos do clado (Buchnera, Sulcia, Sodalis, Wolbachia) via DIAMOND blastx com `--taxonlist`/checagem manual do táxon do melhor hit, quando o filtro automático do EnTAP não for granular o suficiente.
- Para Lepidoptera como *A. gemmatalis*, o risco em endossimbiontes obrigatórios costuma ser menor (dependem menos de simbiose nutricional que Hemiptera sugadores de floema/seiva), mas vírus continuam relevantes — o baculovírus, aliás, é usado como bioinseticida contra essa espécie — assim como bactérias intestinais.

## Recomendação para o curso

Pipeline recomendado para os dois organismos-caso, em ordem:

1. `TransDecoder.LongOrfs` → `TransDecoder.Predict --single_best_only`, com evidência de `diamond blastp` contra Swiss-Prot e `hmmsearch` contra Pfam-A para resgatar ORFs curtas reais.
2. `diamond blastp --more-sensitive` das proteínas preditas contra Swiss-Prot (primeira passada, alta confiança) e, para o que não anotou, contra TrEMBL/nr (segunda passada, cobertura).
3. `eggNOG-mapper` (modo `diamond`) para GO/KEGG/COG via transferência de ortologia — complementa (não substitui) a anotação por melhor-hit.
4. `InterProScan` (série 5, `-goterms -pa`) para domínios Pfam/InterPro e pathways, rodado em lotes no servidor com `--cpu 32`.
5. **EnTAP como camada de integração e QC final**, sobretudo para o filtro de contaminantes/lineage — é a etapa que transforma os outputs brutos das ferramentas acima em uma anotação "de publicação", já filtrada e com melhor-hit informativo escolhido de forma consistente.
6. Trinotate pode ser oferecido como alternativa mais simples de configurar para quem quer um único banco SQLite integrado, mas o curso deve deixar claro que é a opção legada (sem desenvolvimento ativo desde 2024), recomendada apenas quando EnTAP não estiver disponível/configurável a tempo.

Justificativa: nenhuma ferramenta isolada cobre as três necessidades de um transcriptoma de inseto não-modelo ao mesmo tempo — velocidade (DIAMOND), completude ortológica (eggNOG-mapper/InterProScan) e controle de qualidade/contaminação (EnTAP). Essa combinação repete o desenho já usado com sucesso no pipeline de *M. spectabilis* do laboratório (DIAMOND/eggNOG-mapper/Pfam → GO/KEGG), com o filtro EnTAP como boa prática adicional.

## Erro comum — estudo de caso real

No pipeline de *M. spectabilis* (`trinity_maharnava`), a etapa de merge da tabela de anotação do eggNOG-mapper com a tabela de proteínas/GO do projeto tinha um bug de índice de coluna: o script assumia uma posição fixa de coluna no arquivo de output do eggNOG-mapper (ex.: "a coluna 9 sempre é a descrição") para fazer o merge programático via pandas, unindo por posição em vez de por nome de coluna/chave. O problema é que o output do eggNOG-mapper pode variar em número e ordem de colunas dependendo da versão da ferramenta, das flags usadas (`--tax_scope`, `--target_orthologs`, presença ou não de anotação de PFAMs no mesmo arquivo) e de linhas de cabeçalho/comentário (`##`) que o parser não tratou corretamente. O resultado foi um merge silenciosamente incorreto: descrições funcionais e termos GO foram atribuídos ao transcrito errado, sem qualquer erro ou warning visível — o script rodou "com sucesso" e produziu uma tabela de anotação inteira, só que deslocada.

**Lição generalizável de QC para todo o módulo:**

- **Nunca faça merge de tabelas de anotação por posição/índice de coluna.** Carregue os dados usando os nomes de coluna explícitos do cabeçalho (`pandas.read_csv(..., sep="\t", comment="##")` seguido de merge por `on="query"` ou chave equivalente nomeada), nunca por `df.iloc[:, 8]` ou equivalente.
- **Valide o cabeçalho a cada execução**, pois versões diferentes da ferramenta (eggNOG-mapper 2.1.13 vs. 2.1.14, por exemplo) podem alterar sutilmente o número/nome de colunas entre releases.
- **Depois de qualquer merge de anotação, faça uma checagem de sanidade manual**: escolha 5-10 transcritos ao acaso, confira se a descrição funcional atribuída é plausível para o hit de homologia independente (ex.: bate com o gene esperado ou pelo menos com a categoria funcional esperada), antes de prosseguir para análises de enriquecimento (GO/KEGG) que dependem inteiramente da correção desse merge.
- Esse tipo de bug é particularmente perigoso porque não gera erro de execução — o pipeline "funciona", os arquivos são gerados, e o problema só aparece quando alguém nota uma incoerência biológica óbvia, como aconteceu na auditoria que motivou a correção no projeto. Trate merges de anotação com o mesmo rigor de um teste unitário: escreva uma checagem automática que confirme que o número de linhas final bate com o esperado e que uma amostra de chaves (`transcript_id`) realmente corresponde entre as tabelas de origem.

## Aplicação aos organismos-caso

**Anticarsia gemmatalis (Lepidoptera, 42.372 transcritos Trinity):** aplicar TransDecoder com evidência de homologia/Pfam para reduzir o número de ORFs espúrias antes da anotação — transcriptomas de novo tendem a superestimar o número de "genes" por causa de isoformas e transcritos fragmentados. Priorizar DIAMOND contra Swiss-Prot de insetos/artrópodes como primeira passada faz sentido aqui, já que essa espécie tem parentes bem anotados (outros Lepidoptera, ex. *Bombyx mori*, *Helicoverpa* spp.) em bancos públicos, o que deve gerar boa taxa de anotação por similaridade direta.

**Mahanarva spectabilis (Hemiptera, 12.445 proteínas preditas via TransDecoder):** aqui o filtro de contaminante/simbionte é obrigatório, não opcional — dado o clado (Auchenorrhyncha, sugador de floema/xilema), era esperado encontrar sinal de Sulcia/Sodalis-like na anotação, e foi o que de fato ocorreu (58 hits confirmados). Ao reproduzir esse pipeline em sala, vale mostrar explicitamente a tabela final com uma coluna de "lineage/reino do melhor hit", para que os alunos vejam a separação entre proteínas do hospedeiro e proteínas de simbionte lado a lado — e usar esse mesmo dataset para demonstrar ao vivo a lição do bug de merge de coluna descrita acima, já que foi o cenário real onde o problema ocorreu e foi corrigido.

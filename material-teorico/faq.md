# FAQ — Perguntas Frequentes

> ~50 perguntas organizadas por módulo, sintetizadas a partir dos pontos de dúvida, erros comuns e ressalvas identificados na pesquisa de cada módulo. Uso sugerido: distribuir aos alunos como material de consulta pós-aula.

## Módulo 0 — Fundamentos e desenho experimental

1. **Minha espécie não tem genoma de referência, mas um parente próximo tem. Ainda preciso montar de novo?** Depende do objetivo: se é só quantificar expressão em genes conservados, mapeamento cross-species pode bastar; se o objetivo inclui descoberta de genes/isoformas específicos da espécie, monte de novo.
2. **Quantas réplicas biológicas eu preciso?** Para a montagem em si, réplicas não são estritamente necessárias (o foco é maximizar diversidade de transcritos); para expressão diferencial subsequente, o mínimo defensável é n=3 por condição.
3. **Quantos reads por amostra eu preciso sequenciar?** A faixa consolidada é 20-40 milhões de reads pareados por amostra (Francis et al. 2013), com ganho marginal além de ~60M; não há revisão numérica publicada em 2024-2026 que mude isso.
4. **Vale a pena usar long-read (PacBio/ONT) em vez de Illumina?** Long-read resolve isoformas de forma mais direta, mas custa mais por Gb e tem taxa de erro maior (embora HiFi/basecalling recentes tenham reduzido essa diferença). Short-read puro continua cientificamente válido.
5. **Devo combinar amostras de tecidos diferentes numa única montagem?** Sim, para maximizar diversidade de transcritos capturados — mas separe as réplicas biológicas na etapa de quantificação/DE, não na montagem.
6. **Biblioteca stranded vale o custo extra?** Vale, na maioria dos casos — sobretudo sem anotação prévia, pois evita fundir por engano genes vizinhos em fitas opostas.
7. **Onde encontro dados públicos de RNA-Seq de insetos não-modelo?** NCBI SRA e ENA são os principais; bancos especializados como i5k Workspace/InsectBase agregam metadados específicos de insetos.
8. **"Montagem de novo" e "genoma de novo" são a mesma coisa?** Não — montagem de novo de transcriptoma reconstrói apenas os transcritos expressos a partir de RNA-Seq, sem tentar reconstruir o genoma inteiro.

## Módulo 1 — Qualidade e pré-processamento

9. **Por que o QC é mais rigoroso em montagem de novo do que em RNA-Seq com referência?** Porque não há aligner "absorvendo" ruído via soft-clipping — erro residual vira nó extra no grafo de De Bruijn e pode gerar transcrito quimérico/fragmentado.
10. **Duplicação alta no FastQC é sempre problema?** Não em RNA-Seq — genes muito expressos geram muitas cópias idênticas de fragmentos; isso é esperado, diferente de WGS.
11. **Devo sempre rodar SortMeRNA?** Só se o QC (bruto ou pós-trimming) mostrar evidência real de rRNA residual — rodar "por precaução" gasta RAM/tempo sem necessidade.
12. **Normalização in silico do Trinity é sempre recomendada?** Não — só quando o dataset é muito grande (tipicamente >200M reads) e vira gargalo de RAM; pode distorcer isoformas raras se aplicada sem necessidade.
13. **Qual a ordem correta de QC/trimming?** FastQC/MultiQC → rcorrector → remoção de reads unfixable → fastp → SortMeRNA (se necessário) → FastQC/MultiQC de novo → normalização (se necessário).
14. **Por que rodar rcorrector antes do fastp, e não depois?** rcorrector opera melhor em reads ainda não cortadas, e corrige a base errada em vez de simplesmente descartá-la, preservando mais comprimento útil.
15. **É normal ver pico de %GC fora do esperado no FastQC?** Pode ser sinal de contaminação (rRNA, DNA genômico) ou, em insetos com endossimbiontes, biologia real do sistema hospedeiro-simbionte — investigar antes de descartar.
16. **`--trim_poly_g` deve ser sempre ativado no fastp?** Não — é específico de plataformas de dois canais (NovaSeq/NextSeq); confirmar a plataforma antes de ativar.
17. **Esqueci de reavaliar o QC pós-trimming, isso é grave?** É o erro mais comum e mais barato de evitar — sem essa checagem não há como confirmar se o trimming resolveu o problema identificado no QC bruto.

## Módulo 2 — Montagem de novo

18. **Trinity ou rnaSPAdes — qual escolher?** Nenhum é estritamente superior: rnaSPAdes tende a vencer em contiguidade/uso de RAM, Trinity tende a vencer em sensibilidade a isoformas raras. A escolha deve ser empírica (BUSCO/N50) por dataset.
19. **O que significa o "gene" (`g#`) no ID de um transcrito Trinity?** É um agrupamento heurístico de isoformas por similaridade de grafo, não um gene biológico confirmado — pode fundir parálogos ou fragmentar genes de baixa expressão.
20. **Vale a pena usar um pipeline multi-montador (TransPi/ORP) em vez de um único montador?** Sim, é a recomendação de estado da arte — benchmarks mostram maior completude BUSCO e menor duplicação do que qualquer montador isolado, ao custo de mais RAM/CPU.
21. **Quanta RAM preciso para rodar o Trinity?** Regra de bolso ainda válida: ~1 GB por milhão de pares de reads, mas é um piso — heterozigosidade e complexidade do transcriptoma podem exigir bem mais.
22. **Trinity ainda é mantido ativamente?** O financiamento dedicado do ITCR terminou em 2024; o software continua funcional e open-source, mas com ritmo de manutenção reduzido.
23. **Tenho dados PacBio/ONT sem genoma de referência — o que uso?** RNA-Bloom2 é hoje a opção mais bem documentada e validada para montagem long-read reference-free.
24. **É problema usar Trinity puro em vez de um pipeline multi-montador?** Não invalida os resultados (desde que passem por validação funcional a jusante), mas é uma limitação metodológica que vale declarar em qualquer manuscrito.
25. **O que é EvidentialGene e quando uso?** É o "oráculo" que consolida transcritos de múltiplos montadores/k-mers num conjunto primário não-redundante — usado dentro de pipelines como TransPi.

## Módulo 3 — Redução de redundância e avaliação

26. **N50 alto significa montagem boa?** Não necessariamente — N50 mede só tamanho, não correção nem completude; um assembler "preguiçoso" pode inflar N50 concatenando fragmentos não relacionados.
27. **Qual linhagem BUSCO devo usar?** A mais específica disponível para o clado (ex.: `hemiptera_odb12` em vez de `insecta_odb12` ou `arthropoda_odb12`, quando aplicável).
28. **BUSCO odb10 ou odb12?** odb12 é a geração atual (2024-2025); odb10 continua funcional para comparação com literatura publicada, mas não recebe mais expansão de cobertura.
29. **BUSCO "D" (Duplicated) alto é sempre problema?** Pode ser duplicação gênica real em alguns clados, mas junto com C alto é sinal clássico de redução de redundância insuficiente — vale revisar CD-HIT/EvidentialGene.
30. **BUSCO "M" (Missing) alto significa montagem ruim?** Nem sempre — às vezes reflete o tecido amostrado, já que um tecido específico simplesmente não expressa genes de outros sistemas do organismo. Montagem incompleta é só uma das explicações possíveis.
31. **CD-HIT-EST ou EvidentialGene para reduzir redundância?** São complementares: CD-HIT é rápido e cego a contexto biológico (clustering por identidade pura); EvidentialGene é mais sofisticado (usa estrutura de ORF/CDS) mas mais pesado computacionalmente.
32. **TransRate ainda é recomendado?** Ainda aparece em papers recentes, mas o desenvolvimento parou há anos — tratar como métrica complementar histórica, não padrão-ouro atual.
33. **O que é ExN50 e por que é melhor que N50?** Pondera N50 pela contribuição real de expressão de cada transcrito — penaliza menos o "lixo" de baixa cobertura que acumula em toda montagem de novo.
34. **Que % de back-mapping é considerada adequada?** Geralmente ≥80-90%; abaixo disso vale investigar contaminação ou perda real de transcriptoma.
35. **BUSCO alto mas back-mapping baixo — o que fazer?** Investigar contaminação: BUSCOs são genes conservados pouco afetados por contaminação de baixo nível, então essa combinação é consistente com material externo na biblioteca.

## Módulo 4 — Predição de ORFs e anotação funcional

36. **Por que rodar DIAMOND/Pfam antes do `TransDecoder.Predict`?** Para resgatar ORFs curtas reais (ex.: peptídeos secretados pequenos) que o critério puramente estatístico do TransDecoder tende a descartar.
37. **Qual banco usar primeiro na busca de homologia — Swiss-Prot ou nr?** Swiss-Prot primeiro (curado, alta confiança); o que não anotar, buscar contra TrEMBL/nr (cobertura ampla).
38. **EnTAP ou Trinotate?** EnTAP é a opção ativa recomendada para não-modelo (filtro de contaminantes/expressão/frame integrado); Trinotate é a opção clássica/legada, sem desenvolvimento ativo desde 2024.
39. **Como lido com contaminação por endossimbiontes em insetos sugadores de seiva?** Configure filtro de contaminante por lineage no EnTAP, mas mantenha os hits bacterianos numa tabela separada — pode ser achado biológico relevante (simbiose), não erro.
40. **Por que meu merge de tabela de anotação deu resultado biologicamente incoerente sem nenhum erro de execução?** Provável bug de merge por posição/índice de coluna em vez de nome — sempre faça merge por chave nomeada (`on="query"`) e valide uma amostra manualmente depois.
41. **DIAMOND é quantas vezes mais rápido que BLAST?** Até 20.000× mais rápido que BLASTX em reads curtas (não "~1000×", cifra que se refere a benchmarks específicos de blastp proteína-proteína).
42. **`--query-cover`/`--subject-cover` no DIAMOND — por que importam?** Evitam anotar por hits parciais (domínio compartilhado, proteína diferente), especialmente crítico em famílias multigênicas.

## Módulo 5 — Quantificação, DE e enriquecimento

43. **Como faço enriquecimento GO/KEGG se minha espécie não tem OrgDb no Bioconductor?** Construa `TERM2GENE`/`TERM2NAME` manualmente a partir do output do eggNOG-mapper e use `enricher()`/`GSEA()` genéricos do clusterProfiler, em vez de `enrichGO()`/`enrichKEGG()` com `OrgDb=`.
44. **Qual é o erro mais comum em enriquecimento funcional de não-modelo?** Usar um universo/background errado — por exemplo, o genoma inteiro de uma espécie modelo próxima, em vez do conjunto de genes efetivamente testados no próprio experimento.
45. **DESeq2 ou edgeR?** Não há um vencedor claro; edgeR tende a ser mais sensível (menos conservador) com poucas réplicas (n=2-3). A escolha usual do laboratório é DESeq2, pela integração nativa com tximport/lfcShrink.
46. **Quando usar limma-voom em vez de DESeq2/edgeR?** Quando o desenho é complexo (múltiplos fatores/blocos) e há amostras suficientes (tipicamente >10/grupo) para a transformação voom estimar bem a tendência média-variância.
47. **Corset ou tximport para agrupar transcritos em "genes"?** tximport é o caminho mais direto quando já existe o `gene_trans_map` do Trinity; Corset é preferido quando não há esse mapa ou se quer detectar redundância não capturada pelo montador.
48. **Resultado de enriquecimento com poucas réplicas é confiável?** Trate como hipótese a validar (ex.: qRT-PCR), não como conclusão definitiva — evidência recente mostra baixa replicabilidade de DE/enriquecimento em coortes pequenas.

## Módulo 6 — Figuras, interpretação e escrita do paper

49. **Que resolução/formato usar para figuras de submissão?** ≥300 dpi para raster (PNG/TIFF); sempre que possível, exportar também em SVG/PDF (vetorial) para arte final.
50. **É obrigatório usar Docker/Singularity/Nextflow para publicar?** Não há mandato formal identificado em revistas como Genome Biology/GigaScience/BMC Bioinformatics, mas é uma prática crescentemente esperada por revisores — `environment.yml` + registro de versões é o mínimo aceitável.
51. **O que preciso fazer antes de depositar minha montagem como TSA?** Depositar reads brutas no SRA primeiro (com BioProject/BioSample), rodar triagem de contaminação (FCS-GX), e seguir as regras de comprimento mínimo (200 pb) e tratamento de N.
52. **GitHub é suficiente para reprodutibilidade de código?** Não como referência citável permanente — arquive a versão final usada no paper no Zenodo, que gera um DOI imutável.

## Módulo 6b — Ortologia e filogenia (avançado)

53. **Posso rodar OrthoFinder direto no output bruto do Trinity?** Melhor não. Isoformas redundantes inflam artificialmente o proteoma — uma montagem de novo pode gerar várias vezes mais "transcritos" do que genes reais. Use o output pós-CD-HIT, com a isoforma representativa selecionada por gene.
54. **IQ-TREE ou FastTree?** IQ-TREE para conjuntos gerenciáveis (dezenas a centenas de táxons) onde a acurácia do suporte de ramo importa; FastTree/VeryFastTree quando há milhares de táxons ou a árvore é apenas exploratória.
55. **Ortólogos e parálogos são a mesma coisa?** Não — ortólogos divergem por especiação (tendem a reter função equivalente); parálogos divergem por duplicação gênica (frequentemente divergem funcionalmente).

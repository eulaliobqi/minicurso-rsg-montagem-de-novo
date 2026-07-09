# Módulo 1 — Qualidade e Pré-processamento das Reads

## Conceitos-chave

### Por que esta etapa decide o resto do pipeline

Numa montagem de novo não há genoma de referência para "absorver" erro de sequenciamento, adaptador residual ou contaminação por rRNA. Tudo isso vira nó extra no grafo de De Bruijn que o Trinity (ou qualquer outro assembler baseado em k-mer) vai tentar interpretar como transcrito real. Diferente do RNA-Seq com referência, em que um aligner splice-aware consegue "absorver" boa parte do ruído via soft-clipping, na montagem de novo o lixo que entra na etapa de pré-processamento tende a sair como transcrito quimérico, fragmentado ou redundante na saída do Trinity. Por isso este módulo trata QC e trimming com mais rigor do que normalmente se aplicaria a um RNA-Seq de quantificação simples contra referência.

### Inspeção de qualidade: FastQC + MultiQC

FastQC (Andrews, 2010) continua sendo a ferramenta de referência para inspeção por amostra de dados Illumina, mesmo sem nunca ter tido uma publicação formal em periódico. É citada pela documentação oficial do Babraham Institute, e a versão atual estável é a 0.12.1 (lançada em março de 2023). Os módulos que merecem atenção real numa aula prática são:

- **Per base sequence quality (Phred score):** o gráfico mais consultado. Em reads Illumina de 150 pb é normal ver Q30+ nas primeiras ~100 posições e uma queda gradual nas últimas 20-30 pb: isso é degradação química do sequenciamento, não um problema de biblioteca. O corte de qualidade do trimming deve respeitar esse padrão, não zerar a read inteira por causa da cauda.
- **Adapter content:** módulo que detecta contaminação por adaptador Illumina (TruSeq, Nextera). Em bibliotecas de RNA-Seq de insetos com fragmentos curtos ou reads mais longas que o inserto, é comum ver o adaptador aparecer a partir de ~120 pb num kit de 150 ciclos, sinal de que houve "read-through" e o trimming de adaptador é obrigatório, não opcional.
- **Per sequence GC content:** compara a distribuição observada com uma normal teórica. Insetos costumam girar entre 35-40% de GC no transcriptoma; um segundo pico ou uma distribuição bimodal geralmente indica contaminação (rRNA, DNA genômico residual, ou organismo simbionte/endossimbionte, relevante para o caso de *M. spectabilis*, que carrega endossimbiontes conhecidos).
- **Sequence duplication levels:** em RNA-Seq, duplicação alta é esperada (genes muito expressos geram muitas cópias idênticas de fragmentos), diferente de WGS onde indicaria PCR duplicates problemáticos. Não interpretar esse módulo como "falha" sem contexto.
- **Overrepresented sequences:** útil para pegar contaminação por adaptador não identificado automaticamente, ou sequências de rRNA muito abundantes que "vazaram" da depleção laboratorial.

FastQC roda amostra por amostra; **MultiQC** (Ewels et al., 2016) agrega os relatórios de todas as amostras, e de outras ferramentas do pipeline (fastp, Trinity stats, Salmon), num único HTML navegável. Num experimento com 6+ amostras (o mínimo recomendado para DESeq2 com réplicas), abrir um FastQC por vez é inviável; o MultiQC é o que realmente permite decidir rapidamente "essa amostra está fora do padrão das outras".

### SeqKit — estatísticas rápidas como complemento

O FastQC informa qualidade posicional; o **SeqKit** (Shen et al., 2016; atualizado como SeqKit2 em 2024) responde perguntas estruturais que o FastQC não responde bem: número total de reads, distribuição de comprimento, N50 de reads (relevante se houver mistura de comprimentos), contagem de bases N, e conversões rápidas de formato. O comando mais usado no dia a dia é `seqkit stats -a`, que dá em segundos um resumo numérico (n° reads, comprimento min/max/médio, %GC, Q20/Q30) e serve de checagem cruzada rápida antes mesmo de abrir o relatório gráfico do FastQC, útil sobretudo para conferir que os arquivos R1/R2 pareados têm o mesmo número de reads antes de rodar qualquer trimming pareado.

### Trimming e filtragem: fastp

**fastp** (Chen et al., 2018) consolidou-se como a ferramenta padrão de facto para QC + trimming num único passo, substituindo a combinação antiga FastQC (diagnóstico) + Trimmomatic ou Cutadapt (correção) por um único binário em C++ multithread que faz as duas coisas e ainda gera relatório HTML/JSON antes/depois. Os parâmetros que a turma precisa entender, não apenas copiar:

- **Corte por qualidade (`--qualified_quality_phred`, `--unqualified_percent_limit`):** fastp usa uma janela deslizante por padrão, cortando a partir do ponto em que a qualidade média cai abaixo do limiar. Evita o erro clássico de cortar um comprimento fixo arbitrário que tanto desperdiça bases boas quanto deixa lixo passar.
- **Remoção de adaptador (`--detect_adapter_for_pe`):** para dados pareados, fastp consegue detectar o adaptador por sobreposição de overlap entre R1 e R2, sem precisar fornecer a sequência do adaptador manualmente. É mais robusto do que listas fixas de adaptador do Trimmomatic.
- **Comprimento mínimo (`--length_required`):** crítico para montagem de novo. Reads muito curtas pós-trimming (< 36 pb, por exemplo) geram k-mers pouco informativos que inflam artificialmente o grafo de De Bruijn sem contribuir informação real. Devem ser descartadas, não mantidas "para não perder dado".
- **Filtragem poly-G/poly-X (`--trim_poly_g`, `--trim_poly_x`):** poly-G é um artefato específico de sequenciadores Illumina com química de dois canais (NovaSeq, NextSeq); nesses equipamentos, ausência de sinal é lida como G, então o final de reads curtas ou de baixa qualidade vira uma cauda de G's espúria. Vale confirmar qual plataforma gerou os dados antes de ativar esse filtro por padrão.

Em 2025 saiu o **fastp 1.0** (Chen, 2025, iMeta), primeira atualização "maior" da ferramenta desde 2018, com relatório HTML reformulado e suporte a processamento em lote de múltiplos FASTQ em paralelo. Reforça, não substitui, a posição do fastp como ferramenta consolidada.

### Correção de erros de sequenciamento: rcorrector

**rcorrector** (Song & Florea, 2015) corrige erros aleatórios de sequenciamento Illumina usando um grafo de De Bruijn de k-mers confiáveis, com limiar de confiança calculado localmente por posição, ao contrário de corretores desenhados para WGS, que assumem cobertura uniforme (premissa que RNA-Seq viola por definição, já que a cobertura varia com o nível de expressão do gene). rcorrector é a etapa central do **Oyster River Protocol** (MacManes, 2018), um protocolo consolidado e amplamente citado para montagem de novo que recomenda: rcorrector → remoção de reads "unfixable" → fastp → montagem com múltiplos assemblers/k-mers → merge das montagens. A lógica de usar rcorrector antes do trimming é que ele corrige a base errada em vez de simplesmente cortá-la, preservando mais informação (comprimento de read) do que um trimming agressivo por qualidade sozinho conseguiria.

### Remoção de rRNA: SortMeRNA

Mesmo com depleção de rRNA feita em laboratório (kits Ribo-Zero ou equivalente) ou seleção por poli-A, é comum sobrar 1-10% de reads de rRNA numa biblioteca de RNA-Seq. Esse resíduo, se não filtrado, infla os recursos computacionais gastos na montagem com transcritos de rRNA que ninguém vai analisar, e pode até fragmentar transcritos reais por competição de k-mers muito abundantes e conservados. O **SortMeRNA** (Kopylova, Noé & Touzet, 2012) faz esse filtro de forma rápida contra bancos de referência de rRNA (SILVA, Rfam), recomendado principalmente quando o protocolo de biblioteca não incluiu depleção específica, quando o organismo é filogeneticamente distante dos bancos de referência usados na depleção comercial (comum em insetos não-modelo), ou quando o QC pós-trimming ainda mostra picos de overrepresented sequences batendo com rRNA.

### Normalização in silico de leituras

A normalização in silico, implementada no Trinity via `insilico_read_normalization` (inspirada no algoritmo diginorm), reduz a profundidade de cobertura excessiva de regiões muito expressas mantendo a complexidade k-mer do conjunto. O objetivo é permitir que o assembler rode em datasets muito grandes sem gastar RAM/tempo processando redundância que não agrega informação nova sobre a estrutura do transcrito. Ajuda quando o dataset é muito grande (a literatura fala tipicamente em datasets acima de ~200 milhões de reads, onde o custo computacional de montar tudo sem normalizar se torna proibitivo em hardware comum) e a cobertura é muito desigual entre genes. Prejudica quando o objetivo é detectar isoformas raras ou genes de baixíssima expressão: como a normalização é probabilística e baseada em cobertura de k-mer, ela pode descartar desproporcionalmente reads de transcritos raros, cuja cobertura já era baixa antes da normalização. Na prática, isso significa risco maior de colapsar ou perder isoformas alternativas pouco expressas exatamente no tipo de gene que costuma interessar mais numa busca por candidatos biológicos específicos (o que é relevante para os dois organismos-caso deste curso, ambos com objetivos de descoberta de genes/vias específicas, não um catálogo exaustivo do transcriptoma).

### Ordem correta do pipeline e otimização de recursos

A ordem recomendada e testada pelo Oyster River Protocol e por guias recentes de montagem de novo é:

1. FastQC + MultiQC (diagnóstico do dado bruto)
2. rcorrector (correção de erro, opera melhor em reads ainda não cortadas)
3. Remoção de reads "unfixable" que o rcorrector sinaliza mas não consegue corrigir (script auxiliar do próprio Harvard Informatics/Oyster River Protocol)
4. fastp (trimming de qualidade + adaptador + filtro poly-G, com comprimento mínimo adequado)
5. SortMeRNA, se aplicável (remoção de rRNA residual)
6. FastQC + MultiQC novamente sobre o resultado (reavaliação pós-trimming, etapa que costuma ser pulada e não deveria)
7. Normalização in silico, só se o volume de dados justificar

Em termos de recursos no servidor (32 cores, RTX 5070 Ti): fastp e rcorrector escalam bem com threads (`-w`/`-p` no fastp, `-t` no rcorrector). Usar algo entre 8-16 threads por amostra processada em paralelo costuma ser mais eficiente do que jogar todos os 32 cores numa única amostra, já que a maior parte do tempo dessas ferramentas é I/O de disco, não CPU pura. SortMeRNA é o mais pesado em RAM por carregar o índice do banco de referência inteiro em memória; vale rodar com `--idx-dir` apontando para um índice pré-construído e reaproveitado entre amostras, em vez de reconstruir o índice a cada execução.

## Prática usual vs. estado da arte (2024-2026)

| Etapa | Ferramenta/prática usual | Estado da arte 2024-2026 | Fonte (DOI) |
|---|---|---|---|
| QC de leitura bruta | FastQC (por amostra) + MultiQC (agregação) | Sem substituto consolidado identificado nas buscas; FastQC segue como padrão (v0.12.1, 2023) e MultiQC segue em desenvolvimento ativo (v1.34/1.35, 2026). O próprio fastp já embute QC antes/depois no relatório HTML, reduzindo a dependência exclusiva do FastQC em pipelines modernos, mas não o substitui como ferramenta de inspeção universal comparável entre etapas diferentes | FastQC: sem DOI formal (Andrews, 2010, não peer-reviewed) · MultiQC: 10.1093/bioinformatics/btw354 |
| Estatísticas estruturais de sequência | Contagem manual / scripts próprios | SeqKit consolidado (2016) e expandido como SeqKit2 em 2024, dobrando o número de subcomandos (19→38); recomendado como complemento de rotina ao FastQC, não como opcional | SeqKit: 10.1371/journal.pone.0163962 · SeqKit2: 10.1002/imt2.191 |
| Trimming/adaptador | fastp (desde 2018, substituiu a dupla Trimmomatic/Cutadapt) | fastp segue sendo o padrão dominante; versão 1.0 (2025) reforça a posição sem indicar substituto emergente nos periódicos pesquisados. Há, porém, um debate relevante e bem estabelecido (não é "novidade 2024-2026", mas segue citado): para quantificação de expressão contra genoma de referência, trimming agressivo pode ser desnecessário ou até neutro, pois aligners modernos absorvem adaptador/baixa qualidade via soft-clipping. Esse achado **não se aplica** a montagens de novo, onde não há aligner fazendo esse papel: o assembler baseado em k-mer é sensível a erro residual de um jeito que o aligner não é | fastp: 10.1093/bioinformatics/bty560 · fastp 1.0: 10.1002/imt2.70078 · debate trimming pós-mapeamento: 10.1093/nargab/lqaa068 |
| Correção de erro de sequenciamento | rcorrector, dentro do Oyster River Protocol | Não encontrei, nas buscas realizadas, um sucessor publicado em 2024-2026 que tenha deslocado o rcorrector como recomendação de consenso para montagem de novo; ele segue citado como etapa recomendada em guias de 2022 e 2024 | rcorrector: 10.1186/s13742-015-0089-y · Oyster River Protocol: 10.7717/peerj.5428 · guia 2022: 10.1093/bib/bbab563 · guia 2024: 10.1186/s12983-024-00538-y |
| Remoção de rRNA | SortMeRNA | Ferramenta segue em desenvolvimento ativo (release v6.0 disponível via GitHub em 2024-2025; bioconda ainda distribui v4.3.7, conda-forge distribui até v5.0), sem alternativa nova de peso identificada nas buscas; continua sendo o padrão de facto | 10.1093/bioinformatics/bts611 |
| Normalização in silico | `insilico_read_normalization` do Trinity (algoritmo tipo diginorm) | Guia de 2024 (Jackson et al.) reafirma o uso como recomendado apenas para datasets muito grandes, não como etapa universal. Não encontrei nas buscas um consenso novo defendendo abandono da prática, mas o alerta sobre risco de distorcer isoformas raras é recorrente na literatura de metodologia desde a publicação original do Trinity e é repetido nos guias mais recentes | Trinity/normalização: 10.1038/nprot.2013.084 · guia 2024: 10.1186/s12983-024-00538-y |

## Recomendação para o curso

Para os dois organismos-caso deste curso (*A. gemmatalis* e *M. spectabilis*, ambos sem genoma de referência e com objetivo de descoberta de genes específicos: toxinas/endossimbiontes num caso, proteases digestivas no outro), a recomendação é seguir essencialmente o esqueleto do Oyster River Protocol, sem inflar o pipeline com etapas não justificadas pelo volume de dados real:

1. **FastQC + MultiQC** no dado bruto. A decisão de ir ou não adiante depende disso, não é opcional.
2. **rcorrector**, sempre: é barato computacionalmente e reduz erro residual antes do trimming, preservando mais comprimento de read do que cortar tudo por qualidade.
3. **Script de remoção de reads unfixable** do próprio conjunto de ferramentas do Harvard Informatics (companheiro usual do rcorrector).
4. **fastp** com `--detect_adapter_for_pe`, `--length_required 36` (ajustável conforme comprimento original das reads), e `--trim_poly_g` se a plataforma for NovaSeq/NextSeq.
5. **SortMeRNA** apenas se o QC pós-trimming (etapa 6) mostrar sinal de rRNA residual. Não rodar "por precaução" sem evidência: custa tempo e RAM sem necessidade quando a depleção laboratorial já funcionou.
6. **FastQC + MultiQC de novo**, agora sobre as reads processadas. É a etapa mais esquecida e a que efetivamente confirma se o trimming resolveu o problema identificado na etapa 1.
7. **Normalização in silico** apenas se o dataset ultrapassar a faixa em que a RAM do servidor (verificar com `free -h` antes) vira gargalo real para o Trinity. Não aplicar por padrão, justamente pelo risco de comprometer isoformas raras, que é exatamente o tipo de achado que pode interessar num transcriptoma de descoberta como os dois organismos-caso.

## Aplicação aos organismos-caso

Nem *A. gemmatalis* nem *M. spectabilis* têm, na memória deste projeto, métricas de QC de reads brutas registradas (as montagens Trinity existentes partiram de execuções anteriores cujos relatórios FastQC não foram arquivados). Por isso os valores abaixo são ilustrativos de padrões tipicamente observados em bibliotecas de RNA-Seq Illumina de insetos, não números reais desses dois datasets, e devem ser apresentados como tal em aula.

Ao abrir um relatório FastQC de uma biblioteca de inseto na aula, os pontos a destacar são:

- **%GC**: insetos costumam ficar na faixa de 35-40%. Um pico secundário fora dessa faixa no relatório de *M. spectabilis*, por exemplo, seria consistente com a presença de endossimbiontes já confirmados na montagem final desse projeto (hits tipo *Sulcia*/*Sodalis*-like). Ou seja: parte da "contaminação" de rRNA/DNA bacteriano detectável no FastQC bruto não é erro experimental, é biologia real do sistema hospedeiro-simbionte que só será separada mais adiante, na anotação taxonômica dos transcritos montados.
- **Adapter content**: bibliotecas com inserto curto (comum em RNA degradado ou em tecidos difíceis de extrair, como glândula salivar de *M. spectabilis*) tendem a mostrar sinal de adaptador mais cedo no comprimento da read, reforçando a necessidade de detecção automática de adaptador (fastp) em vez de assumir um kit fixo.
- **Per base sequence quality**: quando a queda de qualidade na cauda é abrupta (não gradual), pode indicar problema técnico da corrida de sequenciamento, não apenas o desgaste químico esperado. Vale conferir se todas as amostras do mesmo lote mostram o mesmo padrão (usando o MultiQC) antes de decidir um parâmetro de corte único para todas.
- **Overrepresented sequences**: em transcriptomas de descoberta como estes dois, sequências superrepresentadas que batem com rRNA merecem checagem via SortMeRNA; sequências superrepresentadas que batem com um transcrito real e biologicamente esperado (por exemplo, uma protease digestiva altamente expressa em intestino de lagarta) não são problema, são o próprio sinal biológico de interesse. É importante que a turma aprenda a diferenciar as duas situações antes de sair filtrando "tudo que é overrepresented".

## Erros comuns

- **Trimming agressivo demais**: cortar por um Phred alto demais (ex.: Q30 mínimo por base) ou impor comprimento mínimo alto demais gera perda desnecessária de dado e reduz artificialmente a profundidade de cobertura de genes pouco expressos, justamente os que mais interessam num transcriptoma de descoberta.
- **Esquecer de reavaliar o QC pós-trimming**: rodar fastp e seguir direto para a montagem sem abrir o MultiQC de novo é o erro mais comum e mais barato de evitar. Sem essa checagem não há como saber se o trimming resolveu o problema identificado no QC bruto ou se só reduziu o volume de dado sem resolver nada.
- **Rodar SortMeRNA sem necessidade comprovada**: gasta RAM e tempo de execução relevantes; só justificado quando o QC (bruto ou pós-trimming) mostra evidência real de rRNA residual.
- **Aplicar normalização in silico por padrão, "porque o Trinity oferece"**: pode distorcer a representação de isoformas de baixa expressão sem necessidade real, se o dataset já cabe confortavelmente na RAM disponível sem normalizar.
- **Confundir duplicação alta no FastQC com problema de PCR**: em RNA-Seq isso é esperado para genes altamente expressos; tratar como sinal de falha de biblioteca (e por exemplo rodar deduplicação agressiva) remove sinal biológico real, não ruído.
- **Não confirmar a plataforma de sequenciamento antes de ativar `--trim_poly_g`**: esse filtro é específico de química de dois canais (NovaSeq/NextSeq); ativá-lo cegamente em dados de outras plataformas pode cortar informação real sem necessidade.
- **Ignorar o pareamento R1/R2 depois do trimming**: fastp em modo paired-end mantém a sincronia por padrão, mas scripts customizados de filtragem adicional (como remoção de reads "unfixable" do rcorrector) precisam necessariamente processar R1 e R2 em conjunto. Descartar uma mata do par sem descartar a outra quebra a montagem pareada mais adiante.

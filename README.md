# RNA-Seq Differential Expression Analysis in Coronary Artery Disease

Este projeto teste, analisa perfis de expressão gênica de dados de RNA-Seq para identificar genes diferencialmente expressos em pacientes com Doença Arterial Coronariana (DAC) usando R.

<img width="960" height="86" alt="rna_seq_dac_project_flow" src="https://github.com/user-attachments/assets/2fc846fa-9f34-4273-979b-5c11786f9a9c" />

## Objetivo

Execute um pipeline completo para análise de RNA-Seq:
- Importação e limpeza de dados
- Análise exploratória de dados
- Expressão diferencial com DESeq2
- Análise de enriquecimento funcional (GO, KEGG)
- Visualização interativa com Shiny
- Relatórios automatizados com RMarkdown

## Estrutura Básica

```
rna_seq_dac_project/
│
├── data/
│   ├── raw/
│   └── processed/
│
├── scripts/
│   ├── 01_load_and_clean.R
│   ├── 02_eda_visualization.R
│   ├── 03_deseq_analysis.R
│   ├── 04_pathway_enrichment.R
│   └── 05_generate_report.R
│
├── shiny_app/
│   └── app.R
│
├── results/
│
├── reports/
│   └── analysis_report.Rmd
│
├── .gitignore
├── README.md
└── renv.lock / renv/
```

## Estrutura do Projeto

```
scripts/ scripts Core R
dados/raw/ matriz de contagem bruta de RNA-Seq e metadados
dados/processados/ dados limpos para modelagem
resultados/ saída da análise
relatórios/ relatórios finais (Rmd e HTML)
shiny_app/ painel interativo
```

## Pacotes necessários

- DESeq2
- BiocManager
- dplyr, tidyr, readr
- ggplot2, pheatmap, plotly
- clusterProfiler, enrichR
- shiny, flexdashboard
- rmarkdown

## Execução

Execute scripts sequencialmente a partir de `scripts/` ou carregue como um pipeline.

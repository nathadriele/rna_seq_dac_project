# RNA-Seq Differential Expression Analysis in Coronary Artery Disease

This project analyzes gene expression profiles from RNA-Seq data to identify differentially expressed genes in patients with Coronary Artery Disease (CAD) using R.

<img width="960" height="86" alt="rna_seq_dac_project_flow" src="https://github.com/user-attachments/assets/2fc846fa-9f34-4273-979b-5c11786f9a9c" />

## Objective

Perform a complete pipeline for RNA-Seq analysis:
- Data import and cleaning
- Exploratory data analysis
- Differential expression with DESeq2
- Functional enrichment analysis (GO, KEGG)
- Interactive visualization with Shiny
- Automated reporting with RMarkdown

## Basic Structure
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
└── renv.lock / renv/ (opcional)
```

## Project Structure

```
scripts/                 Core R scripts  
data/raw/                Raw RNA-Seq count matrix and metadata  
data/processed/          Cleaned data for modeling  
results/                 Output from analysis  
reports/                 Final reports (Rmd and HTML)  
shiny_app/               Interactive dashboard  
```

## Required Packages

- DESeq2
- BiocManager
- dplyr, tidyr, readr
- ggplot2, pheatmap, plotly
- clusterProfiler, enrichR
- shiny, flexdashboard
- rmarkdown

## Running

Run scripts sequentially from `scripts/` or load as a pipeline.

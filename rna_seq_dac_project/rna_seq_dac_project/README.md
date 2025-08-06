# RNA-Seq Differential Expression Analysis in Coronary Artery Disease

This project analyzes gene expression profiles from RNA-Seq data to identify differentially expressed genes in patients with Coronary Artery Disease (CAD) using R.

## Objective

Perform a complete pipeline for RNA-Seq analysis:
- Data import and cleaning
- Exploratory data analysis
- Differential expression with DESeq2
- Functional enrichment analysis (GO, KEGG)
- Interactive visualization with Shiny
- Automated reporting with RMarkdown

## Project Structure

scripts/                 Core R scripts  
data/raw/                Raw RNA-Seq count matrix and metadata  
data/processed/          Cleaned data for modeling  
results/                 Output from analysis  
reports/                 Final reports (Rmd and HTML)  
shiny_app/               Interactive dashboard  

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
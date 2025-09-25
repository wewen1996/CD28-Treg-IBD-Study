# CD28-Treg-IBD-Study

## Project Overview

CD28-Treg-IBD-Study is a multi-omics analysis project investigating the role of CD28+ regulatory T cells (Treg) in inflammatory bowel disease (IBD). This project integrates genetics, single-cell RNA sequencing, metabolomics, and microbiomics data, aiming to reveal the role of CD28+ Treg cells in the pathogenesis of IBD and explore potential therapeutic targets.

## Project structure
```
CD28-Treg-IBD-Study/
│
├── README.md                         
├── CITATION.cff                        
├── LICENSE                            
│
├── 01_Genetic_Analysis/                
│   ├── MR_Analysis/                    
│   │   ├── multi_exposure_single_outcome.R      
│   │   ├── single_exposure_multi_outcome.R     
│   │   └── utility_functions.R                  
│   ├── MVMR_Analysis/                  
│   ├── SMR_Analysis/                   
│   ├── COLOC_Analysis/                 
│   │   └── colorR.R                   
│   ├── LDSC_Analysis/                 
│   │   └── LDSC.R                     
│   └── Meta_Analysis/                 
│       └── multi_source_GWAS_meta.R    
│
├── 02_scRNA_Seq_Analysis/             
    ├── preprocessing/                 
    │   └── single_cell_preprocessing.R
    ├── clustering_annotation/         
    │   └── scRNA_analysis_main.R      
    ├── subpopulation_analysis/         
    │   └── Tcell_subset_analysis.R    
    └── visualization/                  
        └── scRNA_visualization.R      

```

## Main analysis module

### 1. Genetic analysis (core part)

This module focuses on studying the genetic association between CD28+Treg cells and IBD, including:

- **Mendelian randomization analysis**: Assessing the causal relationship between CD28+Treg cell characteristics and IBD using GWAS summary data
- **Multivariate MR analysis**: By controlling the influence of other factors, it can more accurately assess causal relationships.
- **Colocalization Analysis**: Exploring the genetic variations shared between the characteristics of CD28+Treg cells and IBD
- **LD Score Regression**: Evaluate the heritability and genetic correlation between CD28+Treg cell characteristics and IBD
- **Meta-analysis**: Integrating results from multiple GWAS studies to improve statistical power

### 2. Single-cell RNA sequencing analysis

This module conducts an in-depth analysis of the transcriptomic characteristics of CD28+Treg cells, including:

- **Data Preprocessing**: Quality control, normalization, batch effect correction
- **Clustering and Annotation**: Cell type identification and annotation
- **Subgroup Analysis**: In-depth analysis of T cell subsets, particularly CD28+Treg cells
- **Visualization**: High-quality visualization results, including UMAP, heatmaps, volcano plots, etc.

## Installation Instructions

### Installation of dependent packages

```r
# Install the necessary R packages
install.packages(c("tidyverse", "devtools", "BiocManager"))

# Install Bioconductor packages
BiocManager::install(c("SummarizedExperiment", "DESeq2", "edgeR", "limma", "org.Hs.eg.db", 
                       "clusterProfiler", "ReactomePA", "DOSE", "enrichplot", "phyloseq", 
                       "microbiome", "mixOmics", "WGCNA", "MetaboAnalystR", "FELLA", 
                       "SingleCellExperiment", "scater", "scran", "batchelor", "celldex", 
                       "SingleR", "DoubletFinder", "scDblFinder"))

# Install CRAN packages
install.packages(c("Seurat", "ggplot2", "ggpubr", "patchwork", "RColorBrewer", "viridis", 
                   "dplyr", "tidyr", "here", "glue", "data.table", "corrplot", "pheatmap", 
                   "ggrepel", "cowplot", "gridExtra", "scales", "factoextra", "FactoMineR", 
                   "network", "sna", "ggnetwork", "intergraph", "igraph", "tidygraph", 
                   "ggraph", "circlize", "ComplexHeatmap", "LDlinkR", "TwoSampleMR", 
                   "MendelianRandomization", "MRInstruments", "MRBase", "mr.raps", 
                   "MVMR", "coloc", "genetics.binaRies", "GenomicSEM", "ldsc"))

# Install GitHub packages
devtools::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("MRCIEU/MRInstruments")
devtools::install_github("MRCIEU/mr.raps")
devtools::install_github("MRCIEU/MVMR")
devtools::install_github("chr1swallace/coloc")
devtools::install_github("andrewparkermorgan/argparse")
devtools::install_github("stephenslab/gtexresults")
devtools::install_github("bcellmanlab/scRNAseqAnalysis")
```

### Data preparation

1. Put the original data in the `data/raw/` directory
2. Put the metadata in the `data/metadata/` directory.
3. Place the external dataset in the `data/external/` directory

## Operation Guide

### Genetic analysis

```bash
# Perform Mendelian randomization analysis
Rscript 01_Genetic_Analysis/MR_Analysis/multi_exposure_single_outcome.R
Rscript 01_Genetic_Analysis/MR_Analysis/single_exposure_multi_outcome.R

# Run colocalization analysis
Rscript 01_Genetic_Analysis/COLOC_Analysis/colorR.R

# Run LD score regression analysis
Rscript 01_Genetic_Analysis/LDSC_Analysis/LDSC.R

# Run a Meta-analysis
Rscript 01_Genetic_Analysis/Meta_Analysis/multi_source_GWAS_meta.R
```

### Single-cell RNA sequencing analysis

```bash
# Run data preprocessing
Rscript 02_scRNA_Seq_Analysis/preprocessing/single_cell_preprocessing.R

# Run clustering and annotation
Rscript 02_scRNA_Seq_Analysis/clustering_annotation/scRNA_analysis_main.R

# Run T cell subset analysis
Rscript 02_scRNA_Seq_Analysis/subpopulation_analysis/Tcell_subset_analysis.R

# Run visualization
Rscript 02_scRNA_Seq_Analysis/visualization/scRNA_visualization.R
```

```

```

## Main results

### Results of genetic analysis

- **Mendelian randomization analysis**: A significant causal relationship was found between CD28+Treg cell characteristics and the risk of IBD.
- **Colocalization analysis**: Genetic loci shared by CD28+Treg cell characteristics and IBD were identified.
- **LD score regression**: A significant genetic correlation was found between CD28+Treg cell characteristics and IBD.
- **Meta-analysis**: Integrating multiple GWAS studies improved the statistical power of CD28+Treg cell-related loci.

### Results of Single-Cell RNA Sequencing Analysis

- **Cell Type Identification**: Successfully identified various immune cell types, including CD4+ T cells, CD8+ T cells, B cells, myeloid cells, etc.
- **T Cell Subset Analysis**: Conducted in-depth analysis of T cell subsets, especially the transcriptome characteristics of CD28+ Treg cells.
- **Differential Expression Analysis**: Discovered significant changes in gene expression of CD28+ Treg cells in IBD patients.
- **Functional Enrichment Analysis**: Found that the differentially expressed genes in CD28+ Treg cells are mainly involved in biological processes such as immune regulation, cell activation, and inflammatory response.


## Key Findings

1. **Causal Relationship between CD28+Treg Cells and IBD**: Genetic analysis indicates that abnormal function of CD28+Treg cells may be one of the pathogenic factors of IBD.
2. **Transcriptomic Characteristics of CD28+Treg Cells**: Single-cell RNA sequencing analysis reveals unique transcriptomic characteristics of CD28+Treg cells in IBD patients.
3. **Functional Heterogeneity of CD28+Treg Cells**: It is found that CD28+Treg cells exhibit functional heterogeneity, and different subsets may play different roles in IBD.


## Discussion and Outlook

### Research Significance

This study explores the role of CD28+Treg cells in IBD through multi-omics integration analysis, providing a new perspective for understanding the pathogenesis of IBD and potential targets for the development of new therapeutic strategies.

### Future Research Directions

1. **Expand Sample Size**: Increase the sample size to improve the statistical power and reliability of the results.
2. **Functional Verification**: Verify the function of CD28+Treg cells in IBD through in vitro and in vivo experiments.
3. **Mechanism Research**: Conduct in-depth research on the molecular mechanisms by which CD28+Treg cells regulate IBD.
4. **Development of Therapeutic Targets**: Based on the findings of this study, develop therapeutic strategies targeting CD28+Treg cells.
5. **Multi-omics Integration**: Further integrate genomics, transcriptomics, proteomics, metabolomics, and microbiome data to comprehensively understand the pathogenesis of IBD.

```

## License

This project is licensed under the MIT License. For details, please refer to the LICENSE file.

## Contact Information

If you have any questions or suggestions, please contact:

- Project Leader: Wenwen Wang and Chen Li
- Email: kjjydawangww@163.com
- GitHub: https://github.com/wewen1996/

## Update Log

- **2024-09-24**: Project initialization, creation of main directory structure and documents
- **2024-09-25**: Addition of main scripts for genetic analysis and single-cell RNA sequencing analysis

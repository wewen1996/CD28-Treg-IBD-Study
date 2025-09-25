# CD28-Treg-IBD-Study 项目依赖包安装脚本
# 此脚本包含了项目中所有分析模块所需的R包

# 检查并安装Bioconductor管理包
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 安装CRAN包
cran_packages <- c(
  # 基础包
  "tidyverse", "data.table", "magrittr", "here", "glue", "yaml", "config",
  "devtools", "remotes", "usethis", "git2r",
  
  # 数据处理
  "dplyr", "tidyr", "purrr", "tibble", "readr", "stringr", "forcats",
  
  # 可视化
  "ggplot2", "ggpubr", "patchwork", "cowplot", "pheatmap", "RColorBrewer",
  "viridis", "corrplot", "igraph", "ggraph", "ComplexHeatmap",
  
  # 统计分析
  "limma", "edgeR", "DESeq2", "lme4", "nlme", "car", "multcomp",
  "emmeans", "lsmeans", "survival", "survminer",
  
  # 遗传学分析 - MR相关
  "TwoSampleMR", "MRInstruments", "MRBase", "MendelianRandomization",
  "mr.raps", "RadialMR", "MVMR", "smr", "coloc", "LDSC",
  
  # 单细胞分析
  "Seurat", "SingleCellExperiment", "scater", "scran", "monocle3",
  "celldex", "SingleR", "scDblFinder", "DoubletFinder",
  
  # 代谢组学分析
  "MetaboAnalystR", "xcms", "CAMERA", "metagenomeSeq",
  
  # 微生物组分析
  "phyloseq", "microbiome", "ANCOMBC", "Maaslin2", "corncob",
  
  # 网络分析
  "WGCNA", "mixOmics", "igraph", "network", "sna",
  
  # 其他
  "knitr", "rmarkdown", "markdown", "bookdown", "kableExtra",
  "sessioninfo", "sessioninfo"
)

# 安装缺失的CRAN包
new_cran_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_cran_packages)) install.packages(new_cran_packages)

# 安装Bioconductor包
bioc_packages <- c(
  # 基础Bioconductor包
  "Biobase", "BiocGenerics", "GenomicRanges", "IRanges", "S4Vectors",
  
  # 测序数据分析
  "SummarizedExperiment", "DESeq2", "edgeR", "limma", "voom",
  
  # 单细胞分析
  "SingleCellExperiment", "scater", "scran", "batchelor", "scuttle",
  
  # 代谢组学分析
  "xcms", "CAMERA", "metabomxtr", "MetaboCoreUtils",
  
  # 遗传学分析
  "GenomeInfoDb", "GenomicFeatures", "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db", "biomaRt",
  
  # 其他
  "BiocStyle", "ReportingTools", "InteractiveComplexHeatmap"
)

# 安装缺失的Bioconductor包
new_bioc_packages <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
if(length(new_bioc_packages)) BiocManager::install(new_bioc_packages)

# 安装GitHub上的包
github_packages <- c(
  "MRCIEU/TwoSampleMR",
  "MRCIEU/MRInstruments",
  "MRCIEU/MRBase",
  "chr1swallace/coloc",
  "GenomicSEM/GenomicSEM",
  "immunogenomics/harmony",
  "satijalab/seurat-data",
  "cole-trapnell-lab/monocle3",
  "statgen/MR-PRESSO",
  "Wenan/MVMR"
)

# 安装GitHub包
for (pkg in github_packages) {
  pkg_name <- strsplit(pkg, "/")[[1]][2]
  if (!pkg_name %in% installed.packages()[,"Package"]) {
    devtools::install_github(pkg)
  }
}

# 检查并加载所有包
load_packages <- function(packages) {
  lapply(packages, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      warning(paste("Package", pkg, "could not be loaded"))
    }
  })
}

# 加载所有包
load_packages(c(cran_packages, bioc_packages))

# 打印会话信息
sessionInfo()

# 打印成功消息
cat("所有必需的R包已成功安装和加载！\n")

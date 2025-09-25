#' CD28-Treg-IBD-Study: 单细胞RNA测序数据预处理
#' 
#' 本脚本用于单细胞RNA测序数据的预处理，包括质量控制、归一化、批次效应校正
#' 和高可变基因鉴定等步骤。

# 加载必要的包
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)
library(DropletUtils)
library(batchelor)
library(scDblFinder)
library(DoubletFinder)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(glue)

# 设置输出目录
output_dir <- here("results", "scRNA_Seq_Analysis", "preprocessing")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------
# 1. 数据加载
# ----------------------

cat("开始数据加载...\n")

# 注意：在实际分析中，应从data目录加载真实数据
# 这里使用示例数据进行演示

# 创建示例数据
create_sample_data <- function(n_cells = 1000, n_genes = 5000, n_samples = 3) {
  # 创建细胞ID
  cell_ids <- paste0("cell_", 1:(n_cells * n_samples))
  
  # 创建基因ID
  gene_ids <- paste0("gene_", 1:n_genes)
  
  # 创建样本信息
  sample_info <- data.frame(
    cell = cell_ids,
    sample = rep(paste0("sample_", 1:n_samples), each = n_cells),
    batch = rep(paste0("batch_", 1:2), length.out = n_cells * n_samples),
    group = rep(c("control", "treatment"), length.out = n_cells * n_samples)
  )
  
  # 创建表达矩阵
  set.seed(1234)
  counts_matrix <- matrix(
    rpois(n_cells * n_samples * n_genes, lambda = 5),
    nrow = n_genes,
    ncol = n_cells * n_samples,
    dimnames = list(gene_ids, cell_ids)
  )
  
  # 添加一些差异表达基因
  de_genes <- sample(gene_ids, 100)
  treatment_cells <- sample_info$group == "treatment"
  count_matrix[de_genes, treatment_cells] <- count_matrix[de_genes, treatment_cells] * 2
  
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(
    counts = count_matrix,
    meta.data = sample_info,
    project = "CD28-Treg-IBD"
  )
  
  return(seurat_obj)
}

# 创建示例数据
seurat_obj <- create_sample_data(n_cells = 1000, n_genes = 5000, n_samples = 3)

# 保存原始数据
saveRDS(seurat_obj, here(output_dir, "raw_seurat_object.rds"))

cat("数据加载完成，创建了包含", ncol(seurat_obj), "个细胞和", nrow(seurat_obj), "个基因的Seurat对象\n")

# ----------------------
# 2. 质量控制
# ----------------------

cat("\n开始质量控制...\n")

# 2.1 计算线粒体比例
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# 2.2 计算核糖体蛋白比例
rb_genes <- grep("^RP[SL]", rownames(seurat_obj), value = TRUE)
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, features = rb_genes)

# 2.3 计算红细胞比例
hb_genes <- grep("^HB[^(P)]", rownames(seurat_obj), value = TRUE)
if (length(hb_genes) > 0) {
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, features = hb_genes)
}

# 2.4 创建质量控制指标的数据框
qc_metrics <- data.frame(
  cell = colnames(seurat_obj),
  nCount_RNA = seurat_obj$nCount_RNA,
  nFeature_RNA = seurat_obj$nFeature_RNA,
  percent.mt = seurat_obj$percent.mt,
  percent.rb = seurat_obj$percent.rb,
  sample = seurat_obj$sample,
  batch = seurat_obj$batch,
  group = seurat_obj$group
)

if (length(hb_genes) > 0) {
  qc_metrics$percent.hb <- seurat_obj$percent.hb
}

# 2.5 保存质量控制指标
write.csv(qc_metrics, here(output_dir, "qc_metrics.csv"), row.names = FALSE)

# 2.6 质量控制可视化
# 创建质量控制指标的小提琴图
p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sample")
p2 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "group")

# 保存小提琴图
ggsave(
  filename = here(output_dir, "qc_violin_plots.png"),
  plot = cowplot::plot_grid(p1, p2, ncol = 1),
  width = 15,
  height = 12,
  dpi = 300
)

# 创建特征数与计数的散点图
p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample") +
  geom_smooth(method = "lm")

# 创建特征数与线粒体比例的散点图
p4 <- FeatureScatter(seurat_obj, feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by = "sample") +
  geom_smooth(method = "lm")

# 保存散点图
ggsave(
  filename = here(output_dir, "qc_scatter_plots.png"),
  plot = cowplot::plot_grid(p3, p4, ncol = 1),
  width = 10,
  height = 12,
  dpi = 300
)

# 2.7 设置质量控制阈值
# 基于数据分布设置阈值
nFeature_RNA_min <- 200
nFeature_RNA_max <- quantile(seurat_obj$nFeature_RNA, 0.95)
nCount_RNA_min <- 500
percent_mt_max <- 10

cat("质量控制阈值:\n")
cat(sprintf("- 最小特征数: %d\n", nFeature_RNA_min))
cat(sprintf("- 最大特征数: %d\n", nFeature_RNA_max))
cat(sprintf("- 最小计数: %d\n", nCount_RNA_min))
cat(sprintf("- 最大线粒体比例: %.1f%%\n", percent_mt_max))

# 2.8 应用质量控制
seurat_obj_filtered <- subset(
  seurat_obj,
  subset = nFeature_RNA > nFeature_RNA_min &
    nFeature_RNA < nFeature_RNA_max &
    nCount_RNA > nCount_RNA_min &
    percent.mt < percent_mt_max
)

# 2.9 统计过滤结果
cat("\n质量控制结果:\n")
cat(sprintf("- 过滤前细胞数: %d\n", ncol(seurat_obj)))
cat(sprintf("- 过滤后细胞数: %d\n", ncol(seurat_obj_filtered)))
cat(sprintf("- 过滤掉的细胞数: %d (%.1f%%)\n", 
            ncol(seurat_obj) - ncol(seurat_obj_filtered),
            (ncol(seurat_obj) - ncol(seurat_obj_filtered)) / ncol(seurat_obj) * 100))

# 2.10 保存过滤后的对象
saveRDS(seurat_obj_filtered, here(output_dir, "filtered_seurat_object.rds"))

# ----------------------
# 3.  doublet检测
# ----------------------

cat("\n开始doublet检测...\n")

# 3.1 使用scDblFinder检测doublet
set.seed(1234)
seurat_obj_filtered <- scDblFinder(seurat_obj_filtered, verbose = TRUE)

# 3.2 使用DoubletFinder检测doublet
# 计算pK值
sweep.res.list <- paramSweep_v3(seurat_obj_filtered, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# 选择最佳pK值
pk <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

# 计算doublet比例
annotations <- seurat_obj_filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075 * ncol(seurat_obj_filtered))  # 假设7.5%的doublet比例
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# 运行DoubletFinder
set.seed(1234)
seurat_obj_filtered <- doubletFinder_v3(
  seurat_obj_filtered,
  pN = 0.25,
  pK = pk,
  nExp = nExp_poi.adj,
  reuse.pANN = FALSE,
  sct = FALSE
)

# 3.3 合并doublet检测结果
# 获取scDblFinder结果
scdblfinder_doublets <- seurat_obj_filtered$scDblFinder.class == "doublet"

# 获取DoubletFinder结果
doubletfinder_doublets <- grepl("Doublet", colnames(seurat_obj_filtered)) %>%
  {colnames(seurat_obj_filtered)[.]} %>%
  {seurat_obj_filtered[[.]] == "Doublet"} %>%
  rowSums() > 0

# 合并结果
seurat_obj_filtered$doublet <- scdblfinder_doublets | doubletfinder_doublets

# 3.4 统计doublet检测结果
cat("Doublet检测结果:\n")
cat(sprintf("- scDblFinder检测到的doublet: %d (%.1f%%)\n", 
            sum(scdblfinder_doublets), sum(scdblfinder_doublets)/ncol(seurat_obj_filtered)*100))
cat(sprintf("- DoubletFinder检测到的doublet: %d (%.1f%%)\n", 
            sum(doubletfinder_doublets), sum(doubletfinder_doublets)/ncol(seurat_obj_filtered)*100))
cat(sprintf("- 合并检测到的doublet: %d (%.1f%%)\n", 
            sum(seurat_obj_filtered$doublet), sum(seurat_obj_filtered$doublet)/ncol(seurat_obj_filtered)*100))

# 3.5 去除doublet
seurat_obj_singlet <- subset(seurat_obj_filtered, subset = doublet == FALSE)

cat(sprintf("- 去除doublet后的细胞数: %d\n", ncol(seurat_obj_singlet)))
cat(sprintf("- 去除的doublet数: %d\n", ncol(seurat_obj_filtered) - ncol(seurat_obj_singlet)))

# 3.6 保存去除doublet后的对象
saveRDS(seurat_obj_singlet, here(output_dir, "singlet_seurat_object.rds"))

# ----------------------
# 4. 数据归一化
# ----------------------

cat("\n开始数据归一化...\n")

# 4.1 使用SCTransform进行归一化
# 这是Seurat v3及以上版本推荐的归一化方法
set.seed(1234)
seurat_obj_normalized <- SCTransform(
  seurat_obj_singlet,
  vars.to.regress = c("percent.mt", "nCount_RNA"),
  verbose = TRUE
)

# 4.2 保存归一化后的对象
saveRDS(seurat_obj_normalized, here(output_dir, "normalized_seurat_object.rds"))

cat("数据归一化完成\n")

# ----------------------
# 5. 批次效应校正
# ----------------------

cat("\n开始批次效应校正...\n")

# 5.1 检查是否有多个批次
if (length(unique(seurat_obj_normalized$batch)) > 1) {
  cat(sprintf("检测到%d个批次，进行批次效应校正...\n", length(unique(seurat_obj_normalized$batch))))
  
  # 5.2 使用Harmony进行批次效应校正
  # 首先进行PCA降维
  seurat_obj_normalized <- RunPCA(seurat_obj_normalized, verbose = FALSE)
  
  # 使用Harmony校正批次效应
  library(harmony)
  set.seed(1234)
  seurat_obj_corrected <- RunHarmony(
    seurat_obj_normalized,
    group.by.vars = "batch",
    verbose = TRUE
  )
  
  # 5.3 批次效应校正前后的可视化比较
  # 校正前的PCA图
  p5 <- DimPlot(seurat_obj_normalized, reduction = "pca", group.by = "batch") +
    ggtitle("批次效应校正前 (PCA)")
  
  # 校正后的Harmony图
  p6 <- DimPlot(seurat_obj_corrected, reduction = "harmony", group.by = "batch") +
    ggtitle("批次效应校正后 (Harmony)")
  
  # 保存批次效应校正前后的比较图
  ggsave(
    filename = here(output_dir, "batch_correction_comparison.png"),
    plot = cowplot::plot_grid(p5, p6, ncol = 2),
    width = 15,
    height = 6,
    dpi = 300
  )
  
  # 5.4 保存校正后的对象
  saveRDS(seurat_obj_corrected, here(output_dir, "batch_corrected_seurat_object.rds"))
  
  cat("批次效应校正完成\n")
} else {
  cat("未检测到多个批次，跳过批次效应校正\n")
  seurat_obj_corrected <- seurat_obj_normalized
  saveRDS(seurat_obj_corrected, here(output_dir, "batch_corrected_seurat_object.rds"))
}

# ----------------------
# 6. 高可变基因鉴定
# ----------------------

cat("\n开始高可变基因鉴定...\n")

# 6.1 使用SCTransform鉴定的高可变基因
# SCTransform已经鉴定了高可变基因，存储在variable.features中
hvg_sct <- VariableFeatures(seurat_obj_corrected)
cat(sprintf("SCTransform鉴定到%d个高可变基因\n", length(hvg_sct)))

# 6.2 使用传统方法鉴定高可变基因（作为比较）
# 首先创建一个临时对象，使用传统方法
seurat_obj_temp <- seurat_obj_singlet
seurat_obj_temp <- NormalizeData(seurat_obj_temp, verbose = FALSE)
seurat_obj_temp <- FindVariableFeatures(
  seurat_obj_temp,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)
hvg_vst <- VariableFeatures(seurat_obj_temp)
cat(sprintf("传统方法鉴定到%d个高可变基因\n", length(hvg_vst)))

# 6.3 比较两种方法鉴定的高可变基因
common_hvg <- intersect(hvg_sct, hvg_vst)
cat(sprintf("两种方法共有的高可变基因: %d (%.1f%%)\n", 
            length(common_hvg), length(common_hvg)/length(hvg_sct)*100))

# 6.4 可视化高可变基因
# 高可变基因热图
p7 <- VariableFeaturePlot(seurat_obj_corrected) +
  ggtitle("高可变基因 (SCTransform)")

# 前20个高可变基因的表达热图
top20_hvg <- head(hvg_sct, 20)
p8 <- DoHeatmap(seurat_obj_corrected, features = top20_hvg, size = 3) +
  NoLegend() +
  ggtitle("前20个高可变基因的表达热图")

# 保存高可变基因可视化结果
ggsave(
  filename = here(output_dir, "variable_features_visualization.png"),
  plot = cowplot::plot_grid(p7, p8, ncol = 1, rel_heights = c(1, 2)),
  width = 12,
  height = 15,
  dpi = 300
)

# 6.5 保存高可变基因列表
write.csv(data.frame(gene = hvg_sct), here(output_dir, "variable_features.csv"), row.names = FALSE)

# ----------------------
# 7. 数据降维和聚类
# ----------------------

cat("\n开始数据降维和聚类...\n")

# 7.1 降维
# 如果进行了批次效应校正，使用harmony结果进行降维
if ("harmony" %in% Reductions(seurat_obj_corrected)) {
  seurat_obj_corrected <- RunUMAP(seurat_obj_corrected, reduction = "harmony", dims = 1:30, verbose = FALSE)
  seurat_obj_corrected <- RunTSNE(seurat_obj_corrected, reduction = "harmony", dims = 1:30, verbose = FALSE)
} else {
  seurat_obj_corrected <- RunPCA(seurat_obj_corrected, verbose = FALSE)
  seurat_obj_corrected <- RunUMAP(seurat_obj_corrected, reduction = "pca", dims = 1:30, verbose = FALSE)
  seurat_obj_corrected <- RunTSNE(seurat_obj_corrected, reduction = "pca", dims = 1:30, verbose = FALSE)
}

# 7.2 聚类
seurat_obj_corrected <- FindNeighbors(seurat_obj_corrected, reduction = ifelse("harmony" %in% Reductions(seurat_obj_corrected), "harmony", "pca"), dims = 1:30, verbose = FALSE)
seurat_obj_corrected <- FindClusters(seurat_obj_corrected, resolution = 0.5, verbose = FALSE)

# 7.3 聚类结果可视化
# UMAP图（按聚类着色）
p9 <- DimPlot(seurat_obj_corrected, reduction = "umap", label = TRUE, label.size = 6) +
  ggtitle("细胞聚类 (UMAP)")

# UMAP图（按样本着色）
p10 <- DimPlot(seurat_obj_corrected, reduction = "umap", group.by = "sample") +
  ggtitle("样本分布 (UMAP)")

# UMAP图（按分组着色）
p11 <- DimPlot(seurat_obj_corrected, reduction = "umap", group.by = "group") +
  ggtitle("分组分布 (UMAP)")

# 如果有批次信息，按批次着色
if (length(unique(seurat_obj_corrected$batch)) > 1) {
  p12 <- DimPlot(seurat_obj_corrected, reduction = "umap", group.by = "batch") +
    ggtitle("批次分布 (UMAP)")
  
  # 保存UMAP可视化结果
  ggsave(
    filename = here(output_dir, "clustering_visualization.png"),
    plot = cowplot::plot_grid(p9, p10, p11, p12, ncol = 2),
    width = 15,
    height = 12,
    dpi = 300
  )
} else {
  # 保存UMAP可视化结果
  ggsave(
    filename = here(output_dir, "clustering_visualization.png"),
    plot = cowplot::plot_grid(p9, p10, p11, ncol = 2),
    width = 15,
    height = 8,
    dpi = 300
  )
}

# 7.4 统计聚类结果
cluster_stats <- table(seurat_obj_corrected$seurat_clusters, seurat_obj_corrected$group) %>%
  as.data.frame.matrix() %>%
  rownames_to_column("cluster") %>%
  mutate(
    total = rowSums(.),
    percent_control = ifelse("control" %in% colnames(.), control / total * 100, NA),
    percent_treatment = ifelse("treatment" %in% colnames(.), treatment / total * 100, NA)
  )

write.csv(cluster_stats, here(output_dir, "cluster_statistics.csv"), row.names = FALSE)

cat("数据降维和聚类完成，共鉴定到", length(unique(seurat_obj_corrected$seurat_clusters)), "个细胞簇\n")

# ----------------------
# 8. 保存最终预处理结果
# ----------------------

cat("\n保存最终预处理结果...\n")

# 8.1 保存最终的Seurat对象
saveRDS(seurat_obj_corrected, here(output_dir, "final_preprocessed_seurat_object.rds"))

# 8.2 创建预处理报告
preprocessing_report <- data.frame(
  step = c(
    "原始数据",
    "质量控制后",
    "去除doublet后",
    "最终预处理后"
  ),
  n_cells = c(
    ncol(seurat_obj),
    ncol(seurat_obj_filtered),
    ncol(seurat_obj_singlet),
    ncol(seurat_obj_corrected)
  ),
  n_genes = c(
    nrow(seurat_obj),
    nrow(seurat_obj_filtered),
    nrow(seurat_obj_singlet),
    nrow(seurat_obj_corrected)
  ),
  n_clusters = c(
    NA,
    NA,
    NA,
    length(unique(seurat_obj_corrected$seurat_clusters))
  )
) %>%
  mutate(
    cells_retained = n_cells / n_cells[1] * 100,
    genes_retained = n_genes / n_genes[1] * 100
  )

write.csv(preprocessing_report, here(output_dir, "preprocessing_report.csv"), row.names = FALSE)

# 打印预处理报告
cat("\n预处理报告:\n")
print(preprocessing_report)

# ----------------------
# 9. 生成分析报告
# ----------------------

cat("\n生成分析报告...\n")

# 创建分析报告
report_content <- glue(
  "# CD28-Treg-IBD-Study: 单细胞RNA测序数据预处理报告\n\n",
  
  "## 分析概述\n",
  "- **分析日期**: {Sys.Date()}\n",
  "- **数据类型**: 单细胞RNA测序数据\n",
  "- **分析步骤**: 数据加载、质量控制、doublet检测、归一化、批次效应校正、高可变基因鉴定、降维和聚类\n\n",
  
  "## 数据概况\n",
  "### 样本信息\n",
  "- **样本数量**: {length(unique(seurat_obj_corrected$sample))}\n",
  "- **批次数量**: {length(unique(seurat_obj_corrected$batch))}\n",
  "- **分组数量**: {length(unique(seurat_obj_corrected$group))}\n",
  "- **分组信息**: {paste(unique(seurat_obj_corrected$group), collapse = ', ')}\n\n",
  
  "### 预处理结果\n",
  "| 步骤 | 细胞数 | 基因数 | 保留比例(细胞) | 保留比例(基因) |\n",
  "|------|--------|--------|----------------|----------------|\n",
  "{paste(sprintf('| %s | %d | %d | %.1f%% | %.1f%% |', \n",
  "        preprocessing_report$step, \n",
  "        preprocessing_report$n_cells, \n",
  "        preprocessing_report$n_genes, \n",
  "        preprocessing_report$cells_retained, \n",
  "        preprocessing_report$genes_retained), collapse = '\\n')}\n\n",
  
  "## 质量控制\n",
  "### 质量控制阈值\n",
  "- **最小特征数**: {nFeature_RNA_min}\n",
  "- **最大特征数**: {nFeature_RNA_max}\n",
  "- **最小计数**: {nCount_RNA_min}\n",
  "- **最大线粒体比例**: {percent_mt_max}%\n\n",
  
  "### 质量控制结果\n",
  "- **过滤前细胞数**: {ncol(seurat_obj)}\n",
  "- **过滤后细胞数**: {ncol(seurat_obj_filtered)}\n",
  "- **过滤掉的细胞数**: {ncol(seurat_obj) - ncol(seurat_obj_filtered)} ({sprintf('%.1f', (ncol(seurat_obj) - ncol(seurat_obj_filtered)) / ncol(seurat_obj) * 100)}%)\n\n",
  
  "## Doublet检测\n",
  "- **scDblFinder检测到的doublet**: {sum(scdblfinder_doublets)} ({sprintf('%.1f', sum(scdblfinder_doublets)/ncol(seurat_obj_filtered)*100)}%)\n",
  "- **DoubletFinder检测到的doublet**: {sum(doubletfinder_doublets)} ({sprintf('%.1f', sum(doubletfinder_doublets)/ncol(seurat_obj_filtered)*100)}%)\n",
  "- **合并检测到的doublet**: {sum(seurat_obj_filtered$doublet)} ({sprintf('%.1f', sum(seurat_obj_filtered$doublet)/ncol(seurat_obj_filtered)*100)}%)\n",
  "- **去除doublet后的细胞数**: {ncol(seurat_obj_singlet)}\n\n",
  
  "## 数据归一化\n",
  "- **归一化方法**: SCTransform\n",
  "- **校正因素**: 线粒体比例、总计数\n\n",
  
  "## 批次效应校正\n",
  "{if (length(unique(seurat_obj_normalized$batch)) > 1) {\n",
  "  sprintf('- **批次数量**: %d\\n- **校正方法**: Harmony', length(unique(seurat_obj_normalized$batch)))\n",
  "} else {\n",
  "  '未检测到多个批次，未进行批次效应校正'\n",
  "}}\n\n",
  
  "## 高可变基因鉴定\n",
  "- **鉴定方法**: SCTransform\n",
  "- **高可变基因数量**: {length(hvg_sct)}\n",
  "- **与传统方法重叠率**: {sprintf('%.1f', length(common_hvg)/length(hvg_sct)*100)}%\n\n",
  
  "## 细胞聚类\n",
  "- **降维方法**: UMAP, t-SNE\n",
  "- **聚类方法**: 基于图形的聚类\n",
  "- **聚类数量**: {length(unique(seurat_obj_corrected$seurat_clusters))}\n\n",
  
  "### 聚类统计\n",
  "{cluster_summary <- table(seurat_obj_corrected$seurat_clusters) %>% as.data.frame()\n",
  "colnames(cluster_summary) <- c('cluster', 'n_cells')\n",
  "cluster_summary$percent <- cluster_summary$n_cells / sum(cluster_summary$n_cells) * 100\n",
  "paste(sprintf('| %s | %d | %.1f%% |', cluster_summary$cluster, cluster_summary$n_cells, cluster_summary$percent), collapse = '\\n')}\n\n",
  
  "## 输出文件\n",
  "- `raw_seurat_object.rds`: 原始Seurat对象\n",
  "- `filtered_seurat_object.rds`: 质量控制后的Seurat对象\n",
  "- `singlet_seurat_object.rds`: 去除doublet后的Seurat对象\n",
  "- `normalized_seurat_object.rds`: 归一化后的Seurat对象\n",
  "- `batch_corrected_seurat_object.rds`: 批次效应校正后的Seurat对象\n",
  "- `final_preprocessed_seurat_object.rds`: 最终预处理后的Seurat对象\n",
  "- `qc_metrics.csv`: 质量控制指标\n",
  "- `preprocessing_report.csv`: 预处理报告\n",
  "- `cluster_statistics.csv`: 聚类统计\n",
  "- `variable_features.csv`: 高可变基因列表\n",
  "- `qc_violin_plots.png`: 质量控制小提琴图\n",
  "- `qc_scatter_plots.png`: 质量控制散点图\n",
  "- `batch_correction_comparison.png`: 批次效应校正比较图\n",
  "- `variable_features_visualization.png`: 高可变基因可视化\n",
  "- `clustering_visualization.png`: 聚类结果可视化\n\n",
  
  "## 下一步分析建议\n",
  "- 细胞类型注释\n",
  "- 差异表达分析\n",
  "- 细胞亚群分析\n",
  "- 轨迹分析\n",
  "- 基因调控网络分析\n\n",
  
  "---\n",
  "*报告生成时间: {Sys.time()}*"
)

# 保存报告
writeLines(
  report_content,
  here(output_dir, "single_cell_preprocessing_report.md")
)

# ----------------------
# 10. 完成预处理
# ----------------------

cat("单细胞RNA测序数据预处理完成！\n")
cat("预处理结果已保存到:", output_dir, "\n")
cat("主要输出文件包括:\n")
cat("- final_preprocessed_seurat_object.rds\n")
cat("- preprocessing_report.csv\n")
cat("- clustering_visualization.png\n")
cat("- single_cell_preprocessing_report.md\n")

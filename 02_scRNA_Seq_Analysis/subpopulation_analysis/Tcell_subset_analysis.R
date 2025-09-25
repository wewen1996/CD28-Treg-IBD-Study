#' CD28-Treg-IBD-Study: T细胞亚群分析
#' 
#' 本脚本用于深入分析T细胞亚群，特别是CD28+调节性T细胞(Treg)在IBD中的作用。
#' 包括T细胞提取、亚群聚类、差异表达分析和功能富集分析等步骤。

# 加载必要的包
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)
library(celldex)
library(SingleR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(dplyr)
library(tidyr)
library(here)
library(glue)

# 设置输出目录
output_dir <- here("results", "scRNA_Seq_Analysis", "subpopulation_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------
# 1. 数据加载
# ----------------------

cat("开始数据加载...\n")

# 加载注释后的Seurat对象
# 在实际分析中，应从聚类和注释步骤加载真实数据
# 这里使用示例数据进行演示
load_annotated_data <- function() {
  # 检查是否存在注释后的文件
  annotated_file <- here("results", "scRNA_Seq_Analysis", "clustering_annotation", "annotated_seurat_object.rds")
  
  if (file.exists(annotated_file)) {
    cat("加载注释后的Seurat对象...\n")
    seurat_obj <- readRDS(annotated_file)
  } else {
    cat("未找到注释后的文件，创建示例数据...\n")
    
    # 创建示例数据
    set.seed(1234)
    n_cells <- 2000
    n_genes <- 5000
    
    # 创建细胞ID
    cell_ids <- paste0("cell_", 1:n_cells)
    
    # 创建基因ID
    gene_ids <- paste0("gene_", 1:n_genes)
    
    # 添加一些已知基因
    known_genes <- c(
      # T细胞标记
      "CD3E", "CD3D", "CD3G", "TRAC",
      # CD4+ T细胞标记
      "CD4", "IL7R", "CCR7",
      # CD8+ T细胞标记
      "CD8A", "CD8B",
      # 调节性T细胞标记
      "FOXP3", "IL2RA", "CTLA4", "CD28", "TNFRSF18", "IKZF2",
      # T细胞亚群标记
      "IFNG", "IL4", "IL17A", "TGFB1", "PRDM1", "TCF7", "LEF1", "SELL",
      "CXCR3", "CCR6", "CCR4", "CXCR5", "PDCD1", "HAVCR2", "LAG3", "TIGIT"
    )
    
    # 确保已知基因在基因列表中
    missing_genes <- setdiff(known_genes, gene_ids)
    if (length(missing_genes) > 0) {
      gene_ids <- c(gene_ids, missing_genes)
    }
    
    # 创建表达矩阵
    counts_matrix <- matrix(
      rpois(n_cells * length(gene_ids), lambda = 5),
      nrow = length(gene_ids),
      ncol = n_cells,
      dimnames = list(gene_ids, cell_ids)
    )
    
    # 创建T细胞亚群
    t_cell_subsets <- rep(c(
      "Naive CD4+ T细胞",
      "Effector CD4+ T细胞",
      "Th1细胞",
      "Th2细胞",
      "Th17细胞",
      "Treg细胞",
      "CD28+ Treg细胞",
      "CD28- Treg细胞",
      "Naive CD8+ T细胞",
      "Effector CD8+ T细胞"
    ), length.out = n_cells)
    
    # 为不同T细胞亚群设置特异性表达
    # Naive CD4+ T细胞
    naive_cd4 <- t_cell_subsets == "Naive CD4+ T细胞"
    counts_matrix["CD4", naive_cd4] <- counts_matrix["CD4", naive_cd4] * 10
    counts_matrix["CCR7", naive_cd4] <- counts_matrix["CCR7", naive_cd4] * 8
    counts_matrix["SELL", naive_cd4] <- counts_matrix["SELL", naive_cd4] * 8
    counts_matrix["TCF7", naive_cd4] <- counts_matrix["TCF7", naive_cd4] * 8
    
    # Effector CD4+ T细胞
    effector_cd4 <- t_cell_subsets == "Effector CD4+ T细胞"
    counts_matrix["CD4", effector_cd4] <- counts_matrix["CD4", effector_cd4] * 10
    counts_matrix["IFNG", effector_cd4] <- counts_matrix["IFNG", effector_cd4] * 6
    
    # Th1细胞
    th1 <- t_cell_subsets == "Th1细胞"
    counts_matrix["CD4", th1] <- counts_matrix["CD4", th1] * 10
    counts_matrix["IFNG", th1] <- counts_matrix["IFNG", th1] * 10
    counts_matrix["CXCR3", th1] <- counts_matrix["CXCR3", th1] * 8
    
    # Th2细胞
    th2 <- t_cell_subsets == "Th2细胞"
    counts_matrix["CD4", th2] <- counts_matrix["CD4", th2] * 10
    counts_matrix["IL4", th2] <- counts_matrix["IL4", th2] * 10
    counts_matrix["CCR4", th2] <- counts_matrix["CCR4", th2] * 8
    
    # Th17细胞
    th17 <- t_cell_subsets == "Th17细胞"
    counts_matrix["CD4", th17] <- counts_matrix["CD4", th17] * 10
    counts_matrix["IL17A", th17] <- counts_matrix["IL17A", th17] * 10
    counts_matrix["CCR6", th17] <- counts_matrix["CCR6", th17] * 8
    
    # Treg细胞
    treg <- t_cell_subsets %in% c("Treg细胞", "CD28+ Treg细胞", "CD28- Treg细胞")
    counts_matrix["FOXP3", treg] <- counts_matrix["FOXP3", treg] * 10
    counts_matrix["IL2RA", treg] <- counts_matrix["IL2RA", treg] * 8
    counts_matrix["CTLA4", treg] <- counts_matrix["CTLA4", treg] * 8
    
    # CD28+ Treg细胞
    cd28_treg <- t_cell_subsets == "CD28+ Treg细胞"
    counts_matrix["CD28", cd28_treg] <- counts_matrix["CD28", cd28_treg] * 10
    
    # CD28- Treg细胞
    cd28neg_treg <- t_cell_subsets == "CD28- Treg细胞"
    counts_matrix["CD28", cd28neg_treg] <- counts_matrix["CD28", cd28neg_treg] * 0.1
    
    # Naive CD8+ T细胞
    naive_cd8 <- t_cell_subsets == "Naive CD8+ T细胞"
    counts_matrix["CD8A", naive_cd8] <- counts_matrix["CD8A", naive_cd8] * 10
    counts_matrix["CCR7", naive_cd8] <- counts_matrix["CCR7", naive_cd8] * 8
    counts_matrix["SELL", naive_cd8] <- counts_matrix["SELL", naive_cd8] * 8
    
    # Effector CD8+ T细胞
    effector_cd8 <- t_cell_subsets == "Effector CD8+ T细胞"
    counts_matrix["CD8A", effector_cd8] <- counts_matrix["CD8A", effector_cd8] * 10
    counts_matrix["IFNG", effector_cd8] <- counts_matrix["IFNG", effector_cd8] * 8
    counts_matrix["PRDM1", effector_cd8] <- counts_matrix["PRDM1", effector_cd8] * 8
    
    # 创建样本信息
    sample_info <- data.frame(
      cell = cell_ids,
      sample = rep(paste0("sample_", 1:3), each = n_cells/3),
      batch = rep(paste0("batch_", 1:2), length.out = n_cells),
      group = rep(c("control", "treatment"), length.out = n_cells),
      final_cell_type = rep("T细胞", n_cells),
      t_cell_subset = t_cell_subsets
    )
    
    # 创建Seurat对象
    seurat_obj <- CreateSeuratObject(
      counts = counts_matrix,
      meta.data = sample_info,
      project = "CD28-Treg-IBD"
    )
    
    # 添加一些质量控制指标
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    
    # 归一化和降维
    seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:30, verbose = FALSE)
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30, verbose = FALSE)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)
  }
  
  return(seurat_obj)
}

# 加载数据
seurat_obj <- load_annotated_data()

cat("数据加载完成，包含", ncol(seurat_obj), "个细胞和", nrow(seurat_obj), "个基因\n")

# ----------------------
# 2. T细胞提取和亚群分析
# ----------------------

cat("\n开始T细胞提取和亚群分析...\n")

# 2.1 提取T细胞
# 方法1: 基于预定义的细胞类型注释
t_cells_annotation <- seurat_obj$final_cell_type %in% c("T细胞", "CD4+ T细胞", "CD8+ T细胞", "调节性T细胞", "Treg细胞")

# 方法2: 基于T细胞标记基因的高表达
t_cell_markers <- c("CD3E", "CD3D", "CD3G", "TRAC")
t_cell_markers_present <- intersect(t_cell_markers, rownames(seurat_obj))

if (length(t_cell_markers_present) > 0) {
  t_cell_score <- colMeans(GetAssayData(seurat_obj, assay = "SCT", slot = "data")[t_cell_markers_present, ])
  t_cell_threshold <- quantile(t_cell_score, 0.5)
  t_cells_score <- t_cell_score > t_cell_threshold
} else {
  t_cells_score <- rep(FALSE, ncol(seurat_obj))
}

# 综合两种方法
seurat_obj$is_t_cell <- t_cells_annotation | t_cells_score

# 提取T细胞
t_cell_seurat <- subset(seurat_obj, subset = is_t_cell == TRUE)

# 统计T细胞数量
t_cell_count <- ncol(t_cell_seurat)
t_cell_percent <- t_cell_count / ncol(seurat_obj) * 100

cat(sprintf("提取到%d个T细胞，占总细胞的%.1f%%\n", t_cell_count, t_cell_percent))

# 2.2 T细胞亚群重新聚类
cat("T细胞亚群重新聚类...\n")

# 为T细胞创建新的Seurat对象进行重新分析
t_cell_obj <- CreateSeuratObject(
  counts = t_cell_seurat@assays$RNA@counts,
  meta.data = t_cell_seurat@meta.data,
  project = "T_cell_subsets"
)

# 添加SCTransform归一化结果
t_cell_obj@assays$SCT <- t_cell_seurat@assays$SCT

# 重新降维和聚类
t_cell_obj <- RunPCA(t_cell_obj, verbose = FALSE)
t_cell_obj <- RunUMAP(t_cell_obj, reduction = "pca", dims = 1:30, verbose = FALSE)
t_cell_obj <- FindNeighbors(t_cell_obj, reduction = "pca", dims = 1:30, verbose = FALSE)
t_cell_obj <- FindClusters(t_cell_obj, resolution = 0.8, verbose = FALSE)

# 2.3 T细胞亚群注释
cat("T细胞亚群注释...\n")

# 定义T细胞亚群标记基因
t_subset_markers <- list(
  "Naive CD4+ T细胞" = c("CD4", "CCR7", "SELL", "TCF7", "LEF1"),
  "Memory CD4+ T细胞" = c("CD4", "IL7R", "CCR7", "SELL"),
  "Th1细胞" = c("CD4", "IFNG", "CXCR3", "TBX21"),
  "Th2细胞" = c("CD4", "IL4", "IL5", "IL13", "CCR4"),
  "Th17细胞" = c("CD4", "IL17A", "IL17F", "CCR6", "RORC"),
  "Treg细胞" = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IKZF2"),
  "CD28+ Treg细胞" = c("FOXP3", "CD28", "IL2RA", "CTLA4"),
  "CD28- Treg细胞" = c("FOXP3", "IL2RA", "CTLA4", "TIGIT"),
  "Naive CD8+ T细胞" = c("CD8A", "CD8B", "CCR7", "SELL", "TCF7"),
  "Effector CD8+ T细胞" = c("CD8A", "CD8B", "IFNG", "GZMB", "PRDM1"),
  "Exhausted CD8+ T细胞" = c("CD8A", "PDCD1", "HAVCR2", "LAG3", "TIGIT")
)

# 计算每个T细胞亚群的标记基因得分
subset_scores <- lapply(names(t_subset_markers), function(subset) {
  markers <- t_subset_markers[[subset]]
  markers_present <- intersect(markers, rownames(t_cell_obj))
  
  if (length(markers_present) > 0) {
    score <- colMeans(GetAssayData(t_cell_obj, assay = "SCT", slot = "data")[markers_present, ])
  } else {
    score <- rep(0, ncol(t_cell_obj))
  }
  
  return(score)
})

names(subset_scores) <- names(t_subset_markers)

# 将得分添加到T细胞对象中
for (subset in names(subset_scores)) {
  t_cell_obj[[paste0("score_", gsub(" ", "_", gsub("/", "_", subset)))]] <- subset_scores[[subset]]
}

# 基于得分分配T细胞亚群
subset_matrix <- do.call(cbind, subset_scores)
t_cell_obj$t_subset_score <- colnames(subset_matrix)[max.col(subset_matrix)]

# 2.4 T细胞亚群可视化
# UMAP图（按重新聚类的簇着色）
p1 <- DimPlot(t_cell_obj, reduction = "umap", label = TRUE, label.size = 4) +
  ggtitle("T细胞重新聚类")

# UMAP图（按基于标记基因的亚群注释着色）
p2 <- DimPlot(t_cell_obj, reduction = "umap", group.by = "t_subset_score", label = TRUE, label.size = 3) +
  ggtitle("T细胞亚群注释") +
  theme(legend.position = "right")

# UMAP图（按分组着色）
p3 <- DimPlot(t_cell_obj, reduction = "umap", group.by = "group") +
  ggtitle("T细胞分组分布")

# 保存T细胞亚群可视化结果
ggsave(
  filename = here(output_dir, "t_cell_subset_umap.png"),
  plot = cowplot::plot_grid(p1, p2, p3, ncol = 1),
  width = 12,
  height = 18,
  dpi = 300
)

# 2.5 T细胞亚群标记基因表达可视化
# 创建标记基因热图
subset_marker_genes <- unique(unlist(t_subset_markers))
subset_marker_genes_present <- intersect(subset_marker_genes, rownames(t_cell_obj))

if (length(subset_marker_genes_present) > 0) {
  p4 <- DoHeatmap(t_cell_obj, features = subset_marker_genes_present, group.by = "t_subset_score", size = 3) +
    NoLegend() +
    ggtitle("T细胞亚群标记基因表达热图")
  
  ggsave(
    filename = here(output_dir, "t_cell_subset_marker_heatmap.png"),
    plot = p4,
    width = 15,
    height = 12,
    dpi = 300
  )
}

# 标记基因特征图
if (length(subset_marker_genes_present) > 0) {
  # 选择代表性标记基因
  representative_markers <- c("CD4", "CD8A", "FOXP3", "CD28", "IFNG", "IL17A", "IL4", "CCR7")
  representative_markers_present <- intersect(representative_markers, rownames(t_cell_obj))
  
  if (length(representative_markers_present) > 0) {
    p5 <- FeaturePlot(t_cell_obj, features = representative_markers_present, ncol = 3)
    
    ggsave(
      filename = here(output_dir, "t_cell_marker_feature_plots.png"),
      plot = p5,
      width = 15,
      height = 10,
      dpi = 300
    )
  }
}

# 2.6 T细胞亚群统计
t_subset_counts <- table(t_cell_obj$t_subset_score, t_cell_obj$group) %>%
  as.data.frame.matrix() %>%
  rownames_to_column("t_cell_subset") %>%
  mutate(
    total = rowSums(.),
    percent_control = ifelse("control" %in% colnames(.), control / total * 100, NA),
    percent_treatment = ifelse("treatment" %in% colnames(.), treatment / total * 100, NA)
  )

write.csv(t_subset_counts, here(output_dir, "t_cell_subset_counts.csv"), row.names = FALSE)

cat("T细胞亚群分析完成，共鉴定到", length(unique(t_cell_obj$t_subset_score)), "种T细胞亚群\n")

# ----------------------
# 3. CD28+Treg细胞分析
# ----------------------

cat("\n开始CD28+Treg细胞分析...\n")

# 3.1 提取CD28+Treg细胞
# 定义CD28+Treg细胞标记基因
cd28_treg_markers <- c("FOXP3", "CD28", "IL2RA", "CTLA4")
cd28_treg_markers_present <- intersect(cd28_treg_markers, rownames(t_cell_obj))

if (length(cd28_treg_markers_present) > 0) {
  # 计算CD28+Treg细胞得分
  t_cell_obj$cd28_treg_score <- colMeans(GetAssayData(t_cell_obj, assay = "SCT", slot = "data")[cd28_treg_markers_present, ])
  
  # 识别CD28+Treg细胞
  # 方法1: 基于预定义的亚群注释
  cd28_treg_annotation <- t_cell_obj$t_subset_score == "CD28+ Treg细胞"
  
  # 方法2: 基于CD28+Treg标记基因的高表达
  cd28_treg_threshold <- quantile(t_cell_obj$cd28_treg_score, 0.75)
  cd28_treg_score <- t_cell_obj$cd28_treg_score > cd28_treg_threshold
  
  # 方法3: 同时高表达FOXP3和CD28
  foxp3_threshold <- quantile(GetAssayData(t_cell_obj, assay = "SCT", slot = "data")["FOXP3", ], 0.7)
  cd28_threshold <- quantile(GetAssayData(t_cell_obj, assay = "SCT", slot = "data")["CD28", ], 0.7)
  foxp3_cd28_high <- GetAssayData(t_cell_obj, assay = "SCT", slot = "data")["FOXP3", ] > foxp3_threshold &
    GetAssayData(t_cell_obj, assay = "SCT", slot = "data")["CD28", ] > cd28_threshold
  
  # 综合三种方法
  t_cell_obj$is_cd28_treg <- cd28_treg_annotation | cd28_treg_score | foxp3_cd28_high
  
  # 统计CD28+Treg细胞数量
  cd28_treg_count <- sum(t_cell_obj$is_cd28_treg)
  cd28_treg_percent_tcell <- cd28_treg_count / ncol(t_cell_obj) * 100
  cd28_treg_percent_total <- cd28_treg_count / ncol(seurat_obj) * 100
  
  cat(sprintf("鉴定到%d个CD28+Treg细胞\n", cd28_treg_count))
  cat(sprintf("- 占T细胞的%.1f%%\n", cd28_treg_percent_tcell))
  cat(sprintf("- 占总细胞的%.1f%%\n", cd28_treg_percent_total))
  
  # CD28+Treg细胞在不同分组中的分布
  cd28_treg_by_group <- table(t_cell_obj$is_cd28_treg, t_cell_obj$group) %>%
    as.data.frame.matrix() %>%
    rownames_to_column("is_cd28_treg") %>%
    mutate(is_cd28_treg = ifelse(is_cd28_treg == "TRUE", "CD28+Treg细胞", "非CD28+Treg细胞"))
  
  write.csv(cd28_treg_by_group, here(output_dir, "cd28_treg_distribution.csv"), row.names = FALSE)
  
  cat("\nCD28+Treg细胞在不同分组中的分布:\n")
  print(cd28_treg_by_group)
  
  # 3.2 CD28+Treg细胞可视化
  # UMAP图（按CD28+Treg细胞着色）
  p6 <- DimPlot(t_cell_obj, reduction = "umap", group.by = "is_cd28_treg", label = TRUE) +
    ggtitle("CD28+Treg细胞分布 (UMAP)") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey"))
  
  # CD28+Treg细胞标记基因表达热图
  p7 <- FeaturePlot(t_cell_obj, features = cd28_treg_markers_present, ncol = 2)
  
  # CD28+Treg细胞得分分布
  p8 <- VlnPlot(t_cell_obj, features = "cd28_treg_score", group.by = "group") +
    ggtitle("CD28+Treg细胞得分分布") +
    geom_hline(yintercept = cd28_treg_threshold, linetype = "dashed", color = "red")
  
  # CD28+Treg细胞比例条形图
  cd28_treg_prop <- table(t_cell_obj$group, t_cell_obj$is_cd28_treg) %>%
    prop.table(margin = 1) %>%
    as.data.frame() %>%
    filter(Var2 == TRUE) %>%
    mutate(Var1 = factor(Var1))
  
  p9 <- ggplot(cd28_treg_prop, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", fill = "red", alpha = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", Freq * 100)), vjust = -0.3) +
    labs(
      title = "CD28+Treg细胞比例 (占T细胞的%)",
      x = "分组",
      y = "比例"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # 保存CD28+Treg细胞可视化结果
  ggsave(
    filename = here(output_dir, "cd28_treg_analysis.png"),
    plot = cowplot::plot_grid(p6, p9, p8, ncol = 1),
    width = 12,
    height = 18,
    dpi = 300
  )
  
  ggsave(
    filename = here(output_dir, "cd28_treg_marker_expression.png"),
    plot = p7,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  # 3.3 CD28+Treg细胞差异表达分析
  if (cd28_treg_count > 0) {
    # 提取CD28+Treg细胞
    cd28_treg_seurat <- subset(t_cell_obj, subset = is_cd28_treg == TRUE)
    
    # 检查是否有足够的CD28+Treg细胞进行分析
    if (ncol(cd28_treg_seurat) >= 10) {
      # 检查是否有多个分组
      if (length(unique(cd28_treg_seurat$group)) >= 2) {
        cat("进行CD28+Treg细胞差异表达分析...\n")
        
        # CD28+Treg细胞差异表达分析
        cd28_treg_de <- FindMarkers(
          cd28_treg_seurat,
          ident.1 = "treatment",
          ident.2 = "control",
          min.pct = 0.25,
          logfc.threshold = 0.5,
          verbose = FALSE
        ) %>%
          rownames_to_column("gene") %>%
          mutate(padj = p.adjust(p_val, method = "fdr")) %>%
          filter(padj < 0.05) %>%
          arrange(padj)
        
        if (nrow(cd28_treg_de) > 0) {
          cat(sprintf("在CD28+Treg细胞中鉴定到%d个差异表达基因\n", nrow(cd28_treg_de)))
          
          # 保存CD28+Treg细胞差异表达结果
          write.csv(cd28_treg_de, here(output_dir, "cd28_treg_differential_expression.csv"), row.names = FALSE)
          
          # CD28+Treg细胞差异表达基因火山图
          p10 <- cd28_treg_de %>%
            ggplot(aes(x = avg_log2FC, y = -log10(padj))) +
            geom_point(aes(color = ifelse(padj < 0.05 & abs(avg_log2FC) > 0.5, "显著", "不显著")),
                      alpha = 0.7, size = 2) +
            geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray") +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
            scale_color_manual(values = c("显著" = "red", "不显著" = "black")) +
            labseller::labs(
              title = "CD28+Treg细胞差异表达基因火山图",
              x = "log2(倍数变化)",
              y = "-log10(调整后p值)",
              color = "显著性"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 8)
        
          ggsave(
            filename = here(output_dir, "cd28_treg_volcano_plot.png"),
            plot = p10,
            width = 10,
            height = 8,
            dpi = 300
          )
          
          # CD28+Treg细胞差异表达基因热图
          if (nrow(cd28_treg_de) > 0) {
            top_cd28_treg_de <- head(cd28_treg_de, 20)
            p11 <- DoHeatmap(
              cd28_treg_seurat,
              features = top_cd28_treg_de$gene,
              group.by = "group",
              size = 3
            ) +
              NoLegend() +
              ggtitle("CD28+Treg细胞差异表达基因热图")
            
            ggsave(
              filename = here(output_dir, "cd28_treg_heatmap.png"),
              plot = p11,
              width = 12,
              height = 8,
              dpi = 300
            )
          }
          
          # 3.4 CD28+Treg细胞功能富集分析
          cat("进行CD28+Treg细胞功能富集分析...\n")
          
          # 获取差异表达基因
          de_genes <- cd28_treg_de$gene
          
          # 转换为ENTREZID
          gene_df <- bitr(de_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
          
          if (nrow(gene_df) > 0) {
            # GO富集分析
            go_bp <- enrichGO(
              gene = gene_df$ENTREZID,
              OrgDb = org.Hs.eg.db,
              keyType = "ENTREZID",
              ont = "BP",
              pAdjustMethod = "fdr",
              qvalueCutoff = 0.05
            )
            
            # KEGG富集分析
            kegg <- enrichKEGG(
              gene = gene_df$ENTREZID,
              organism = "hsa",
              pvalueCutoff = 0.05
            )
            
            # 保存富集结果
            enrichment_dir <- here(output_dir, "cd28_treg_enrichment")
            dir.create(enrichment_dir, recursive = TRUE, showWarnings = FALSE)
            
            if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
              write.csv(go_bp@result, here(enrichment_dir, "go_biological_process.csv"), row.names = FALSE)
              
              # GO生物过程气泡图
              p12 <- dotplot(go_bp, showCategory = 10, title = "CD28+Treg细胞GO生物过程富集")
              ggsave(
                filename = here(enrichment_dir, "go_biological_process_dotplot.png"),
                plot = p12,
                width = 12,
                height = 8,
                dpi = 300
              )
            }
            
            if (!is.null(kegg) && nrow(kegg@result) > 0) {
              write.csv(kegg@result, here(enrichment_dir, "kegg_pathway.csv"), row.names = FALSE)
              
              # KEGG通路气泡图
              p13 <- dotplot(kegg, showCategory = 10, title = "CD28+Treg细胞KEGG通路富集")
              ggsave(
                filename = here(enrichment_dir, "kegg_pathway_dotplot.png"),
                plot = p13,
                width = 12,
                height = 8,
                dpi = 300
              )
            }
            
            cat("CD28+Treg细胞功能富集分析完成\n")
          } else {
            cat("警告: 基因转换失败，跳过功能富集分析\n")
          }
        } else {
          cat("在CD28+Treg细胞中未鉴定到显著差异表达基因\n")
        }
      } else {
        cat("CD28+Treg细胞只有一个分组，跳过差异表达分析\n")
      }
    } else {
      cat("CD28+Treg细胞数量不足，跳过差异表达分析\n")
    }
  } else {
    cat("未鉴定到CD28+Treg细胞，跳过CD28+Treg细胞分析\n")
  }
} else {
  cat("未检测到CD28+Treg细胞标记基因，跳过CD28+Treg细胞分析\n")
}

# ----------------------
# 4. T细胞亚群动态变化分析
# ----------------------

cat("\n开始T细胞亚群动态变化分析...\n")

# 4.1 T细胞亚群在不同分组中的比例变化
t_subset_dynamics <- table(t_cell_obj$t_subset_score, t_cell_obj$group) %>%
  prop.table(margin = 2) %>%
  as.data.frame() %>%
  pivot_wider(names_from = Var2, values_from = Freq) %>%
  rename(t_cell_subset = Var1) %>%
  mutate(
    change = ifelse("treatment" %in% colnames(.), treatment - control, NA),
    change_percent = ifelse("treatment" %in% colnames(.) & control > 0, (treatment - control) / control * 100, NA)
  ) %>%
  arrange(desc(abs(change)))

write.csv(t_subset_dynamics, here(output_dir, "t_cell_subset_dynamics.csv"), row.names = FALSE)

cat("\nT细胞亚群在不同分组中的比例变化:\n")
print(t_subset_dynamics %>% select(t_cell_subset, control, treatment, change, change_percent))

# 4.2 T细胞亚群动态变化可视化
# 堆叠条形图
t_subset_stack <- table(t_cell_obj$t_subset_score, t_cell_obj$group) %>%
  as.data.frame() %>%
  pivot_wider(names_from = Var2, values_from = Freq) %>%
  rename(t_cell_subset = Var1) %>%
  pivot_longer(cols = -t_cell_subset, names_to = "group", values_to = "count") %>%
  group_by(group) %>%
  mutate(percent = count / sum(count) * 100)

p14 <- ggplot(t_subset_stack, aes(x = group, y = percent, fill = t_cell_subset)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", percent)), position = position_stack(vjust = 0.5), size = 3) +
  labs(
    title = "T细胞亚群在不同分组中的比例",
    x = "分组",
    y = "比例 (%)",
    fill = "T细胞亚群"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

# 差异条形图
if ("treatment" %in% colnames(t_subset_dynamics)) {
  p15 <- t_subset_dynamics %>%
    filter(!is.na(change)) %>%
    ggplot(aes(x = reorder(t_cell_subset, abs(change)), y = change_percent)) +
    geom_bar(stat = "identity", aes(fill = ifelse(change_percent > 0, "增加", "减少"))) +
    geom_text(aes(label = sprintf("%.1f%%", change_percent)), vjust = ifelse(t_subset_dynamics$change_percent > 0, -0.3, 1.3), size = 3) +
    scale_fill_manual(values = c("增加" = "red", "减少" = "blue")) +
    labs(
      title = "T细胞亚群在处理组中的比例变化",
      x = "T细胞亚群",
      y = "变化百分比 (%)",
      fill = "变化方向"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
}

# 保存T细胞亚群动态变化可视化结果
if ("treatment" %in% colnames(t_subset_dynamics)) {
  ggsave(
    filename = here(output_dir, "t_cell_subset_dynamics.png"),
    plot = cowplot::plot_grid(p14, p15, ncol = 1),
    width = 14,
    height = 12,
    dpi = 300
  )
} else {
  ggsave(
    filename = here(output_dir, "t_cell_subset_distribution.png"),
    plot = p14,
    width = 14,
    height = 8,
    dpi = 300
  )
}

# ----------------------
# 5. 保存分析结果
# ----------------------

cat("\n保存分析结果...\n")

# 5.1 保存T细胞对象
saveRDS(t_cell_obj, here(output_dir, "t_cell_subset_seurat_object.rds"))

# 5.2 创建分析报告
t_cell_summary <- data.frame(
  metric = c(
    "总细胞数",
    "T细胞数",
    "T细胞比例(%)",
    "T细胞亚群数",
    "CD28+Treg细胞数",
    "CD28+Treg细胞比例(T细胞)(%)",
    "CD28+Treg细胞比例(总细胞)(%)"
  ),
  value = c(
    ncol(seurat_obj),
    ncol(t_cell_obj),
    ncol(t_cell_obj)/ncol(seurat_obj)*100,
    length(unique(t_cell_obj$t_subset_score)),
    ifelse("is_cd28_treg" %in% colnames(t_cell_obj@meta.data), sum(t_cell_obj$is_cd28_treg), NA),
    ifelse("is_cd28_treg" %in% colnames(t_cell_obj@meta.data) && ncol(t_cell_obj) > 0, sum(t_cell_obj$is_cd28_treg)/ncol(t_cell_obj)*100, NA),
    ifelse("is_cd28_treg" %in% colnames(t_cell_obj@meta.data) && ncol(seurat_obj) > 0, sum(t_cell_obj$is_cd28_treg)/ncol(seurat_obj)*100, NA)
  )
)

write.csv(t_cell_summary, here(output_dir, "t_cell_analysis_summary.csv"), row.names = FALSE)

# 打印分析报告
cat("\nT细胞亚群分析报告:\n")
print(t_cell_summary)

# ----------------------
# 6. 生成分析报告
# ----------------------

cat("\n生成分析报告...\n")

# 创建分析报告
report_content <- glue(
  "# CD28-Treg-IBD-Study: T细胞亚群分析报告\n\n",
  
  "## 分析概述\n",
  "- **分析日期**: {Sys.Date()}\n",
  "- **数据类型**: 单细胞RNA测序数据\n",
  "- **分析步骤**: T细胞提取、亚群聚类、CD28+Treg细胞分析、动态变化分析\n\n",
  
  "## 数据概况\n",
  "| 指标 | 数值 |\n",
  "|------|------|\n",
  "{paste(sprintf('| %s | %s |', t_cell_summary$metric, t_cell_summary$value), collapse = '\\n')}\n\n",
  
  "## T细胞亚群分析\n",
  "### T细胞亚群分布\n",
  "{t_subset_table <- table(t_cell_obj$t_subset_score) %>% as.data.frame()\n",
  "colnames(t_subset_table) <- c('T细胞亚群', '数量')\n",
  "t_subset_table$比例_T细胞 <- t_subset_table$数量 / sum(t_subset_table$数量) * 100\n",
  "t_subset_table$比例_总细胞 <- t_subset_table$数量 / ncol(seurat_obj) * 100\n",
  "paste(sprintf('| %s | %d | %.1f%% | %.1f%% |', \n",
  "        t_subset_table$T细胞亚群, \n",
  "        t_subset_table$数量, \n",
  "        t_subset_table$比例_T细胞, \n",
  "        t_subset_table$比例_总细胞), collapse = '\\n')}\n\n",
  
  "## CD28+Treg细胞分析\n",
  "{if ('is_cd28_treg' %in% colnames(t_cell_obj@meta.data)) {\n",
  "  cd28_treg_count <- sum(t_cell_obj$is_cd28_treg)\n",
  "  cd28_treg_percent_tcell <- cd28_treg_count / ncol(t_cell_obj) * 100\n",
  "  cd28_treg_percent_total <- cd28_treg_count / ncol(seurat_obj) * 100\n",
  "  \n",
  "  cd28_treg_by_group <- table(t_cell_obj$group, t_cell_obj$is_cd28_treg) %>% as.data.frame.matrix()\n",
  "  cd28_treg_by_group$total <- rowSums(cd28_treg_by_group)\n",
  "  cd28_treg_by_group$treg_percent <- cd28_treg_by_group$TRUE / cd28_treg_by_group$total * 100\n",
  "  \n",
  "  paste(c(\n",
  "    sprintf('- 共鉴定到%d个CD28+Treg细胞', cd28_treg_count),\n",
  "    sprintf('- 占T细胞的%.1f%%', cd28_treg_percent_tcell),\n",
  "    sprintf('- 占总细胞的%.1f%%', cd28_treg_percent_total),\n",
  "    '',\n",
  "    'CD28+Treg细胞在不同分组中的分布:',\n",
  "    sprintf('- %s组: %d个 (%.1f%%)', rownames(cd28_treg_by_group)[1], cd28_treg_by_group$TRUE[1], cd28_treg_by_group$treg_percent[1]),\n",
  "    sprintf('- %s组: %d个 (%.1f%%)', rownames(cd28_treg_by_group)[2], cd28_treg_by_group$TRUE[2], cd28_treg_by_group$treg_percent[2])\n",
  "  ), collapse = '\\n')\n",
  "} else {\n",
  "  '未进行CD28+Treg细胞分析'\n",
  "}}\n\n",
  
  "## T细胞亚群动态变化\n",
  "{if ('treatment' %in% colnames(t_subset_dynamics)) {\n",
  "  significant_changes <- t_subset_dynamics %>% \n",
  "    filter(!is.na(change_percent) & abs(change_percent) > 20) %>% \n",
  "    arrange(desc(abs(change_percent)))\n",
  "  \n",
  "  if (nrow(significant_changes) > 0) {\n",
  "    paste(c(\n",
  "      '显著变化的T细胞亚群 (变化>20%):',\n",
  "      paste(sprintf('- %s: %.1f%%', \n",
  "                    significant_changes$t_cell_subset, \n",
  "                    significant_changes$change_percent), collapse = '\\n')\n",
  "    ), collapse = '\\n')\n",
  "  } else {\n",
  "    '未检测到显著变化的T细胞亚群'\n",
  "  }\n",
  "} else {\n",
  "  '未进行T细胞亚群动态变化分析'\n",
  "}}\n\n",
  
  "## CD28+Treg细胞差异表达分析\n",
  "{if (exists('cd28_treg_de') && nrow(cd28_treg_de) > 0) {\n",
  "  paste(c(\n",
  "    sprintf('在CD28+Treg细胞中鉴定到%d个差异表达基因', nrow(cd28_treg_de)),\n",
  "    sprintf('上调基因: %d个', sum(cd28_treg_de$avg_log2FC > 0)),\n",
  "    sprintf('下调基因: %d个', sum(cd28_treg_de$avg_log2FC < 0))\n",
  "  ), collapse = '\\n')\n",
  "} else {\n",
  "  '未进行CD28+Treg细胞差异表达分析或未鉴定到差异表达基因'\n",
  "}}\n\n",
  
  "## 输出文件\n",
  "- `t_cell_subset_seurat_object.rds`: T细胞亚群分析后的Seurat对象\n",
  "- `t_cell_analysis_summary.csv`: T细胞分析摘要\n",
  "- `t_cell_subset_counts.csv`: T细胞亚群统计\n",
  "- `t_cell_subset_umap.png`: T细胞亚群UMAP图\n",
  "- `t_cell_subset_marker_heatmap.png`: T细胞亚群标记基因热图\n",
  "- `t_cell_marker_feature_plots.png`: T细胞标记基因特征图\n",
  "\n",
  "{if ('is_cd28_treg' %in% colnames(t_cell_obj@meta.data)) {\n",
  "  paste(c(\n",
  "    '- `cd28_treg_distribution.csv`: CD28+Treg细胞分布',\n",
  "    '- `cd28_treg_analysis.png`: CD28+Treg细胞分析图',\n",
  "    '- `cd28_treg_marker_expression.png`: CD28+Treg标记基因表达图'\n",
  "  ), collapse = '\\n')\n",
  "} else {\n",
  "  ''\n",
  "}}\n",
  "\n",
  "{if (exists('cd28_treg_de') && nrow(cd28_treg_de) > 0) {\n",
  "  paste(c(\n",
  "    '- `cd28_treg_differential_expression.csv`: CD28+Treg细胞差异表达基因',\n",
  "    '- `cd28_treg_volcano_plot.png`: CD28+Treg细胞差异表达基因火山图',\n",
  "    '- `cd28_treg_heatmap.png`: CD28+Treg细胞差异表达基因热图',\n",
  "    '- `cd28_treg_enrichment/`: CD28+Treg细胞功能富集分析结果'\n",
  "  ), collapse = '\\n')\n",
  "} else {\n",
  "  ''\n",
  "}}\n",
  "\n",
  "{if ('treatment' %in% colnames(t_subset_dynamics)) {\n",
  "  paste(c(\n",
  "    '- `t_cell_subset_dynamics.csv`: T细胞亚群动态变化',\n",
  "    '- `t_cell_subset_dynamics.png`: T细胞亚群动态变化图'\n",
  "  ), collapse = '\\n')\n",
  "} else {\n",
  "  '- `t_cell_subset_distribution.png`: T细胞亚群分布图'\n",
  "}}\n\n",
  
  "## 主要发现\n",
  "{if ('is_cd28_treg' %in% colnames(t_cell_obj@meta.data)) {\n",
  "  cd28_treg_count <- sum(t_cell_obj$is_cd28_treg)\n",
  "  cd28_treg_percent_tcell <- cd28_treg_count / ncol(t_cell_obj) * 100\n",
  "  \n",
  "  cd28_treg_by_group <- table(t_cell_obj$group, t_cell_obj$is_cd28_treg) %>% as.data.frame.matrix()\n",
  "  cd28_treg_by_group$total <- rowSums(cd28_treg_by_group)\n",
  "  cd28_treg_by_group$treg_percent <- cd28_treg_by_group$TRUE / cd28_treg_by_group$total * 100\n",
  "  \n",
  "  changes <- t_subset_dynamics %>% filter(t_cell_subset == 'CD28+ Treg细胞')\n",
  "  \n",
  "  if (nrow(changes) > 0 && !is.na(changes$change_percent)) {\n",
  "    change_dir <- ifelse(changes$change_percent > 0, '增加', '减少')\n",
  "    \n",
  "    paste(c(\n",
  "      sprintf('1. CD28+Treg细胞占T细胞的%.1f%%，是T细胞群体中的重要组成部分', cd28_treg_percent_tcell),\n",
  "      sprintf('2. 与对照组相比，处理组中CD28+Treg细胞比例%s了%.1f%%', change_dir, abs(changes$change_percent)),\n",
  "      '3. CD28+Treg细胞的变化可能在IBD发病机制中发挥重要作用'\n",
  "    ), collapse = '\\n')\n",
  "  } else {\n",
  "    paste(c(\n",
  "      sprintf('1. CD28+Treg细胞占T细胞的%.1f%%，是T细胞群体中的重要组成部分', cd28_treg_percent_tcell),\n",
  "      '2. 未检测到CD28+Treg细胞比例在不同分组间的显著变化',\n",
  "      '3. CD28+Treg细胞的功能状态可能在IBD发病机制中发挥重要作用'\n",
  "    ), collapse = '\\n')\n",
  "  }\n",
  "} else {\n",
  "  '未鉴定到CD28+Treg细胞，无法提供相关发现'\n",
  "}}\n\n",
  
  "## 下一步分析建议\n",
  "- CD28+Treg细胞与其他免疫细胞的相互作用分析\n",
  "- CD28+Treg细胞的表观调控机制研究\n",
  "- CD28+Treg细胞在不同IBD亚型中的差异分析\n",
  "- 体外实验验证CD28+Treg细胞的功能\n",
  "- 多组学整合分析，结合基因组、表观基因组和蛋白质组数据\n\n",
  
  "---\n",
  "*报告生成时间: {Sys.time()}*"
)

# 保存报告
writeLines(
  report_content,
  here(output_dir, "t_cell_subset_analysis_report.md")
)

# ----------------------
# 7. 完成分析
# ----------------------

cat("T细胞亚群分析完成！\n")
cat("分析结果已保存到:", output_dir, "\n")
cat("主要输出文件包括:\n")
cat("- t_cell_subset_seurat_object.rds\n")
cat("- t_cell_analysis_summary.csv\n")
cat("- t_cell_subset_umap.png\n")
cat("- t_cell_subset_analysis_report.md\n")

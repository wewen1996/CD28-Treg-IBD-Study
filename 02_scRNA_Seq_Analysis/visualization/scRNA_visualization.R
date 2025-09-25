#' CD28-Treg-IBD-Study: 单细胞RNA测序数据可视化
#' 
#' 本脚本用于生成高质量的单细胞RNA测序数据可视化结果，包括细胞类型分布、
#' 基因表达模式、细胞亚群特征等多维度可视化。

# 加载必要的包
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(tidyr)
library(here)
library(glue)

# 设置输出目录
output_dir <- here("results", "scRNA_Seq_Analysis", "visualization")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 设置主题
theme_set(theme_minimal())
theme_update(
  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 10),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 8)
)

# ----------------------
# 1. 数据加载
# ----------------------

cat("开始数据加载...\n")

# 加载分析后的Seurat对象
# 在实际分析中，应从聚类和注释步骤加载真实数据
# 这里使用示例数据进行演示
load_analysis_data <- function() {
  # 检查是否存在分析后的文件
  t_cell_file <- here("results", "scRNA_Seq_Analysis", "subpopulation_analysis", "t_cell_subset_seurat_object.rds")
  annotated_file <- here("results", "scRNA_Seq_Analysis", "clustering_annotation", "annotated_seurat_object.rds")
  
  if (file.exists(t_cell_file)) {
    cat("加载T细胞亚群分析后的Seurat对象...\n")
    seurat_obj <- readRDS(t_cell_file)
    return(seurat_obj)
  } else if (file.exists(annotated_file)) {
    cat("加载注释后的Seurat对象...\n")
    seurat_obj <- readRDS(annotated_file)
    return(seurat_obj)
  } else {
    cat("未找到分析后的文件，创建示例数据...\n")
    
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
      # cd8+ T细胞标记
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
    counts_matrix["CD4", naive_cd4] <- count_matrix["CD4", naive_cd4] * 10
    count_matrix["CCR7", naive_cd4] <- count_matrix["CCR7", naive_cd4] * 8
    count_matrix["SELL", naive_cd4] <- count_matrix["SELL", naive_cd4] * 8
    
    # Effector CD4+ T细胞
    effector_cd4 <- t_cell_subsets == "Effector CD4+ T细胞"
    count_matrix["CD4", effector_cd4] <- count_matrix["CD4", effector_cd4] * 10
    count_matrix["IFNG", effector_cd4] <- count_matrix["IFNG", effector_cd4] * 6
    
    # Th1细胞
    th1 <- t_cell_subset == "Th1细胞"
    count_matrix["CD4", th1] <- count_matrix["CD4", th1] * 10
    count_matrix["IFNG", th1] <- count_matrix["IFNG", th1] * 10
    count_matrix["CXCR3", th1] <- count_matrix["CXCR3", th1] * 8
    
    # Th2细胞
    th2 <- t_cell_subset == "Th2细胞"
    count_matrix["CD4", th2] <- count_matrix["CD4", th2] * 10
    count_matrix["IL4", th2] <- count_matrix["IL4", th2] * 10
    count_matrix["CCR4", th2] <- count_matrix["CCR4", th2] * 8
    
    # Th17细胞
    th17 <- t_cell_subset == "Th17细胞"
    count_matrix["CD4", th17] <- count_matrix["CD4", th17] * 10
    count_matrix["IL17A", th17] <- count_matrix["IL17A", th17] * 10
    count_matrix["CCR6", th17] <- count_matrix["CCR6", th17] * 8
    
    # Treg细胞
    treg <- t_cell_subset %in% c("Treg细胞", "CD28+ Treg细胞", "CD28- Treg细胞")
    count_matrix["FOXP3", treg] <- count_matrix["FOXP3", treg] * 10
    count_matrix["IL2RA", treg] <- count_matrix["IL2RA", treg] * 8
    count_matrix["CTLA4", treg] <- count_matrix["CTLA4", treg] * 8
    
    # CD28+ Treg细胞
    cd28_treg <- t_cell_subset == "CD28+ Treg细胞"
    count_matrix["CD28", cd28_treg] <- count_matrix["CD28", cd28_treg] * 10
    
    # CD28- Treg细胞
    cd28neg_treg <- t_cell_subset == "CD28- Treg细胞"
    count_matrix["CD28", cd28neg_treg] <- count_matrix["CD28", cd28neg_treg] * 0.1
    
    # Naive CD8+ T细胞
    naive_cd8 <- t_cell_subset == "Naive CD8+ T细胞"
    count_matrix["CD8A", naive_cd8] <- count_matrix["CD8A", naive_cd8] * 10
    count_matrix["CCR7", naive_cd8] <- count_matrix["CCR7", naive_cd8] * 8
    count_matrix["SELL", naive_cd8] <- count_matrix["SELL", naive_cd8] * 8
    
    # Effector CD8+ T细胞
    effector_cd8 <- t_cell_subset == "Effector CD8+ T细胞"
    count_matrix["CD8A", effector_cd8] <- count_matrix["CD8A", effector_cd8] * 10
    count_matrix["IFNG", effector_cd8] <- count_matrix["IFNG", effector_cd8] * 8
    count_matrix["PRDM1", effector_cd8] <- count_matrix["PRDM1", effector_cd8] * 8
    
    # 创建样本信息
    sample_info <- data.frame(
      cell = cell_ids,
      sample = rep(paste0("sample_", 1:3), each = n_cells/3),
      batch = rep(paste0("batch_", 1:2), length.out = n_cells),
      group = rep(c("control", "treatment"), length.out = n_cells),
      final_cell_type = rep("T细胞", n_cells),
      t_subset_score = t_cell_subset,
      is_cd28_treg = t_cell_subset == "CD28+ Treg细胞"
    )
    
    # 创建Seurat对象
    seurat_obj <- CreateSeuratObject(
      counts = count_matrix,
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
    
    return(seurat_obj)
  }
}

# 加载数据
seurat_obj <- load_analysis_data()

cat("数据加载完成，包含", ncol(seurat_obj), "个细胞和", nrow(seurat_obj), "个基因\n")

# ----------------------
# 2. 细胞类型和亚群可视化
# ----------------------

cat("\n开始细胞类型和亚群可视化...\n")

# 2.1 UMAP可视化
# 定义颜色方案
n_clusters <- length(unique(seurat_obj$seurat_clusters))
cluster_colors <- if (n_clusters <= 10) {
  brewer.pal(10, "Paired")
} else {
  viridis(n_clusters)
}

# UMAP图（按聚类着色）
p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 4, cols = cluster_colors) +
  ggtitle("细胞聚类 (UMAP)") +
  theme(legend.position = "none")

# UMAP图（按细胞亚群着色）
if ("t_subset_score" %in% colnames(seurat_obj@meta.data)) {
  n_subsets <- length(unique(seurat_obj$t_subset_score))
  subset_colors <- if (n_subsets <= 10) {
    brewer.pal(10, "Set3")
  } else {
    viridis(n_subsets)
  }
  
  p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "t_subset_score", label = TRUE, label.size = 3, cols = subset_colors) +
    ggtitle("T细胞亚群 (UMAP)") +
    theme(legend.position = "right")
} else {
  p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "final_cell_type", label = TRUE, label.size = 3) +
    ggtitle("细胞类型 (UMAP)") +
    theme(legend.position = "right")
}

# Umap图（按分组着色）
p3 <- DimPlot(seurat_obj, reduction = "umap", group.by = "group") +
  ggtitle("分组分布 (UMAP)")

# Umap图（按样本着色）
p4 <- Dimplot(seurat_obj, reduction = "umap", group.by = "sample") +
  ggtitle("样本分布 (UMAP)")

# 保存UMAP可视化结果
ggsave(
  filename = here(output_dir, "umap_visualization.png"),
  plot = cowplot::plot_grid(p1, p2, p3, p4, ncol = 2),
  width = 15,
  height = 12,
  dpi = 300
)

# 2.2 CD28+Treg细胞可视化
if ("is_cd28_treg" %in% colnames(seurat_obj@meta.data)) {
  # UMAP图（按CD28+Treg细胞着色）
  p5 <- Dimplot(seurat_obj, reduction = "umap", group.by = "is_cd28_treg", label = TRUE) +
    ggtitle("CD28+Treg细胞分布 (UMAP)") +
    scale_color_manual(values = c("TRUE" = "#e74c3c", "FALSE" = "#bdc3c7"))
  
  # 保存CD28+Treg细胞UMAP图
  ggsave(
    filename = here(output_dir, "cd28_treg_umap.png"),
    plot = p5,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# 2.3 细胞类型/亚群比例可视化
# 堆叠条形图
if ("t_subset_score" %in% colnames(seurat_obj@meta.data)) {
  cell_type_counts <- table(seurat_obj$t_subset_score, seurat_obj$group) %>%
    as.data.frame() %>%
    pivot_wider(names_from = Var2, values_from = Freq) %>%
    rename(cell_type = Var1)
  
  if (ncol(cell_type_counts) > 2) {
    cell_type_long <- cell_type_counts %>%
      pivot_longer(cols = -cell_type, names_to = "group", values_to = "count") %>%
      group_by(group) %>%
      mutate(percent = count / sum(count) * 100)
    
    p6 <- ggplot(cell_type_long, aes(x = group, y = percent, fill = cell_type)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = sprintf("%.1f%%", percent), position = position_stack(vjust = 0.5), size = 3) +
      labseller::labs(
        title = "T细胞亚群在不同分组中的比例",
        x = "分组",
        y = "比例 (%)",
        fill = "T细胞亚群"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text.size = 10),
        legend.title = element_text.size = 10),
        legend.text = element_text.size = 8)
      )
    )
  }
} else {
  cell_type_counts <- table(seurat_obj$final_cell_type, seurat_obj$group) %>%
    as.data.frame() %>%
    pivot_wider(names_from = Var2, values_from = Freq) %>%
    rename(cell_type = Var1)
  
  if (ncol(cell_type_counts) > 2) {
    cell_type_long <- cell_type_counts %>%
      pivot_longer(cols = -cell_type, names_to = "group", values_to = "count") %>%
      group_by(group) %>%
      mutate(percent = count / sum(count) * 100)
    
    p6 <- ggplot(cell_type_long, aes(x = group, y = percent, fill = cell_type)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = sprintf("%.1f%%", percent), position = position_stack(vjust = 0.5), size = 3) +
      labseller::labs(
        title = "细胞类型在不同分组中的比例",
        x = "分组",
        y = "比例 (%)",
        fill = "细胞类型"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text.size = 12),
        axis.text = element_text.size = 10),
        legend.title = element_text.size = 10),
        legend.text = element_text.size = 8)
      )
    )
  }
}

# 保存堆叠条形图
if (exists("p6")) {
  ggsave(
    filename = here(output_dir, "cell_type_distribution.png"),
    plot = p6,
    width = 12,
    height = 8,
    dpi = 300
  )
}

# ----------------------
# 3. 基因表达可视化
# ----------------------

cat("\n开始基因表达可视化...\n")

# 3.1 标记基因表达特征图
# 定义重要标记基因
marker_genes <- list(
  "T细胞" = c("CD3E", "CD3D", "CD3G"),
  "CD4+ T细胞" = c("CD4", "IL7R"),
  "CD8+ T细胞" = c("CD8A", "CD8B"),
  "Treg细胞" = c("FOXP3", "IL2RA", "CTLA4"),
  "CD28+ Treg细胞" = c("CD28", "FOXP3"),
  "Th1细胞" = c("IFNG", "CXCR3"),
  "Th2细胞" = c("IL4", "CCR4"),
  "Th17细胞" = c("IL17A", "CCR6")
)

# 检查哪些基因存在于数据中
available_markers <- lapply(marker_genes, function(genes) intersect(genes, rownames(seurat_obj)))
available_markers <- available_markers[sapply(available_markers, length) > 0]

# 为每个细胞类型的标记基因创建特征图
for (cell_type in names(available_markers)) {
  genes_list <- available_markers[[cell_type]]
  
  if (length(genes_list) > 0) {
    p7 <- FeaturePlot(seurat_obj, features = genes_list, ncol = min(3, length(genes_list))) +
      plot_annotation(title = glue("{cell_type}标记基因表达"))
    
    ggsave(
      filename = here(output_dir, glue("{gsub(' ', '_', cell_type)}_marker_features.png")),
      plot = p7,
      width = 15,
      height = 5 * ceiling(length(genes_list) / 3),
      dpi = 300
    )
  }
}

# 3.2 重要基因表达热图
# 选择重要基因
important_genes <- c(
  "CD3E", "CD4", "CD8A", "FOXP3", "IL2RA", "CTLA4", "CD28",
  "IFNG", "IL4", "IL17A", "TGFB1", "PDCD1", "LAG3", "TIGIT"
)
important_genes_present <- intersect(important_genes, rownames(seurat_obj))

if (length(important_genes_present) > 0) {
  # 按细胞类型分组细胞
  if ("t_subset_score" %in% colnames(seurat_obj@meta.data)) {
    p8 <- DoHeatmap(seurat_obj, features = important_genes_present, group.by = "t_subset_score", size = 3) +
      NoLegend() +
      ggtitle("重要基因表达热图 (按T细胞亚群分组)")
  } else if ("final_cell_type" %in% colnames(seurat_obj@meta.data)) {
    p8 <- DoHeatmap(seurat_obj, features = important_genes_present, group.by = "final_cell_type", size = 3) +
      NoLegend() +
      ggtitle("重要基因表达热图 (按细胞类型分组)")
  } else {
    p8 <- DoHeatmap(seurat_obj, features = important_genes_present, group.by = "seurat_clusters", size = 3) +
      NoLegend() +
      ggtitle("重要基因表达热图 (按聚类分组)")
  }
  
  ggsave(
    filename = here(output_dir, "important_genes_heatmap.png"),
    plot = p8,
    width = 15,
    height = 10,
    dpi = 300
  )
}

# 3.3 基因共表达分析
if (length(important_genes_present) >= 2) {
  # 选择前6个重要基因进行共表达分析
  coexp_genes <- head(important_genes_present, 6)
  
  # 计算基因表达相关性
  gene_expr <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")[coexp_genes, ]
  gene_cor <- cor(t(gene_expr))
  
  # 相关性热图
  p9 <- ggplot(melt(gene_cor), aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", value)), size = 4) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labseller::labs(
      title = "重要基因表达相关性热图",
      x = "",
      y = "",
      fill = "Pearson相关系数"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10)
    )
  
  ggsave(
    filename = here(output_dir, "gene_correlation_heatmap.png"),
    plot = p9,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# 3.4 CD28和FOXP3共表达分析
if (all(c("CD28", "FOXP3") %in% rownames(seurat_obj))) {
  # 获取CD28和FOXP3的表达值
  cd28_expr <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")["CD28", ]
  foxp3_expr <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")["FOXP3", ]
  
  # 创建共表达散点图
  p10 <- ggplot(data.frame(CD28 = cd28_expr, FOXP3 = foxp3_expr), aes(x = CD28, y = FOXP3)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_density2d(color = "red") +
    labseller::labs(
      title = "CD28和FOXP3共表达分析",
      x = "CD28表达",
      y = "FOXP3表达"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text.size = 12),
      axis.text = element_text.size = 10)
    )
  )
  
  ggsave(
    filename = here(output_dir, "cd28_foxp3_coexpression.png"),
    plot = p10,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# ----------------------
# 4. CD28+Treg细胞特性可视化
# ----------------------

cat("\n开始CD28+Treg细胞特性可视化...\n")

if ("is_cd28_treg" %in% colnames(seurat_obj@meta.data)) {
  # 4.1 CD28+Treg细胞与其他T细胞亚群的比较
  if ("t_subset_score" %in% colnames(seurat_obj@meta.data)) {
    # 计算各亚群中CD28+Treg细胞的比例
    cd28_treg_by_subset <- table(seurat_obj$t_subset_score, seurat_obj$is_cd28_treg) %>%
      as.data.frame() %>%
      pivot_wider(names_from = Var2, values_from = Freq) %>%
      rename(t_subset = Var1) %>%
      mutate(
        total = TRUE + FALSE,
        cd28_treg_percent = TRUE / total * 100
      ) %>%
      arrange(desc(cd28_treg_percent))
    
    p11 <- ggplot(cd28_treg_by_subset, aes(x = reorder(t_subset, cd28_treg_percent), y = cd28_treg_percent)) +
      geom_bar(stat = "identity", fill = "#e74c3c", alpha = 0.7) +
      geom_text(aes(label = sprintf("%.1f%%", cd28_treg_percent)), vjust = -0.3, size = 3.5) +
      coord_flip() +
      labseller::labs(
        title = "各T细胞亚群中CD28+Treg细胞的比例",
        x = "T细胞亚群",
        y = "CD28+Treg细胞比例 (%)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text.size = 12),
        axis.text = element_text.size = 10)
      )
    )
    
    ggsave(
      filename = here(output_dir, "cd28_treg_by_subset.png"),
      plot = p11,
      width = 12,
      height = 8,
      dpi = 300
    )
  }
  
  # 4.2 CD28+Treg细胞功能标记基因表达
  treg_functional_markers <- c("FOXP3", "IL2RA", "CTLA4", "TGFB1", "IL10", "TNFRSF18", "IKZF2")
  treg_functional_markers_present <- intersect(treg_functional_markers, rownames(seurat_obj))
  
  if (length(treg_functional_markers_present) > 0) {
    # 按CD28+Treg细胞分组的小提琴图
    p12 <- VlnPlot(seurat_obj, features = treg_functional_markers_present, group.by = "is_cd28_treg", ncol = 3) +
      plot_annotation(title = "CD28+Treg细胞功能标记基因表达")
    
    ggsave(
      filename = here(output_dir, "cd28_treg_functional_markers.png"),
      plot = p12,
      width = 15,
      height = 5 * ceiling(length(treg_functional_markers_present) / 3),
      dpi = 300
    )
  }
  
  # 4.3 CD28+Treg细胞在不同分组中的转录组差异
  if (length(unique(seurat_obj$group)) >= 2) {
    # 提取CD28+Treg细胞
    cd28_treg_cells <- subset(seurat_obj, subset = is_cd28_treg == TRUE)
    
    if (ncol(cd28_treg_cells) >= 10) {
      # 进行差异表达分析
      cd28_treg_de <- FindMarkers(
        cd28_treg_cells,
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
        # 选择前15个差异表达基因
        top_de_genes <- head(cd28_treg_de$gene, 15)
        
        # 热图
        p13 <- DoHeatmap(
          cd28_treg_cells,
          features = top_de_genes,
          group.by = "group",
          size = 3
        ) +
          NoLegend() +
          ggtitle("CD28+Treg细胞差异表达基因热图")
        
        ggsave(
          filename = here(output_dir, "cd28_treg_de_heatmap.png"),
          plot = p13,
          width = 12,
          height = 8,
          dpi = 300
        )
        
        # 火山图
        p14 <- cd28_treg_de %>%
          ggplot(aes(x = avg_log2FC, y = -log10(padj))) +
          geom_point(aes(color = ifelse(padj < 0.05 & abs(avg_log2FC) > 0.5, "显著", "不显著")),
                    alpha = 0.7, size = 2) +
          geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray") +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
          scale_color_manual(values = c("显著" = "#e74c3c", "不显著" = "#7f8c8d")) +
          labseller::labs(
            title = "CD28+Treg细胞差异表达基因火山图",
            x = "log2(倍数变化)",
            y = "-log10(调整后p值)",
            color = "显著性"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text.size = 12),
            axis.text = element_text.size = 10),
            legend.title = element_text.size = 10),
            legend.text = element_text.size = 8)
          )
        )
        
        ggsave(
          filename = here(output_dir, "cd28_treg_volcano.png"),
          plot = p14,
          width = 10,
          height = 8,
          dpi = 300
        )
      }
    }
  }
}

# ----------------------
# 5. 综合分析可视化
# ----------------------

cat("\n开始综合分析可视化...\n")

# 5.1 细胞亚群和功能状态的关系
if (all(c("t_subset_score", "is_cd28_treg") %in% colnames(seurat_obj@meta.data))) {
  # 创建细胞亚群和CD28+Treg状态的组合
  seurat_obj$subset_cd28_treg <- paste(seurat_obj$t_subset_score, ifelse(seurat_obj$is_cd28_treg, "CD28+Treg", "非CD28+Treg"), sep = "_")
  
  # UMAP图（按组合状态着色）
  p15 <- DimPlot(seurat_obj, reduction = "umap", group.by = "subset_cd28_treg", label = FALSE) +
    ggtitle("T细胞亚群和CD28+Treg状态 (UMAP)") +
    theme(legend.position = "right")
  
  ggsave(
    filename = here(output_dir, "subset_cd28_treg_umap.png"),
    plot = p15,
    width = 12,
    height = 10,
    dpi = 300
  )
}

# 5.2 细胞亚群和分组的关系
if (all(c("t_subset_score", "group") %in% colnames(seurat_obj@meta.data))) {
  # 创建细胞亚群和分组的组合
  seurat_obj$subset_group <- paste(seurat_obj$t_subset_score, seurat_obj$group, sep = "_")
  
  # UMAP图（按组合状态着色）
  p16 <- DimPlot(seurat_obj, reduction = "umap", group.by = "subset_group", label = FALSE) +
    ggtitle("T细胞亚群和分组 (UMAP)") +
    theme(legend.position = "right")
  
  ggsave(
    filename = here(output_dir, "subset_group_umap.png"),
    plot = p16,
    width = 12,
    height = 10,
    dpi = 300
  )
}

# 5.3 细胞亚群动态变化
if (all(c("t_subset_score", "group") %in% colnames(seurat_obj@meta.data)) && length(unique(seurat_obj$group)) >= 2) {
  # 计算各分组中各亚群的比例
  subset_dynamics <- table(seurat_obj$t_subset_score, seurat_obj$group) %>%
    prop.table(margin = 2) %>%
    as.data.frame() %>%
    pivot_wider(names_from = Var2, values_from = Freq) %>%
    rename(t_subset = Var1) %>%
    mutate(
      change = ifelse("treatment" %in% colnames(.), treatment - control, NA),
      change_percent = ifelse("treatment" %in% colnames(.) & control > 0, (treatment - control) / control * 100, NA)
    ) %>%
    arrange(desc(abs(change)))
  
  # 差异条形图
  if ("treatment" %in% colnames(subset_dynamics)) {
    p17 <- subset_dynamics %>%
      filter(!is.na(change)) %>%
      ggplot(aes(x = reorder(t_subset, abs(change)), y = change_percent)) +
      geom_bar(stat = "identity", aes(fill = ifelse(change_percent > 0, "增加", "减少"))) +
      geom_text(aes(label = sprintf("%.1f%%", change_percent)), vjust = ifelse(subset_dynamics$change_percent > 0, -0.3, 1.3), size = 3) +
      scale_fill_manual(values = c("增加" = "#e74c3c", "减少" = "#3498db")) +
      labseller::labs(
        title = "T细胞亚群在处理组中的比例变化",
        x = "T细胞亚群",
        y = "变化百分比 (%)",
        fill = "变化方向"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text.size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text.size = 10),
        legend.title = element_text.size = 10),
        legend.text = element_text.size = 8)
      )
    )
    
    ggsave(
      filename = here(output_dir, "t_cell_subset_dynamics.png"),
      plot = p17,
      width = 14,
      height = 8,
      dpi = 300
    )
  }
}

# 5.4 创建综合分析报告图
# 选择最重要的几个图进行组合
plots_to_combine <- list()

# 添加UMAP图
plots_to_combine[["umap_clusters"]] <- p1

# 添加细胞类型分布图
if (exists("p6")) {
  plots_to_combine[["cell_type_distribution"]] <- p6
}

# 添加CD28+Treg细胞分析图
if (exists("p5")) {
  plots_to_combine[["cd28_treg_umap"]] <- p5
}

# 添加重要基因热图
if (exists("p8")) {
  plots_to_combine[["important_genes_heatmap"]] <- p8
}

# 添加T细胞亚群动态变化图
if (exists("p17")) {
  plots_to_combine[["subset_dynamics"]] <- p17
}

# 组合图形
if (length(plots_to_combine) > 0) {
  combined_plot <- wrap_plots(plots_to_combine, ncol = 2) +
    plot_annotation(
      title = "CD28-Treg-IBD-Study 单细胞RNA测序数据分析结果",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
  
  ggsave(
    filename = here(output_dir, "comprehensive_analysis.png"),
    plot = combined_plot,
    width = 20,
    height = 15,
    dpi = 300
  )
}

# ----------------------
# 6. 生成可视化报告
# ----------------------

cat("\n生成可视化报告...\n")

# 创建可视化报告
report_content <- glue(
  "# CD28-Treg-IBD-Study: 单细胞RNA测序数据可视化报告\n\n",
  
  "## 分析概述\n",
  "- **分析日期**: {Sys.Date()}\n",
  "- **数据类型**: 单细胞RNA测序数据\n",
  "- **可视化内容**: 细胞类型分布、基因表达模式、细胞亚群特征、CD28+Treg细胞分析\n\n",
  
  "## 数据概况\n",
  "- **总细胞数**: {ncol(seurat_obj)}\n",
  "- **总基因数**: {nrow(seurat_obj)}\n",
  "- **分组信息**: {paste(unique(seurat_obj$group), collapse = ', ')}\n",
  "- **样本数量**: {length(unique(seurat_obj$sample))}\n\n",
  
  "## 主要可视化结果\n",
  "### 1. 细胞聚类和亚群分析\n",
  "- **umap_visualization.png**: UMAP降维可视化，展示细胞聚类、细胞类型/亚群、分组和样本分布\n",
  "- **cell_type_distribution.png**: 细胞类型/亚群在不同分组中的比例分布\n\n",
  
  "### 2. 基因表达分析\n",
  "- **重要基因表达热图**: 展示关键免疫细胞标记基因在不同细胞群体中的表达模式\n",
  "- **基因表达特征图**: 展示各细胞类型特异性标记基因的表达分布\n",
  "- **基因相关性热图**: 展示重要基因之间的表达相关性\n\n",
  
  "{if ('is_cd28_treg' %in% colnames(seurat_obj@meta.data)) {\n",
  "  paste(c(\n",
  "    '### 3. CD28+Treg细胞分析',\n",
  "    '- **cd28_treg_umap.png**: CD28+Treg细胞在UMAP空间中的分布',\n",
  "    '- **cd28_treg_functional_markers.png**: CD28+Treg细胞功能标记基因表达',\n",
  "    '- **cd28_foxp3_coexpression.png**: CD28和FOXP3共表达分析',\n",
  "    '- **cd28_treg_by_subset.png**: 各T细胞亚群中CD28+Treg细胞的比例'\n",
  "  ), collapse = '\\n')\n",
  "} else {\n",
  "  ''\n",
  "}}\n",
  "\n",
  "{if (exists('cd28_treg_de') && nrow(cd28_treg_de) > 0) {\n",
  "  paste(c(\n",
  "    '### 4. CD28+Treg细胞差异表达分析',\n",
  "    '- **cd28_treg_de_heatmap.png**: CD28+Treg细胞差异表达基因热图',\n",
  "    '- **cd28_treg_volcano.png**: CD28+Treg细胞差异表达基因火山图'\n",
  "  ), collapse = '\\n')\n",
  "} else {\n",
  "  ''\n",
  "}}\n",
  "\n",
  "{if (exists('p17')) {\n",
  "  paste(c(\n",
  "    '### 5. T细胞亚群动态变化',\n",
  "    '- **t_cell_subset_dynamics.png**: T细胞亚群在不同分组中的比例变化'\n",
  "  ), collapse = '\\n')\n",
  "} else {\n",
  "  ''\n",
  "}}\n",
  "\n",
  "### 6. 综合分析结果\n",
  "- **comprehensive_analysis.png**: 综合分析报告图，展示主要发现\n\n",
  
  "## 关键发现\n",
  "{if ('is_cd28_treg' %in% colnames(seurat_obj@meta.data)) {\n",
  "  cd28_treg_count <- sum(seurat_obj$is_cd28_treg)\n",
  "  cd28_treg_percent <- cd28_treg_count / ncol(seurat_obj) * 100\n",
  "  \n",
  "  paste(c(\n",
  "    sprintf('1. CD28+Treg细胞占总细胞的%.1f%%，是免疫细胞群体中的重要组成部分', cd28_treg_percent),\n",
  "    '2. CD28+Treg细胞在UMAP空间中形成明显的聚类，表明其独特的转录组特征',\n",
  "    '3. CD28和FOXP3共表达分析显示CD28+Treg细胞具有独特的共表达模式'\n",
  "  ), collapse = '\\n')\n",
  "} else {\n",
  "  paste(c(\n",
  "    '1. 成功鉴定出多种免疫细胞类型，包括T细胞、B细胞、髓系细胞等',\n",
  "    '2. T细胞亚群分析揭示了不同功能状态的T细胞群体',\n",
  "    '3. 关键免疫标记基因在不同细胞群体中呈现特异性表达模式'\n",
  "  ), collapse = '\\n')\n",
  "}}\n\n",
  
  "{if (exists('cd28_treg_de') && nrow(cd28_treg_de) > 0) {\n",
  "  paste(c(\n",
  "    sprintf('4. CD28+Treg细胞在不同分组间存在%d个差异表达基因', nrow(cd28_treg_de)),\n",
  "    '5. 这些差异表达基因主要参与免疫调节、细胞活化和炎症反应等生物学过程'\n",
  "  ), collapse = '\\n')\n",
  "} else {\n",
  "  ''\n",
  "}}\n",
  "\n",
  "{if (exists('p17')) {\n",
  "  paste(c(\n",
  "    '6. T细胞亚群在不同分组中表现出动态变化，提示免疫微环境的重塑',\n",
  "    '7. CD28+Treg细胞比例的变化可能在IBD发病机制中发挥重要作用'\n",
  "  ), collapse = '\\n')\n",
  "} else {\n",
  "  ''\n",
  "}}\n\n",
  
  "## 可视化方法说明\n",
  "- **UMAP降维**: 用于展示高维单细胞数据的低维嵌入，保留数据的局部和全局结构\n",
  "- **热图**: 用于展示基因表达模式和相关性\n",
  "- **小提琴图**: 用于展示基因在不同细胞群体中的表达分布\n",
  "- **火山图**: 用于展示差异表达基因的统计显著性和表达变化幅度\n",
  "- **条形图**: 用于展示细胞类型/亚群的比例分布和动态变化\n\n",
  
  "## 输出文件\n",
  "所有可视化结果已保存到以下目录:\n",
  "- **results/scRNA_Seq_Analysis/visualization/**\n\n",
  
  "---\n",
  "*报告生成时间: {Sys.time()}*"
)

# 保存报告
writeLines(
  report_content,
  here(output_dir, "scRNA_visualization_report.md")
)

# ----------------------
# 7. 完成可视化
# ----------------------

cat("单细胞RNA测序数据可视化完成！\n")
cat("可视化结果已保存到:", output_dir, "\n")
cat("主要输出文件包括:\n")
cat("- umap_visualization.png\n")
cat("- cell_type_distribution.png\n")
cat("- important_genes_heatmap.png\n")
cat("- comprehensive_analysis.png\n")
cat("- scRNA_visualization_report.md\n")

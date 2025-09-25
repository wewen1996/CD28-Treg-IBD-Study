#' CD28-Treg-IBD-Study: 单细胞RNA测序数据分析主脚本
#' 
#' 本脚本用于单细胞RNA测序数据的聚类和注释，包括细胞类型识别、
#' 差异表达分析和功能富集分析等步骤。

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
output_dir <- here("results", "scRNA_Seq_Analysis", "clustering_annotation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------
# 1. 数据加载
# ----------------------

cat("开始数据加载...\n")

# 加载预处理后的Seurat对象
# 在实际分析中，应从预处理步骤加载真实数据
# 这里使用示例数据进行演示
load_preprocessed_data <- function() {
  # 检查是否存在预处理后的文件
  preprocessed_file <- here("results", "scRNA_Seq_Analysis", "preprocessing", "final_preprocessed_seurat_object.rds")
  
  if (file.exists(preprocessed_file)) {
    cat("加载预处理后的Seurat对象...\n")
    seurat_obj <- readRDS(preprocessed_file)
  } else {
    cat("未找到预处理后的文件，创建示例数据...\n")
    
    # 创建示例数据
    set.seed(1234)
    n_cells <- 2000
    n_genes <- 5000
    n_samples <- 3
    
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
    counts_matrix <- matrix(
      rpois(n_cells * n_samples * n_genes, lambda = 5),
      nrow = n_genes,
      ncol = n_cells * n_samples,
      dimnames = list(gene_ids, cell_ids)
    )
    
    # 添加一些已知基因
    known_genes <- c(
      # T细胞标记
      "CD3E", "CD3D", "CD4", "CD8A", "FOXP3", "IL2RA", "CTLA4", "CD28",
      # B细胞标记
      "CD19", "MS4A1", "CD79A", "IGHM",
      # 髓系细胞标记
      "CD14", "CD16", "CD68", "CD1C", "CLEC9A",
      # 内皮细胞标记
      "PECAM1", "VWF",
      # 上皮细胞标记
      "EPCAM", "KRT8", "KRT18"
    )
    
    # 确保已知基因在基因列表中
    missing_genes <- setdiff(known_genes, gene_ids)
    if (length(missing_genes) > 0) {
      gene_ids <- c(gene_ids, missing_genes)
      counts_matrix <- rbind(
        counts_matrix,
        matrix(rpois(n_cells * n_samples * length(missing_genes), lambda = 5),
               nrow = length(missing_genes),
               ncol = n_cells * n_samples,
               dimnames = list(missing_genes, cell_ids))
      )
    }
    
    # 为不同细胞类型设置特异性表达
    cell_types <- rep(c("T细胞", "B细胞", "髓系细胞", "内皮细胞", "上皮细胞"), 
                     length.out = n_cells * n_samples)
    sample_info$cell_type <- cell_types
    
    # T细胞高表达T细胞标记
    t_cells <- cell_types == "T细胞"
    counts_matrix["CD3E", t_cells] <- counts_matrix["CD3E", t_cells] * 10
    counts_matrix["CD3D", t_cells] <- counts_matrix["CD3D", t_cells] * 10
    counts_matrix["CD4", t_cells] <- counts_matrix["CD4", t_cells] * 8
    
    # B细胞高表达B细胞标记
    b_cells <- cell_types == "B细胞"
    counts_matrix["CD19", b_cells] <- counts_matrix["CD19", b_cells] * 10
    counts_matrix["MS4A1", b_cells] <- counts_matrix["MS4A1", b_cells] * 10
    
    # 髓系细胞高表达髓系细胞标记
    myeloid_cells <- cell_types == "髓系细胞"
    counts_matrix["CD14", myeloid_cells] <- counts_matrix["CD14", myeloid_cells] * 10
    counts_matrix["CD68", myeloid_cells] <- count_matrix["CD68", myeloid_cells] * 8
    
    # 内皮细胞高表达内皮细胞标记
    endothelial_cells <- cell_types == "内皮细胞"
    counts_matrix["PECAM1", endothelial_cells] <- count_matrix["PECAM1", endothelial_cells] * 10
    counts_matrix["VWF", endothelial_cells] <- count_matrix["VWF", endothelial_cells] * 8
    
    # 上皮细胞高表达上皮细胞标记
    epithelial_cells <- cell_types == "上皮细胞"
    counts_matrix["EPCAM", epithelial_cells] <- count_matrix["EPCAM", epithelial_cells] * 10
    counts_matrix["KRT8", epithelial_cells] <- count_matrix["KRT8", epithelial_cells] * 8
    
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
  }
  
  return(seurat_obj)
}

# 加载数据
seurat_obj <- load_preprocessed_data()

cat("数据加载完成，包含", ncol(seurat_obj), "个细胞和", nrow(seurat_obj), "个基因\n")

# ----------------------
# 2. 细胞类型注释
# ----------------------

cat("\n开始细胞类型注释...\n")

# 2.1 加载参考数据集
cat("加载参考数据集...\n")
ref_data <- celldex::HumanPrimaryCellAtlasData()

# 2.2 使用SingleR进行细胞类型注释
cat("使用SingleR进行细胞类型注释...\n")
set.seed(1234)
singleR_pred <- SingleR(
  test = GetAssayData(seurat_obj, assay = "SCT", slot = "data"),
  ref = ref_data,
  labels = ref_data$label.main,
  de.method = "wilcox"
)

# 将注释结果添加到Seurat对象中
seurat_obj$singleR_main <- singleR_pred$labels
seurat_obj$singleR_pruned <- singleR_pred$pruned.labels

# 2.3 基于已知标记基因的细胞类型注释
cat("基于已知标记基因的细胞类型注释...\n")

# 定义细胞类型标记基因
cell_type_markers <- list(
  "T细胞" = c("CD3E", "CD3D", "CD3G", "TRAC"),
  "CD4+ T细胞" = c("CD4", "IL7R", "CCR7"),
  "CD8+ T细胞" = c("CD8A", "CD8B"),
  "调节性T细胞" = c("FOXP3", "IL2RA", "CTLA4", "CD28"),
  "B细胞" = c("CD19", "MS4A1", "CD79A", "IGHM"),
  "浆细胞" = c("SDC1", "IGHG1", "MZB1"),
  "单核细胞" = c("CD14", "LYZ", "CD68"),
  "巨噬细胞" = c("CD68", "CD163", "MRC1"),
  "树突状细胞" = c("CD1C", "CLEC9A", "CD207"),
  "NK细胞" = c("NKG7", "KLRD1", "CD16"),
  "内皮细胞" = c("PECAM1", "VWF", "CD34"),
  "上皮细胞" = c("EPCAM", "KRT8", "KRT18"),
  "成纤维细胞" = c("ACTA2", "COL1A1", "FAP")
)

# 计算每个细胞类型的标记基因得分
cell_type_scores <- lapply(names(cell_type_markers), function(cell_type) {
  markers <- cell_type_markers[[cell_type]]
  markers_present <- intersect(markers, rownames(seurat_obj))
  
  if (length(markers_present) > 0) {
    score <- colMeans(GetAssayData(seurat_obj, assay = "SCT", slot = "data")[markers_present, ])
  } else {
    score <- rep(0, ncol(seurat_obj))
  }
  
  return(score)
})

names(cell_type_scores) <- names(cell_type_markers)

# 将得分添加到Seurat对象中
for (cell_type in names(cell_type_scores)) {
  seurat_obj[[paste0("score_", gsub(" ", "_", cell_type))]] <- cell_type_scores[[cell_type]]
}

# 基于得分分配细胞类型
cell_type_matrix <- do.call(cbind, cell_type_scores)
seurat_obj$marker_based_type <- colnames(cell_type_matrix)[max.col(cell_type_matrix)]

# 2.4 整合注释结果
cat("整合注释结果...\n")

# 创建最终的细胞类型注释
# 优先使用SingleR的修剪标签，如果没有则使用基于标记基因的注释
seurat_obj$final_cell_type <- ifelse(
  !is.na(seurat_obj$singleR_pruned),
  seurat_obj$singleR_pruned,
  seurat_obj$marker_based_type
)

# 2.5 细胞类型注释可视化
# UMAP图（按SingleR注释着色）
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "singleR_main", label = TRUE, label.size = 4) +
  ggtitle("SingleR细胞类型注释") +
  theme(legend.position = "right")

# Umap图（按基于标记基因的注释着色）
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "marker_based_type", label = TRUE, label.size = 4) +
  ggtitle("基于标记基因的细胞类型注释") +
  theme(legend.position = "right")

# Umap图（按最终注释着色）
p3 <- DimPlot(seurat_obj, reduction = "umap", group.by = "final_cell_type", label = TRUE, label.size = 4) +
  ggtitle("最终细胞类型注释") +
  theme(legend.position = "right")

# 保存细胞类型注释可视化结果
ggsave(
  filename = here(output_dir, "cell_type_annotation_umap.png"),
  plot = cowplot::plot_grid(p1, p2, p3, ncol = 1),
  width = 12,
  height = 18,
  dpi = 300
)

# 2.6 标记基因表达可视化
# 创建标记基因热图
marker_genes <- unique(unlist(cell_type_markers))
marker_genes_present <- intersect(marker_genes, rownames(seurat_obj))

if (length(marker_genes_present) > 0) {
  p4 <- DoHeatmap(seurat_obj, features = marker_genes_present, group.by = "final_cell_type", size = 3) +
    NoLegend() +
    ggtitle("细胞类型标记基因表达热图")
  
  ggsave(
    filename = here(output_dir, "marker_genes_heatmap.png"),
    plot = p4,
    width = 15,
    height = 10,
    dpi = 300
  )
}

# 2.7 细胞类型统计
cell_type_counts <- table(seurat_obj$final_cell_type, seurat_obj$group) %>%
  as.data.frame.matrix() %>%
  rownames_to_column("cell_type") %>%
  mutate(
    total = rowSums(.),
    percent_control = ifelse("control" %in% colnames(.), control / total * 100, NA),
    percent_treatment = ifelse("treatment" %in% colnames(.), treatment / total * 100, NA)
  )

write.csv(cell_type_counts, here(output_dir, "cell_type_counts.csv"), row.names = FALSE)

cat("细胞类型注释完成，共鉴定到", length(unique(seurat_obj$final_cell_type)), "种细胞类型\n")

# ----------------------
# 3. 差异表达分析
# ----------------------

cat("\n开始差异表达分析...\n")

# 3.1 按细胞类型进行差异表达分析
cell_types <- unique(seurat_obj$final_cell_type)
de_results <- list()

for (cell_type in cell_types) {
  cat(sprintf("分析%s的差异表达基因...\n", cell_type))
  
  # 提取该细胞类型的细胞
  cell_type_cells <- subset(seurat_obj, subset = final_cell_type == cell_type)
  
  # 检查是否有足够的细胞进行分析
  if (ncol(cell_type_cells) < 10) {
    cat(sprintf("警告: %s的细胞数量不足，跳过差异表达分析\n", cell_type))
    next
  }
  
  # 检查是否有多个分组
  if (length(unique(cell_type_cells$group)) < 2) {
    cat(sprintf("警告: %s只有一个分组，跳过差异表达分析\n", cell_type))
    next
  }
  
  # 进行差异表达分析
  de_genes <- FindMarkers(
    cell_type_cells,
    ident.1 = "treatment",
    ident.2 = "control",
    min.pct = 0.25,
    logfc.threshold = 0.5,
    verbose = FALSE
  ) %>%
    rownames_to_column("gene") %>%
    mutate(
      cell_type = cell_type,
      padj = p.adjust(p_val, method = "fdr")
    ) %>%
    filter(padj < 0.05) %>%
    arrange(padj)
  
  if (nrow(de_genes) > 0) {
    de_results[[cell_type]] <- de_genes
    cat(sprintf("在%s中鉴定到%d个差异表达基因\n", cell_type, nrow(de_genes)))
  } else {
    cat(sprintf("在%s中未鉴定到显著差异表达基因\n", cell_type))
  }
}

# 3.2 整合差异表达结果
if (length(de_results) > 0) {
  all_de_genes <- bind_rows(de_results)
  
  # 保存差异表达结果
  write.csv(all_de_genes, here(output_dir, "differential_expression_results.csv"), row.names = FALSE)
  
  # 差异表达基因统计
  de_summary <- all_de_genes %>%
    group_by(cell_type) %>%
    summarize(
      n_genes = n(),
      n_upregulated = sum(avg_log2FC > 0),
      n_downregulated = sum(avg_log2FC < 0)
    )
  
  write.csv(de_summary, here(output_dir, "differential_expression_summary.csv"), row.names = FALSE)
  
  cat("\n差异表达分析结果汇总:\n")
  print(de_summary)
} else {
  cat("未鉴定到任何差异表达基因\n")
}

# 3.3 差异表达基因可视化
if (length(de_results) > 0) {
  # 选择每个细胞类型的前10个差异表达基因
  top_de_genes <- all_de_genes %>%
    group_by(cell_type) %>%
    slice_head(n = 10) %>%
    pull(gene) %>%
    unique()
  
  # 热图
  if (length(top_de_genes) > 0) {
    p5 <- DoHeatmap(
      seurat_obj,
      features = top_de_genes,
      group.by = c("final_cell_type", "group"),
      size = 3
    ) +
      NoLegend() +
      ggtitle("差异表达基因热图")
    
    ggsave(
      filename = here(output_dir, "differential_expression_heatmap.png"),
      plot = p5,
      width = 15,
      height = 10,
      dpi = 300
    )
  }
  
  # 火山图
  # 为每个细胞类型创建火山图
  for (cell_type in names(de_results)) {
    de_data <- de_results[[cell_type]]
    
    if (nrow(de_data) > 0) {
      p6 <- de_data %>%
        ggplot(aes(x = avg_log2FC, y = -log10(padj))) +
        geom_point(aes(color = ifelse(padj < 0.05 & abs(avg_log2FC) > 0.5, "显著", "不显著")),
                  alpha = 0.7, size = 2) +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
        scale_color_manual(values = c("显著" = "red", "不显著" = "black")) +
        labs(
          title = sprintf("%s的差异表达基因火山图", cell_type),
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
        )
      
      ggsave(
        filename = here(output_dir, sprintf("%s_volcano_plot.png", gsub(" ", "_", cell_type))),
        plot = p6,
        width = 10,
        height = 8,
        dpi = 300
      )
    }
  }
}

# ----------------------
# 4. 功能富集分析
# ----------------------

cat("\n开始功能富集分析...\n")

# 4.1 进行GO和KEGG富集分析
if (length(de_results) > 0) {
  enrichment_results <- list()
  
  for (cell_type in names(de_results)) {
    de_data <- de_results[[cell_type]]
    
    if (nrow(de_data) > 0) {
      cat(sprintf("分析%s的功能富集...\n", cell_type))
      
      # 获取差异表达基因
      de_genes <- de_data$gene
      
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
        
        go_mf <- enrichGO(
          gene = gene_df$ENTREZID,
          OrgDb = org.Hs.eg.db,
          keyType = "ENTREZID",
          ont = "MF",
          pAdjustMethod = "fdr",
          qvalueCutoff = 0.05
        )
        
        go_cc <- enrichGO(
          gene = gene_df$ENTREZID,
          OrgDb = org.Hs.eg.db,
          keyType = "ENTREZID",
          ont = "CC",
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
        enrichment_results[[cell_type]] <- list(
          go_bp = go_bp,
          go_mf = go_mf,
          go_cc = go_cc,
          kegg = kegg
        )
        
        # 打印富集结果统计
        cat(sprintf("GO生物过程: %d个显著条目\n", ifelse(!is.null(go_bp), nrow(go_bp@result), 0)))
        cat(sprintf("GO分子功能: %d个显著条目\n", ifelse(!is.null(go_mf), nrow(go_mf@result), 0)))
        cat(sprintf("GO细胞组分: %d个显著条目\n", ifelse(!is.null(go_cc), nrow(go_cc@result), 0)))
        cat(sprintf("KEGG通路: %d个显著条目\n", ifelse(!is.null(kegg), nrow(kegg@result), 0)))
      } else {
        cat(sprintf("警告: %s的基因转换失败，跳过功能富集分析\n", cell_type))
      }
    }
  }
  
  # 4.2 保存富集分析结果
  if (length(enrichment_results) > 0) {
    # 创建富集结果目录
    enrichment_dir <- here(output_dir, "enrichment_results")
    dir.create(enrichment_dir, recursive = TRUE, showWarnings = FALSE)
    
    # 保存每个细胞类型的富集结果
    for (cell_type in names(enrichment_results)) {
      results <- enrichment_results[[cell_type]]
      cell_type_dir <- here(enrichment_dir, gsub(" ", "_", cell_type))
      dir.create(cell_type_dir, recursive = TRUE, showWarnings = FALSE)
      
      # 保存GO生物过程结果
      if (!is.null(results$go_bp) && nrow(results$go_bp@result) > 0) {
        write.csv(results$go_bp@result, here(cell_type_dir, "go_biological_process.csv"), row.names = FALSE)
      }
      
      # 保存GO分子功能结果
      if (!is.null(results$go_mf) && nrow(results$go_mf@result) > 0) {
        write.csv(results$go_mf@result, here(cell_type_dir, "go_molecular_function.csv"), row.names = FALSE)
      }
      
      # 保存GO细胞组分结果
      if (!is.null(results$go_cc) && nrow(results$go_cc@result) > 0) {
        write.csv(results$go_cc@result, here(cell_type_dir, "go_cell_component.csv"), row.names = FALSE)
      }
      
      # 保存KEGG通路结果
      if (!is.null(results$kegg) && nrow(results$kegg@result) > 0) {
        write.csv(results$kegg@result, here(cell_type_dir, "kegg_pathway.csv"), row.names = FALSE)
      }
    }
    
    # 4.3 富集结果可视化
    for (cell_type in names(enrichment_results)) {
      results <- enrichment_results[[cell_type]]
      cell_type_dir <- here(enrichment_dir, gsub(" ", "_", cell_type))
      
      # GO生物过程气泡图
      if (!is.null(results$go_bp) && nrow(results$go_bp@result) > 0) {
        p7 <- dotplot(results$go_bp, showCategory = 10, title = sprintf("%s的GO生物过程富集", cell_type))
        ggsave(
          filename = here(cell_type_dir, "go_biological_process_dotplot.png"),
          plot = p7,
          width = 12,
          height = 8,
          dpi = 300
        )
      }
      
      # KEGG通路气泡图
      if (!is.null(results$kegg) && nrow(results$kegg@result) > 0) {
        p8 <- dotplot(results$kegg, showCategory = 10, title = sprintf("%s的KEGG通路富集", cell_type))
        ggsave(
          filename = here(cell_type_dir, "kegg_pathway_dotplot.png"),
          plot = p8,
          width = 12,
          height = 8,
          dpi = 300
        )
      }
    }
  } else {
    cat("未进行功能富集分析\n")
  }
} else {
  cat("未进行功能富集分析\n")
}

# ----------------------
# 5. Treg细胞分析
# ----------------------

cat("\n开始Treg细胞分析...\n")

# 5.1 提取Treg细胞
# 定义Treg细胞标记基因
treg_markers <- c("FOXP3", "IL2RA", "CTLA4", "CD28", "TNFRSF18", "IKZF2")
treg_markers_present <- intersect(treg_markers, rownames(seurat_obj))

if (length(treg_markers_present) > 0) {
  # 计算Treg细胞得分
  seurat_obj$treg_score <- colMeans(GetAssayData(seurat_obj, assay = "SCT", slot = "data")[treg_markers_present, ])
  
  # 识别Treg细胞
  # 方法1: 基于预定义的细胞类型注释
  treg_cells_annotation <- seurat_obj$final_cell_type %in% c("调节性T细胞", "Treg")
  
  # 方法2: 基于Treg标记基因的高表达
  treg_threshold <- quantile(seurat_obj$treg_score, 0.75)
  treg_cells_score <- seurat_obj$treg_score > treg_threshold
  
  # 综合两种方法
  seurat_obj$is_treg <- treg_cells_annotation | treg_cells_score
  
  # 统计Treg细胞数量
  treg_count <- sum(seurat_obj$is_treg)
  treg_percent <- treg_count / ncol(seurat_obj) * 100
  
  cat(sprintf("鉴定到%d个Treg细胞，占总细胞的%.1f%%\n", treg_count, treg_percent))
  
  # Treg细胞在不同分组中的分布
  treg_by_group <- table(seurat_obj$is_treg, seurat_obj$group) %>%
    as.data.frame.matrix() %>%
    rownames_to_column("is_treg") %>%
    mutate(is_treg = ifelse(is_treg == "TRUE", "Treg细胞", "非Treg细胞"))
  
  write.csv(treg_by_group, here(output_dir, "treg_cell_distribution.csv"), row.names = FALSE)
  
  cat("\nTreg细胞在不同分组中的分布:\n")
  print(treg_by_group)
  
  # 5.2 Treg细胞可视化
  # Umap图（按Treg细胞着色）
  p9 <- DimPlot(seurat_obj, reduction = "umap", group.by = "is_treg", label = TRUE) +
    ggtitle("Treg细胞分布 (UMAP)") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey"))
  
  # Treg细胞标记基因表达热图
  p10 <- FeaturePlot(seurat_obj, features = treg_markers_present, ncol = 3)
  
  # Treg细胞得分分布
  p11 <- VlnPlot(seurat_obj, features = "treg_score", group.by = "group") +
    ggtitle("Treg细胞得分分布") +
    geom_hline(yintercept = treg_threshold, linetype = "dashed", color = "red")
  
  # Treg细胞比例条形图
  treg_prop <- table(seurat_obj$group, seurat_obj$is_treg) %>%
    prop.table(margin = 1) %>%
    as.data.frame() %>%
    filter(Var2 == TRUE) %>%
    mutate(Var1 = factor(Var1))
  
  p12 <- ggplot(treg_prop, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", fill = "red", alpha = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", Freq * 100)), vjust = -0.3) +
    labs(
      title = "Treg细胞比例",
      x = "分组",
      y = "比例"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # 保存Treg细胞可视化结果
  ggsave(
    filename = here(output_dir, "treg_cell_analysis.png"),
    plot = cowplot::plot_grid(p9, p12, p11, ncol = 1),
    width = 12,
    height = 18,
    dpi = 300
  )
  
  ggsave(
    filename = here(output_dir, "treg_marker_expression.png"),
    plot = p10,
    width = 15,
    height = 10,
    dpi = 300
  )
  
  # 5.3 Treg细胞差异表达分析
  if (treg_count > 0) {
    # 提取Treg细胞
    treg_seurat <- subset(seurat_obj, subset = is_treg == TRUE)
    
    # 检查是否有足够的Treg细胞进行分析
    if (ncol(treg_seurat) >= 10) {
      # 检查是否有多个分组
      if (length(unique(treg_seurat$group)) >= 2) {
        cat("进行Treg细胞差异表达分析...\n")
        
        # Treg细胞差异表达分析
        treg_de <- FindMarkers(
          treg_seurat,
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
        
        if (nrow(treg_de) > 0) {
          cat(sprintf("在Treg细胞中鉴定到%d个差异表达基因\n", nrow(treg_de)))
          
          # 保存Treg细胞差异表达结果
          write.csv(treg_de, here(output_dir, "treg_cell_differential_expression.csv"), row.names = FALSE)
          
          # Treg细胞差异表达基因火山图
          p13 <- treg_de %>%
            ggplot(aes(x = avg_log2FC, y = -log10(padj))) +
            geom_point(aes(color = ifelse(padj < 0.05 & abs(avg_log2FC) > 0.5, "显著", "不显著")),
                      alpha = 0.7, size = 2) +
            geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray") +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
            scale_color_manual(values = c("显著" = "red", "不显著" = "black")) +
            labseller::labs(
              title = "Treg细胞差异表达基因火山图",
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
            filename = here(output_dir, "treg_cell_volcano_plot.png"),
            plot = p13,
            width = 10,
            height = 8,
            dpi = 300
          )
          
          # Treg细胞差异表达基因热图
          if (nrow(treg_de) > 0) {
            top_treg_de <- head(treg_de, 20)
            p14 <- DoHeatmap(
              treg_seurat,
              features = top_treg_de$gene,
              group.by = "group",
              size = 3
            ) +
              NoLegend() +
              ggtitle("Treg细胞差异表达基因热图")
            
            ggsave(
              filename = here(output_dir, "treg_cell_heatmap.png"),
              plot = p14,
              width = 12,
              height = 8,
              dpi = 300
            )
          }
        } else {
          cat("在Treg细胞中未鉴定到显著差异表达基因\n")
        }
      } else {
        cat("Treg细胞只有一个分组，跳过差异表达分析\n")
      }
    } else {
      cat("Treg细胞数量不足，跳过差异表达分析\n")
    }
  } else {
    cat("未鉴定到Treg细胞，跳过Treg细胞分析\n")
  }
} else {
  cat("未检测到Treg细胞标记基因，跳过Treg细胞分析\n")
}

# ----------------------
# 6. 保存分析结果
# ----------------------

cat("\n保存分析结果...\n")

# 6.1 保存最终的Seurat对象
saveRDS(seurat_obj, here(output_dir, "annotated_seurat_object.rds"))

# 6.2 创建分析报告
analysis_summary <- data.frame(
  metric = c(
    "总细胞数",
    "总基因数",
    "细胞类型数",
    "Treg细胞数",
    "Treg细胞比例(%)",
    "差异表达基因总数"
  ),
  value = c(
    ncol(seurat_obj),
    nrow(seurat_obj),
    length(unique(seurat_obj$final_cell_type)),
    ifelse("is_treg" %in% colnames(seurat_obj@meta.data), sum(seurat_obj$is_treg), NA),
    ifelse("is_treg" %in% colnames(seurat_obj@meta.data), sum(seurat_obj$is_treg)/ncol(seurat_obj)*100, NA),
    if (length(de_results) > 0) nrow(all_de_genes) else 0
  )
)

write.csv(analysis_summary, here(output_dir, "analysis_summary.csv"), row.names = FALSE)

# 打印分析报告
cat("\n分析报告:\n")
print(analysis_summary)

# ----------------------
# 7. 生成分析报告
# ----------------------

cat("\n生成分析报告...\n")

# 创建分析报告
report_content <- glue(
  "# CD28-Treg-IBD-Study: 单细胞RNA测序数据分析报告\n\n",
  
  "## 分析概述\n",
  "- **分析日期**: {Sys.Date()}\n",
  "- **数据类型**: 单细胞RNA测序数据\n",
  "- **分析步骤**: 细胞类型注释、差异表达分析、功能富集分析、Treg细胞分析\n\n",
  
  "## 数据概况\n",
  "| 指标 | 数值 |\n",
  "|------|------|\n",
  "{paste(sprintf('| %s | %s |', analysis_summary$metric, analysis_summary$value), collapse = '\\n')}\n\n",
  
  "## 细胞类型注释\n",
  "### 注释方法\n",
  "- **SingleR**: 使用HumanPrimaryCellAtlasData作为参考数据集\n",
  "- **标记基因**: 基于已知细胞类型标记基因的表达模式\n\n",
  
  "### 细胞类型分布\n",
  "{cell_type_table <- table(seurat_obj$final_cell_type) %>% as.data.frame()\n",
  "colnames(cell_type_table) <- c('细胞类型', '数量')\n",
  "cell_type_table$比例 <- cell_type_table$数量 / sum(cell_type_table$数量) * 100\n",
  "paste(sprintf('| %s | %d | %.1f%% |', cell_type_table$细胞类型, cell_type_table$数量, cell_type_table$比例), collapse = '\\n')}\n\n",
  
  "## 差异表达分析\n",
  "{if (length(de_results) > 0) {\n",
  "  paste(sprintf('- 共鉴定到%d个差异表达基因\\n- 涉及%d种细胞类型', nrow(all_de_genes), length(de_results)), collapse = '\\n')\n",
  "} else {\n",
  "  '未鉴定到任何差异表达基因'\n",
  "}}\n\n",
  
  "{if (length(de_results) > 0) {\n",
  "  paste(sprintf('| 细胞类型 | 差异表达基因数 | 上调基因数 | 下调基因数 |', \n",
  "                de_summary$cell_type, \n",
  "                de_summary$n_genes, \n",
  "                de_summary$n_upregulated, \n",
  "                de_summary$n_downregulated), collapse = '\\n')\n",
  "}}\n\n",
  
  "## Treg细胞分析\n",
  "{if ('is_treg' %in% colnames(seurat_obj@meta.data)) {\n",
  "  treg_count <- sum(seurat_obj$is_treg)\n",
  "  treg_percent <- treg_count / ncol(seurat_obj) * 100\n",
  "  \n",
  "  treg_by_group <- table(seurat_obj$group, seurat_obj$is_treg) %>% as.data.frame.matrix()\n",
  "  treg_by_group$total <- rowSums(treg_by_group)\n",
  "  treg_by_group$treg_percent <- treg_by_group$TRUE / treg_by_group$total * 100\n",
  "  \n",
  "  paste(c(\n",
  "    sprintf('- 共鉴定到%d个Treg细胞，占总细胞的%.1f%%', treg_count, treg_percent),\n",
  "    '',\n",
  "    'Treg细胞在不同分组中的分布:',\n",
  "    sprintf('- %s组: %d个Treg细胞 (%.1f%%)', rownames(treg_by_group)[1], treg_by_group$TRUE[1], treg_by_group$treg_percent[1]),\n",
  "    sprintf('- %s组: %d个Treg细胞 (%.1f%%)', rownames(treg_by_group)[2], treg_by_group$TRUE[2], treg_by_group$treg_percent[2])\n",
  "  ), collapse = '\\n')\n",
  "} else {\n",
  "  '未进行Treg细胞分析'\n",
  "}}\n\n",
  
  "## 功能富集分析\n",
  "{if (length(enrichment_results) > 0) {\n",
  "  paste(sprintf('对%d种细胞类型的差异表达基因进行了功能富集分析', length(enrichment_results)), collapse = '\\n')\n",
  "} else {\n",
  "  '未进行功能富集分析'\n",
  "}}\n\n",
  
  "## 输出文件\n",
  "- `annotated_seurat_object.rds`: 注释后的Seurat对象\n",
  "- `cell_type_counts.csv`: 细胞类型统计\n",
  "- `analysis_summary.csv`: 分析摘要\n",
  "- `cell_type_annotation_umap.png`: 细胞类型注释UMAP图\n",
  "- `marker_genes_heatmap.png`: 标记基因热图\n",
  "\n",
  "{if (length(de_results) > 0) {\n",
  "  paste(c(\n",
  "    '- `differential_expression_results.csv`: 差异表达基因结果',\n",
  "    '- `differential_expression_summary.csv`: 差异表达基因摘要',\n",
  "    '- `differential_expression_heatmap.png`: 差异表达基因热图',\n",
  "    '- 各细胞类型的差异表达基因火山图'\n",
  "  ), collapse = '\\n')\n",
  "} else {\n",
  "  ''\n",
  "}}\n",
  "\n",
  "{if (length(enrichment_results) > 0) {\n",
  "  paste(c(\n",
  "    '- `enrichment_results/`: 功能富集分析结果目录',\n",
  "    '- 各细胞类型的GO和KEGG富集分析结果',\n",
  "    '- 各细胞类型的GO和KEGG富集分析可视化'\n",
  "  ), collapse = '\\n')\n",
  "} else {\n",
  "  ''\n",
  "}}\n",
  "\n",
  "{if ('is_treg' %in% colnames(seurat_obj@meta.data)) {\n",
  "  paste(c(\n",
  "    '- `treg_cell_distribution.csv`: Treg细胞分布',\n",
  "    '- `treg_cell_analysis.png`: Treg细胞分析图',\n",
  "    '- `treg_marker_expression.png`: Treg标记基因表达图'\n",
  "  ), collapse = '\\n')\n",
  "} else {\n",
  "  ''\n",
  "}}\n\n",
  
  "## 下一步分析建议\n",
  "- 细胞轨迹分析，探究细胞分化路径\n",
  "- 基因调控网络分析，识别关键调控因子\n",
  "- Treg细胞亚型分析，探究功能异质性\n",
  "- 多组学整合分析，结合基因组、表观基因组数据\n",
  "- 实验验证关键发现\n\n",
  
  "---\n",
  "*报告生成时间: {Sys.time()}*"
)

# 保存报告
writeLines(
  report_content,
  here(output_dir, "single_cell_analysis_report.md")
)

# ----------------------
# 8. 完成分析
# ----------------------

cat("单细胞RNA测序数据分析完成！\n")
cat("分析结果已保存到:", output_dir, "\n")
cat("主要输出文件包括:\n")
cat("- annotated_seurat_object.rds\n")
cat("- analysis_summary.csv\n")
cat("- cell_type_annotation_umap.png\n")
cat("- single_cell_analysis_report.md\n")

#' CD28-Treg-IBD-Study: LD评分回归分析
#' 
#' 本脚本用于通过LD评分回归(LDSC)分析CD28+Treg细胞特征和IBD之间的遗传相关性，
#' 并评估遗传力和遗传相关性。

# 加载必要的包
library(ldsc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(here)
library(glue)
library(data.table)

# 设置输出目录
output_dir <- here("results", "LDSC_Analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------
# 1. 数据准备
# ----------------------

cat("开始数据准备...\n")

# 注意：在实际分析中，应从data目录加载真实数据
# 这里使用示例数据进行演示

# 设置LD评分文件路径（实际分析中应使用真实的LD评分文件）
ldsc_dir <- here("data", "external", "ldsc")
dir.create(ldsc_dir, recursive = TRUE, showWarnings = FALSE)

# 模拟创建LD评分文件（仅用于演示，实际分析中应使用真实数据）
# 在实际分析中，应从LDSC官网下载预计算的LD评分文件
simulate_ld_scores <- function(output_dir) {
  # 创建模拟的LD评分文件
  n_snps <- 10000
  n_bins <- 100
  
  for (chr in 1:22) {
    # 创建SNP信息
    snp_info <- data.frame(
      SNP = paste0("rs", sample(100000:999999, n_snps)),
      CHR = chr,
      BP = sample(1:100000000, n_snps),
      CM = sample(1:1000, n_snps),
      MAF = runif(n_snps, 0.01, 0.5),
      ALT_FREQ = runif(n_snps, 0.01, 0.5),
      REF = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
      ALT = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE)
    )
    
    # 确保REF和ALT不同
    snp_info$ALT[snp_info$REF == snp_info$ALT] <- sample(c("A", "T", "C", "G"), sum(snp_info$REF == snp_info$ALT), replace = TRUE)
    
    # 创建LD评分
    ld_scores <- matrix(runif(n_snps * n_bins, 0, 2), nrow = n_snps, ncol = n_bins)
    
    # 保存SNP信息
    write.table(
      snp_info,
      file = here(output_dir, paste0("chr", chr, ".snp")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    # 保存LD评分
    write.table(
      ld_scores,
      file = here(output_dir, paste0("chr", chr, ".l2.ldscore.gz")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }
  
  cat("模拟LD评分文件创建完成\n")
}

# 模拟创建LD评分文件
simulate_ld_scores(ldsc_dir)

# 创建示例GWAS汇总统计数据
create_gwas_summary <- function(pheno_name, n_snps = 10000) {
  set.seed(1234)
  
  gwas_data <- data.frame(
    SNP = paste0("rs", sample(100000:999999, n_snps)),
    CHR = sample(1:22, n_snps, replace = TRUE),
    BP = sample(1:100000000, n_snps),
    A1 = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
    A2 = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
    N = sample(5000:50000, n_snps, replace = TRUE),
    Z = rnorm(n_snps, mean = 0, sd = 1),
    P = pnorm(-abs(rnorm(n_snps, mean = 0, sd = 1))) * 2
  )
  
  # 确保A1和A2不同
  gwas_data$A2[gwas_data$A1 == gwas_data$A2] <- sample(c("A", "T", "C", "G"), sum(gwas_data$A1 == gwas_data$A2), replace = TRUE)
  
  # 计算BETA和SE
  gwas_data$BETA <- gwas_data$Z / sqrt(gwas_data$N)
  gwas_data$SE <- 1 / sqrt(gwas_data$N)
  
  # 保存GWAS汇总统计数据
  output_file <- here(output_dir, paste0(pheno_name, ".sumstats.gz"))
  write.table(
    gwas_data,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  return(list(
    data = gwas_data,
    file = output_file
  ))
}

# 创建CD28+Treg细胞特征的GWAS汇总统计数据
cd28_treg_gwas <- create_gwas_summary("cd28_treg_ratio")

# 创建IBD的GWAS汇总统计数据
ibd_gwas <- create_gwas_summary("ibd")

# 创建其他免疫相关性状的GWAS汇总统计数据（用于遗传相关性分析）
other_traits <- list(
  "crohn_disease" = "克罗恩病",
  "ulcerative_colitis" = "溃疡性结肠炎",
  "inflammatory_bowel_disease" = "炎症性肠病",
  "rheumatoid_arthritis" = "类风湿关节炎",
  "type_1_diabetes" = "1型糖尿病",
  "asthma" = "哮喘"
)

other_gwas <- lapply(names(other_traits), function(trait) {
  create_gwas_summary(trait)
})

names(other_gwas) <- names(other_traits)

# ----------------------
# 2. 遗传力分析
# ----------------------

cat("\n开始遗传力分析...\n")

# 2.1 函数：运行LDSC遗传力分析
run_ldsc_heritability <- function(sumstats_file, ldsc_dir, output_prefix) {
  # 设置输出目录
  res_dir <- file.path(output_dir, output_prefix)
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 运行LDSC遗传力分析
  h2_results <- ldsc::ldsc(
    h2 = sumstats_file,
    ref_ld_chr = file.path(ldsc_dir, "chr@.l2.ldscore.gz"),
    w_ld_chr = file.path(ldsc_dir, "chr@.l2.ldscore.gz"),
    out = file.path(res_dir, output_prefix),
    overwrite = TRUE
  )
  
  return(list(
    results = h2_results,
    output_dir = res_dir
  ))
}

# 2.2 分析CD28+Treg细胞特征的遗传力
cd28_treg_h2 <- run_ldsc_heritability(
  sumstats_file = cd28_treg_gwas$file,
  ldsc_dir = ldsc_dir,
  output_prefix = "cd28_treg_ratio_h2"
)

# 2.3 分析IBD的遗传力
ibd_h2 <- run_ldsc_heritability(
  sumstats_file = ibd_gwas$file,
  ldsc_dir = ldsc_dir,
  output_prefix = "ibd_h2"
)

# 2.4 分析其他性状的遗传力
other_h2 <- lapply(names(other_gwas), function(trait) {
  run_ldsc_heritability(
    sumstats_file = other_gwas[[trait]]$file,
    ldsc_dir = ldsc_dir,
    output_prefix = paste0(trait, "_h2")
  )
})

names(other_h2) <- names(other_gwas)

# 2.5 汇总遗传力结果
h2_summary <- data.frame(
  trait = c("CD28+Treg细胞比例", "IBD", unlist(other_traits)),
  h2 = c(
    cd28_treg_h2$results$h2,
    ibd_h2$results$h2,
    sapply(other_h2, function(x) x$results$h2)
  ),
  h2_se = c(
    cd28_treg_h2$results$h2_se,
    ibd_h2$results$h2_se,
    sapply(other_h2, function(x) x$results$h2_se)
  ),
  h2_z = c(
    cd28_treg_h2$results$h2 / cd28_treg_h2$results$h2_se,
    ibd_h2$results$h2 / ibd_h2$results$h2_se,
    sapply(other_h2, function(x) x$results$h2 / x$results$h2_se)
  ),
  h2_p = c(
    2 * pnorm(-abs(cd28_treg_h2$results$h2 / cd28_treg_h2$results$h2_se)),
    2 * pnorm(-abs(ibd_h2$results$h2 / ibd_h2$results$h2_se)),
    sapply(other_h2, function(x) 2 * pnorm(-abs(x$results$h2 / x$results$h2_se)))
  )
) %>%
  mutate(
    h2_lower = h2 - 1.96 * h2_se,
    h2_upper = h2 + 1.96 * h2_se,
    significance = ifelse(h2_p < 0.05, "显著", "不显著")
  )

# 保存遗传力结果
write.csv(
  h2_summary,
  here(output_dir, "heritability_summary.csv"),
  row.names = FALSE
)

# 打印遗传力结果
cat("遗传力分析结果:\n")
print(h2_summary %>% select(trait, h2, h2_se, h2_p, significance))

# ----------------------
# 3. 遗传相关性分析
# ----------------------

cat("\n开始遗传相关性分析...\n")

# 3.1 函数：运行LDSC遗传相关性分析
run_ldsc_genetic_correlation <- function(sumstats1, sumstats2, ldsc_dir, output_prefix) {
  # 设置输出目录
  res_dir <- file.path(output_dir, output_prefix)
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 运行LDSC遗传相关性分析
  rg_results <- ldsc::ldsc(
    rg = paste(sumstats1, sumstats2, sep = ","),
    ref_ld_chr = file.path(ldsc_dir, "chr@.l2.ldscore.gz"),
    w_ld_chr = file.path(ldsc_dir, "chr@.l2.ldscore.gz"),
    out = file.path(res_dir, output_prefix),
    overwrite = TRUE
  )
  
  return(list(
    results = rg_results,
    output_dir = res_dir
  ))
}

# 3.2 分析CD28+Treg细胞特征与IBD的遗传相关性
cd28_treg_ibd_rg <- run_ldsc_genetic_correlation(
  sumstats1 = cd28_treg_gwas$file,
  sumstats2 = ibd_gwas$file,
  ldsc_dir = ldsc_dir,
  output_prefix = "cd28_treg_ibd_rg"
)

# 3.3 分析CD28+Treg细胞特征与其他免疫相关性状的遗传相关性
cd28_treg_other_rg <- lapply(names(other_gwas), function(trait) {
  run_ldsc_genetic_correlation(
    sumstats1 = cd28_treg_gwas$file,
    sumstats2 = other_gwas[[trait]]$file,
    ldsc_dir = ldsc_dir,
    output_prefix = paste0("cd28_treg_", trait, "_rg")
  )
})

names(cd28_treg_other_rg) <- names(other_gwas)

# 3.4 汇总遗传相关性结果
rg_summary <- data.frame(
  trait1 = rep("CD28+Treg细胞比例", length(other_gwas) + 1),
  trait2 = c("IBD", unlist(other_traits)),
  rg = c(
    cd28_treg_ibd_rg$results$rg,
    sapply(cd28_treg_other_rg, function(x) x$results$rg)
  ),
  rg_se = c(
    cd28_treg_ibd_rg$results$rg_se,
    sapply(cd28_treg_other_rg, function(x) x$results$rg_se)
  ),
  rg_z = c(
    cd28_treg_ibd_rg$results$rg / cd28_treg_ibd_rg$results$rg_se,
    sapply(cd28_treg_other_rg, function(x) x$results$rg / x$results$rg_se)
  ),
  rg_p = c(
    2 * pnorm(-abs(cd28_treg_ibd_rg$results$rg / cd28_treg_ibd_rg$results$rg_se)),
    sapply(cd28_treg_other_rg, function(x) 2 * pnorm(-abs(x$results$rg / x$results$rg_se)))
  )
) %>%
  mutate(
    rg_lower = rg - 1.96 * rg_se,
    rg_upper = rg + 1.96 * rg_se,
    significance = case_when(
      rg_p < 0.001 ~ "***",
      rg_p < 0.01 ~ "**",
      rg_p < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# 保存遗传相关性结果
write.csv(
  rg_summary,
  here(output_dir, "genetic_correlation_summary.csv"),
  row.names = FALSE
)

# 打印遗传相关性结果
cat("\n遗传相关性分析结果:\n")
print(rg_summary %>% select(trait1, trait2, rg, rg_se, rg_p, significance))

# ----------------------
# 4. 结果可视化
# ----------------------

cat("\n开始结果可视化...\n")

# 4.1 遗传力森林图
h2_forest <- h2_summary %>%
  ggplot(aes(x = h2, y = trait)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_point(size = 3, aes(color = significance)) +
  geom_errorbarh(aes(xmin = h2_lower, xmax = h2_upper, color = significance), height = 0.2) +
  geom_text(aes(label = sprintf("%.3f (%.3f, %.3f)", h2, h2_lower, h2_upper)), 
            hjust = -0.1, size = 3.5) +
  scale_color_manual(values = c("显著" = "#e74c3c", "不显著" = "#3498db")) +
  labs(
    title = "各性状的遗传力估计",
    x = "遗传力 (h² ± 95% CI)",
    y = "性状",
    color = "显著性"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  expand_limits(x = c(min(h2_summary$h2_lower) * 1.2, max(h2_summary$h2_upper) * 1.2))

# 保存遗传力森林图
ggsave(
  filename = here(output_dir, "heritability_forest_plot.png"),
  plot = h2_forest,
  width = 12,
  height = 8,
  dpi = 300
)

# 4.2 遗传相关性森林图
rg_forest <- rg_summary %>%
  ggplot(aes(x = rg, y = trait2)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_point(size = 3, aes(color = significance)) +
  geom_errorbarh(aes(xmin = rg_lower, xmax = rg_upper, color = significance), height = 0.2) +
  geom_text(aes(label = sprintf("%.3f (%.3f, %.3f) %s", rg, rg_lower, rg_upper, significance)), 
            hjust = -0.1, size = 3.5) +
  scale_color_manual(values = c("***" = "#c0392b", "**" = "#e74c3c", "*" = "#e67e22", "ns" = "#3498db")) +
  labs(
    title = "CD28+Treg细胞比例与其他性状的遗传相关性",
    x = "遗传相关性 (rg ± 95% CI)",
    y = "性状",
    color = "显著性"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  expand_limitseller::xlim(c(min(rg_summary$rg_lower) * 1.2, max(rg_summary$rg_upper) * 1.2))

# 保存遗传相关性森林图
ggsave(
  filename = here(output_dir, "genetic_correlation_forest_plot.png"),
  plot = rg_forest,
  width = 12,
  height = 8,
  dpi = 300
)

# 4.3 遗传相关性热图
# 创建遗传相关性矩阵
traits <- c("CD28+Treg细胞比例", "IBD", unlist(other_traits))
n_traits <- length(traits)

# 初始化遗传相关性矩阵
rg_matrix <- matrix(NA, nrow = n_traits, ncol = n_traits, dimnames = list(traits, traitseller::s))
p_matrix <- matrix(NA, nrow = n_traits, ncol = n_traits, dimnames = list(traits, traiteller::s))

# 填充对角线（遗传力）
for (i in 1:n_traits) {
  rg_matrix[i, i] <- h2_summary$h2[i]
  p_matrix[i, i] <- h2_summary$h2_p[i]
}

# 填充CD28+Treg细胞比例与其他性状的遗传相关性
for (i in 2:n_traits) {
  rg_matrix[1, i] <- rg_summary$rg[i-1]
  rg_matrix[i, 1] <- rg_summary$rg[i-1]
  p_matrix[1, i] <- rg_summary$rg_p[i-1]
  p_matrix[i, 1] <- rg_summary$rg_p[i-1]
}

# 为其他性状对填充模拟的遗传相关性
for (i in 2:n_traits) {
  for (j in i:n_traits) {
    if (i == j) next
    
    # 模拟遗传相关性
    rg <- runif(1, -0.5, 0.8)
    rg_se <- runif(1, 0.05, 0.3)
    rg_z <- rg / rg_se
    rg_p <- 2 * pnorm(-abs(rg_z))
    
    rg_matrix[i, j] <- rg
    rg_matrix[j, i] <- rg
    p_matrix[i, j] <- rg_p
    p_matrix[j, i] <- rg_p
  }
}

# 创建热图数据
heatmap_data <- reshape2::melt(rg_matrix) %>%
  rename(trait1 = Var1, trait2 = Var2, rg = value) %>%
  left_join(
    reshape2::melt(p_matrix) %>% rename(trait1 = Var1, trait2 = Var2, p = value),
    by = c("trait1", "trait2")
  ) %>%
  mutate(
    significance = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    label = ifelse(trait1 == trait2, 
                   sprintf("h²=%.2f", rg), 
                   sprintf("%.2f%s", rg, significance))
  )

# 创建热图
rg_heatmap <- ggplot(heatmap_data, aes(x = trait2, y = trait1, fill = rg)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3.5) +
  scale_fill_gradient2(
    low = "#2980b9",
    mid = "white",
    high = "#e74c3c",
    midpoint = 0,
    limits = c(-1, 1),
    name = "遗传相关性"
  ) +
  labs(
    title = "免疫相关性状的遗传相关性矩阵",
    x = "性状",
    y = "性状"
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

# 保存热图
ggsave(
  filename = here(output_dir, "genetic_correlation_heatmap.png"),
  plot = rg_heatmap,
  width = 14,
  height = 12,
  dpi = 300
)

# ----------------------
# 5. 分区遗传力分析
# ----------------------

cat("\n开始分区遗传力分析...\n")

# 5.1 函数：运行分区遗传力分析
run_partitioned_heritability <- function(sumstats_file, ldsc_dir, annot_dir, output_prefix) {
  # 设置输出目录
  res_dir <- file.path(output_dir, output_prefix)
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 运行分区遗传力分析
  partitioned_results <- ldsc::ldsc(
    h2 = sumstats_file,
    ref_ld_chr = file.path(ldsc_dir, "chr@.l2.ldscore.gz"),
    w_ld_chr = file.path(ldsc_dir, "chr@.l2.ldscore.gz"),
    annot = file.path(annot_dir, "chr@.annot.gz"),
    print_coefficients = TRUE,
    out = file.path(res_dir, output_prefix),
    overwrite = TRUE
  )
  
  return(list(
    results = partitioned_results,
    output_dir = res_dir
  ))
}

# 5.2 模拟功能注释文件（仅用于演示）
simulate_annotations <- function(output_dir) {
  # 创建模拟的功能注释文件
  n_snps <- 10000
  annotations <- c("coding", "promoter", "enhancer", "ctcf_binding", "tf_binding", "repressed")
  
  for (chr in 1:22) {
    # 创建注释矩阵（每行是一个SNP，每列是一个注释）
    annot_matrix <- matrix(
      sample(c(0, 1), n_snps * length(annotations), replace = TRUE, prob = c(0.9, 0.1)),
      nrow = n_snps,
      ncol = length(annotations),
      dimnames = list(NULL, annotations)
    )
    
    # 保存注释文件
    write.table(
      annot_matrix,
      file = here(output_dir, paste0("chr", chr, ".annot.gz")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  cat("模拟功能注释文件创建完成\n")
}

# 创建功能注释目录
annot_dir <- here(ldsc_dir, "annotations")
dir.create(annot_dir, recursive = TRUE, showWarnings = FALSE)

# 模拟功能注释文件
simulate_annotations(annot_dir)

# 5.3 分析CD28+Treg细胞特征的分区遗传力
cd28_treg_partitioned <- run_partitioned_heritability(
  sumstats_file = cd28_treg_gwas$file,
  ldsc_dir = ldsc_dir,
  annot_dir = annot_dir,
  output_prefix = "cd28_treg_partitioned_h2"
)

# 5.4 分析IBD的分区遗传力
ibd_partitioned <- run_partitioned_heritability(
  sumstats_file = ibd_gwas$file,
  ldsc_dir = ldsc_dir,
  annot_dir = annot_dir,
  output_prefix = "ibd_partitioned_h2"
)

# 5.5 汇总分区遗传力结果
partitioned_summary <- data.frame(
  annotation = colnames(cd28_treg_partitioned$results$coef),
  cd28_treg_coef = cd28_treg_partitioned$results$coef,
  cd28_treg_se = cd28_treg_partitioned$results$coef_se,
  cd28_treg_p = cd28_treg_partitioned$results$p,
  ibd_coef = ibd_partitioned$results$coef,
  ibd_se = ibd_partitioned$results$coef_se,
  ibd_p = ibd_partitioned$results$p
) %>%
  mutate(
    cd28_treg_significance = case_when(
      cd28_treg_p < 0.001 ~ "***",
      cd28_treg_p < 0.01 ~ "**",
      cd28_treg_p < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    ibd_significance = case_when(
      ibd_p < 0.001 ~ "***",
      ibd_p < 0.01 ~ "**",
      ibd_p < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# 保存分区遗传力结果
write.csv(
  partitioned_summary,
  here(output_dir, "partitioned_heritability_summary.csv"),
  row.names = FALSE
)

# 打印分区遗传力结果
cat("\n分区遗传力分析结果:\n")
print(partitioned_summary %>% select(annotation, cd28_treg_coef, cd28_treg_p, cd28_treg_significance, ibd_coef, ibd_p, ibd_significance))

# 5.6 分区遗传力可视化
partitioned_plot <- partitioned_summary %>%
  pivot_longer(
    cols = c(cd28_treg_coef, ibd_coef),
    names_to = "trait",
    values_to = "coef"
  ) %>%
  pivot_longer(
    cols = c(cd28_treg_se, ibd_se),
    names_to = "trait_se",
    values_to = "se"
  ) %>%
  filter(substr(trait, 1, nchar(trait)-5) == substr(trait_se, 1, nchar(trait_se)-3)) %>%
  pivot_longer(
    cols = c(cd28_treg_significance, ibd_significance),
    names_to = "trait_sig",
    values_to = "significance"
  ) %>%
  filter(substr(trait, 1, nchar(trait)-5) == substr(trait_sig, 1, nchar(trait_sig)-13)) %>%
  mutate(
    trait = ifelse(trait == "cd28_treg_coef", "CD28+Treg细胞比例", "IBD"),
    ci_lower = coef - 1.96 * se,
    ci_upper = coef + 1.96 * se
  ) %>%
  ggplot(aes(x = coef, y = annotation, color = trait)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
  geom_text(aes(label = significance), hjust = -0.1, size = 4) +
  scale_color_manual(values = c("#3498db", "#e74c3c")) +
  labseller::labs(
    title = "功能注释区域的遗传力富集",
    x = "遗传力富集系数",
    y = "功能注释",
    color = "性状"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

# 保存分区遗传力图
ggsave(
  filename = here(output_dir, "partitioned_heritability_plot.png"),
  plot = partitioned_plot,
  width = 12,
  height = 8,
  dpi = 300
)

# ----------------------
# 6. 生成分析报告
# ----------------------

cat("\n生成分析报告...\n")

# 创建分析报告
report_content <- glue(
  "# CD28-Treg-IBD-Study: LD评分回归分析报告\n\n",
  
  "## 分析概述\n",
  "- **分析日期**: {Sys.Date()}\n",
  "- **主要性状**: CD28+Treg细胞比例\n",
  "- **比较性状**: IBD及其他免疫相关疾病\n",
  "- **分析方法**: LD评分回归 (LDSC)\n\n",
  
  "## 数据概况\n",
  "- **CD28+Treg细胞比例GWAS**: {nrow(cd28_treg_gwas$data)}个SNP\n",
  "- **IBD GWAS**: {nrow(ibd_gwas$data)}个SNP\n",
  "- **其他免疫相关性状**: {length(other_gwas)}个\n",
  "- **LD评分文件**: 模拟数据，覆盖22条染色体\n\n",
  
  "## 遗传力分析结果\n",
  "### CD28+Treg细胞比例的遗传力\n",
  "- **遗传力 (h²)**: {sprintf('%.3f', cd28_treg_h2$results$h2)} ± {sprintf('%.3f', cd28_treg_h2$results$h2_se)}\n",
  "- **Z值**: {sprintf('%.3f', cd28_treg_h2$results$h2 / cd28_treg_h2$results$h2_se)}\n",
  "- **P值**: {sprintf('%.3e', 2 * pnorm(-abs(cd28_treg_h2$results$h2 / cd28_treg_h2$results$h2_se)))} {ifelse(2 * pnorm(-abs(cd28_treg_h2$results$h2 / cd28_treg_h2$results$h2_se)) < 0.05, '**显著**', '**不显著**')}\n\n",
  
  "### IBD的遗传力\n",
  "- **遗传力 (h²)**: {sprintf('%.3f', ibd_h2$results$h2)} ± {sprintf('%.3f', ibd_h2$results$h2_se)}\n",
  "- **Z值**: {sprintf('%.3f', ibd_h2$results$h2 / ibd_h2$results$h2_se)}\n",
  "- **P值**: {sprintf('%.3e', 2 * pnorm(-abs(ibd_h2$results$h2 / ibd_h2$results$h2_se)))} {ifelse(2 * pnorm(-abs(ibd_h2$results$h2 / ibd_h2$results$h2_se)) < 0.05, '**显著**', '**不显著**')}\n\n",
  
  "## 遗传相关性分析结果\n",
  "### CD28+Treg细胞比例与IBD的遗传相关性\n",
  "- **遗传相关性 (rg)**: {sprintf('%.3f', cd28_treg_ibd_rg$results$rg)} ± {sprintf('%.3f', cd28_treg_ibd_rg$results$rg_se)}\n",
  "- **Z值**: {sprintf('%.3f', cd28_treg_ibd_rg$results$rg / cd28_treg_ibd_rg$results$rg_se)}\n",
  "- **P值**: {sprintf('%.3e', 2 * pnorm(-abs(cd28_treg_ibd_rg$results$rg / cd28_treg_ibd_rg$results$rg_se)))} {ifelse(2 * pnorm(-abs(cd28_treg_ibd_rg$results$rg / cd28_treg_ibd_rg$results$rg_se)) < 0.05, '**显著**', '**不显著**')}\n\n",
  
  "### CD28+treg细胞比例与其他免疫相关性状的遗传相关性\n",
  "{paste(lapply(names(other_gwas), function(trait) {\n",
  "  rg <- cd28_treg_other_rg[[trait]]$results$rg\n",
  "  rg_se <- cd28_treg_other_rg[[trait]]$results$rg_se\n",
  "  rg_p <- 2 * pnorm(-abs(rg / rg_se))\n",
  "  sig <- ifelse(rg_p < 0.05, '**显著**', '**不显著**')\n",
  "  sprintf('- **%s**: %.3f ± %.3f (p = %.3e) %s', other_traits[[trait]], rg, rg_se, rg_p, sig)\n",
  "}), collapse = '\\n')}\n\n",
  
  "## 分区遗传力分析结果\n",
  "### CD28+Treg细胞比例的功能富集\n",
  "{significant_annotations_cd28 <- partitioned_summary %>% filter(cd28_treg_p < 0.05) %>% pull(annotation)\n",
  "if (length(significant_annotations_cd28) > 0) {\n",
  "  paste(sprintf('- **%s区域**存在显著遗传力富集', significant_annotations_cd28), collapse = '\\n')\n",
  "} else {\n",
  "  '未发现显著的遗传力富集区域'\n",
  "}}\n\n",
  
  "### IBD的功能富集\n",
  "{significant_annotations_ibd <- partitioned_summary %>% filter(ibd_p < 0.05) %>% pull(annotation)\n",
  "if (length(significant_annotations_ibd) > 0) {\n",
  "  paste(sprintf('- **%s区域**存在显著遗传力富集', significant_annotations_ibd), collapse = '\\n')\n",
  "} else {\n",
  "  '未发现显著的遗传力富集区域'\n",
  "}}\n\n",
  
  "## 主要发现\n",
  "{if (2 * pnorm(-abs(cd28_treg_ibd_rg$results$rg / cd28_treg_ibd_rg$results$rg_se)) < 0.05) {\n",
  "  sprintf('CD28+Treg细胞比例与IBD存在显著的遗传相关性 (rg = %.3f, p = %.3e)', \n",
  "          cd28_treg_ibd_rg$results$rg, \n",
  "          2 * pnorm(-abs(cd28_treg_ibd_rg$results$rg / cd28_treg_ibd_rg$results$rg_se)))\n",
  "} else {\n",
  "  '未发现CD28+Treg细胞比例与IBD之间存在显著的遗传相关性'\n",
  "}}\n\n",
  
  "{if (cd28_treg_h2$results$h2 > 0 && 2 * pnorm(-abs(cd28_treg_h2$results$h2 / cd28_treg_h2$results$h2_se)) < 0.05) {\n",
  "  sprintf('CD28+Treg细胞比例具有显著的遗传力 (h² = %.3f, p = %.3e)，表明遗传因素对该性状有重要影响', \n",
  "          cd28_treg_h2$results$h2, \n",
  "          2 * pnorm(-abs(cd28_treg_h2$results$h2 / cd28_treg_h2$results$h2_se)))\n",
  "} else {\n",
  "  'CD28+Treg细胞比例的遗传力估计不显著'\n",
  "}}\n\n",
  
  "## 局限性\n",
  "- 本分析基于模拟数据，实际应用中应使用真实GWAS数据和LD评分文件\n",
  "- 未考虑样本重叠对遗传相关性估计的影响\n",
  "- 功能注释文件为模拟生成，实际分析中应使用真实的功能注释数据\n",
  "- 未进行跨 ancestry 分析\n\n",
  
  "## 建议\n",
  "- 使用真实GWAS数据和完整的LD信息进行分析\n",
  "- 对遗传相关性显著的性状对进行共定位分析，探究共享的遗传机制\n",
  "- 使用多种方法验证遗传相关性结果\n",
  "- 结合功能基因组学数据，深入理解遗传相关性的生物学意义\n\n",
  
  "## 输出文件\n",
  "- `heritability_summary.csv`: 遗传力分析汇总结果\n",
  "- `genetic_correlation_summary.csv`: 遗传相关性分析汇总结果\n",
  "- `partitioned_heritability_summary.csv`: 分区遗传力分析汇总结果\n",
  "- `heritability_forest_plot.png`: 遗传力森林图\n",
  "- `genetic_correlation_forest_plot.png`: 遗传相关性森林图\n",
  "- `genetic_correlation_heatmap.png`: 遗传相关性热图\n",
  "- `partitioned_heritability_plot.png`: 分区遗传力图\n\n",
  
  "---\n",
  "*报告生成时间: {sys.time()}*"
)

# 保存报告
writeLines(
  report_content,
  here(output_dir, "ldsc_analysis_report.md")
)

# ----------------------
# 7. 完成分析
# ----------------------

cat("LD评分回归分析完成！\n")
cat("分析结果已保存到:", output_dir, "\n")
cat("主要输出文件包括:\n")
cat("- heritability_summary.csv\n")
cat("- genetic_correlation_summary.csv\n")
cat("- partitioned_heritability_summary.csv\n")
cat("- heritability_forest_plot.png\n")
cat("- genetic_correlation_forest_plot.png\n")
cat("- genetic_correlation_heatmap.png\n")
cat("- partitioned_heritability_plot.png\n")
cat("- ldsc_analysis_report.md\n")

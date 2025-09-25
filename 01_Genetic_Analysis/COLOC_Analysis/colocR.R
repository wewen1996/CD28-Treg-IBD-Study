#' CD28-Treg-IBD-Study: 共定位分析
#' 
#' 本脚本用于分析CD28+Treg细胞特征和IBD之间的共定位信号，
#' 以识别可能同时影响这两种表型的遗传变异。

# 加载必要的包
library(coloc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(here)
library(glue)
library(data.table)

# 设置输出目录
output_dir <- here("results", "COLOC_Analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------
# 1. 数据准备
# ----------------------

# 注意：在实际分析中，应从data目录加载真实数据
# 这里使用示例数据进行演示

set.seed(1234)
n_snps <- 1000

# 创建示例暴露数据（CD28+Treg细胞比例）
cd28_treg_data <- data.frame(
  SNP = paste0("rs", sample(100000:999999, n_snps)),
  CHR = sample(1:22, n_snps, replace = TRUE),
  POS = sample(1:100000000, n_snps),
  beta = rnorm(n_snps, mean = 0, sd = 0.1),
  se = runif(n_snps, min = 0.01, max = 0.05),
  pval = rbeta(n_snps, shape1 = 1, shape2 = 100)
)

# 创建示例结局数据（IBD）
ibd_data <- data.frame(
  SNP = cd28_treg_data$SNP,
  CHR = cd28_treg_data$CHR,
  POS = cd28_treg_data$POS,
  beta = rnorm(n_snps, mean = 0, sd = 0.1),
  se = runif(n_snps, min = 0.01, max = 0.05),
  pval = rbeta(n_snps, shape1 = 1, shape2 = 100)
)

# 为了模拟共定位信号，在特定区域设置显著关联
# 选择第6号染色体上的一个区域
cd28_treg_data$pval[cd28_treg_data$CHR == 6 & cd28_treg_data$POS > 30000000 & cd28_treg_data$POS < 32000000] <- 
  rbeta(sum(cd28_treg_data$CHR == 6 & cd28_treg_data$POS > 30000000 & cd28_treg_data$POS < 32000000), 
        shape1 = 0.1, shape2 = 10)

ibd_data$pval[ibd_data$CHR == 6 & ibd_data$POS > 30000000 & ibd_data$POS < 32000000] <- 
  rbeta(sum(ibd_data$CHR == 6 & ibd_data$POS > 30000000 & ibd_data$POS < 32000000), 
        shape1 = 0.1, shape2 = 10)

# 计算log10(p值)
cd28_treg_data$neg_log10_p <- -log10(cd28_treg_data$pval)
ibd_data$neg_log10_p <- -log10(ibd_data$pval)

# ----------------------
# 2. 全基因组关联分析可视化
# ----------------------

cat("生成全基因组关联分析曼哈顿图...\n")

# 2.1 CD28+Treg细胞特征的曼哈顿图
manhattan_cd28 <- cd28_treg_data %>%
  # 为每个染色体创建x轴位置
  group_by(CHR) %>%
  mutate(
    CHR_POS = POS,
    BPcum = cumsum(max(CHR_POS) - min(CHR_POS) + 1e6) - (max(CHR_POS) - min(CHR_POS) + 1e6)/2 + CHR_POS - min(CHR_POS)
  ) %>%
  ungroup() %>%
  # 创建曼哈顿图
  ggplot(aes(x = BPcum, y = neg_log10_p)) +
  # 按染色体着色
  geom_point(aes(color = as.factor(CHR)), alpha = 0.7, size = 1) +
  # 添加显著线
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(1e-5), color = "blue", linetype = "dashed") +
  # 设置x轴标签
  scale_x_continuous(
    labels = unique(cd28_treg_data$CHR),
    breaks = cd28_treg_data %>%
      group_by(CHR) %>%
      summarize(mean_BPcum = mean(BPcum)) %>%
      pull(mean_BPcum)
  ) +
  # 设置颜色
  scale_color_manual(values = rep(c("#2c3e50", "#3498db"), 11)) +
  # 添加标题和标签
  labs(
    title = "CD28+Treg细胞比例 GWAS",
    x = "染色体",
    y = "-log10(P值)",
    color = "染色体"
  ) +
  # 设置主题
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )

# 2.2 IBD的曼哈顿图
manhattan_ibd <- ibd_data %>%
  # 为每个染色体创建x轴位置
  group_by(CHR) %>%
  mutate(
    CHR_POS = POS,
    BPcum = cumsum(max(CHR_POS) - min(CHR_POS) + 1e6) - (max(CHR_POS) - min(CHR_POS) + 1e6)/2 + CHR_POS - min(CHR_POS)
  ) %>%
  ungroup() %>%
  # 创建曼哈顿图
  ggplot(aes(x = BPcum, y = neg_log10_p)) +
  # 按染色体着色
  geom_point(aes(color = as.factor(CHR)), alpha = 0.7, size = 1) +
  # 添加显著线
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(1e-5), color = "blue", linetype = "dashed") +
  # 设置x轴标签
  scale_x_continuous(
    labels = unique(ibd_data$CHR),
    breaks = ibd_data %>%
      group_by(CHR) %>%
      summarize(mean_BPcum = mean(BPcum)) %>%
      pull(mean_BPcum)
  ) +
  # 设置颜色
  scale_color_manual(values = rep(c("#2c3e50", "#3498db"), 11)) +
  # 添加标题和标签
  labs(
    title = "炎症性肠病(IBD) GWAS",
    x = "染色体",
    y = "-log10(P值)",
    color = "染色体"
  ) +
  # 设置主题
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )

# 合并两个曼哈顿图
combined_manhattan <- ggarrange(manhattan_cd28, manhattan_ibd, ncol = 1, nrow = 2, align = "v")

# 保存曼哈顿图
ggsave(
  filename = here(output_dir, "gwas_manhattan_plots.png"),
  plot = combined_manhattan,
  width = 14,
  height = 10,
  dpi = 300
)

# ----------------------
# 3. 识别显著关联区域
# ----------------------

cat("识别显著关联区域...\n")

# 定义显著阈值
significance_threshold <- 5e-8
suggestive_threshold <- 1e-5

# 识别CD28+Treg细胞特征的显著关联区域
cd28_significant <- cd28_treg_data %>%
  filter(pval < significance_threshold) %>%
  group_by(CHR) %>%
  summarize(
    start = min(POS),
    end = max(POS),
    n_snps = n()
  ) %>%
  filter(n_snps >= 1) %>%
  mutate(region = paste0("chr", CHR, ":", start, "-", end))

# 识别IBD的显著关联区域
ibd_significant <- ibd_data %>%
  filter(pval < significance_threshold) %>%
  group_by(CHR) %>%
  summarize(
    start = min(POS),
    end = max(POS),
    n_snps = n()
  ) %>%
  filter(n_snps >= 1) %>%
  mutate(region = paste0("chr", CHR, ":", start, "-", end))

# 打印显著区域信息
cat("CD28+Treg细胞特征的显著关联区域:\n")
print(cd28_significant)

cat("\nIBD的显著关联区域:\n")
print(ibd_significant)

# ----------------------
# 4. 共定位分析
# ----------------------

cat("\n开始共定位分析...\n")

# 4.1 准备共定位分析数据
prepare_coloc_data <- function(exposure_data, outcome_data, region_chr, region_start, region_end) {
  # 提取区域内的SNP
  exposure_region <- exposure_data %>%
    filter(CHR == region_chr, POS >= region_start, POS <= region_end) %>%
    select(SNP, beta, se, pval) %>%
    rename(beta.exposure = beta, se.exposure = se, pval.exposure = pval)
  
  outcome_region <- outcome_data %>%
    filter(CHR == region_chr, POS >= region_start, POS <= region_end) %>%
    select(SNP, beta, se, pval) %>%
    rename(beta.outcome = beta, se.outcome = se, pval.outcome = pval)
  
  # 合并数据
  coloc_data <- inner_join(exposure_region, outcome_region, by = "SNP")
  
  # 准备coloc输入格式
  coloc_input <- list(
    p1 = coloc_data$pval.exposure,
    p2 = coloc_data$pval.outcome,
    beta1 = coloc_data$beta.exposure,
    beta2 = coloc_data$beta.outcome,
    varbeta1 = (coloc_data$se.exposure)^2,
    varbeta2 = (coloc_data$se.outcome)^2,
    snp = coloc_data$SNP,
    type = "quant",  # 假设是数量性状
    MAF = rep(0.2, nrow(coloc_data))  # 假设MAF为0.2（实际分析中应使用真实MAF）
  )
  
  return(list(
    coloc_input = coloc_input,
    coloc_data = coloc_data
  ))
}

# 4.2 对每个显著区域进行共定位分析
coloc_results <- list()

# 如果有显著区域，则进行共定位分析
if (nrow(cd28_significant) > 0 && nrow(ibd_significant) > 0) {
  # 找出重叠的染色体
  common_chromosomes <- intersect(cd28_significant$CHR, ibd_significant$CHR)
  
  for (chr in common_chromosomes) {
    cat(sprintf("分析染色体 %d 上的共定位信号...\n", chr))
    
    # 获取该染色体上的所有显著区域
    cd28_regions <- cd28_significant %>% filter(CHR == chr)
    ibd_regions <- ibd_significant %>% filter(CHR == chr)
    
    # 对每个区域对进行分析
    for (i in 1:nrow(cd28_regions)) {
      for (j in 1:nrow(ibd_regions)) {
        # 检查区域是否重叠
        if (cd28_regions$start[i] <= ibd_regions$end[j] && cd28_regions$end[i] >= ibd_regions$start[j]) {
          # 定义重叠区域
          overlap_start <- max(cd28_regions$start[i], ibd_regions$start[j])
          overlap_end <- min(cd28_regions$end[i], ibd_regions$end[j])
          
          cat(sprintf("分析重叠区域 chr%s:%s-%s...\n", chr, overlap_start, overlap_end))
          
          # 准备共定位数据
          prepared_data <- prepare_coloc_data(
            exposure_data = cd28_treg_data,
            outcome_data = ibd_data,
            region_chr = chr,
            region_start = overlap_start,
            region_end = overlap_end
          )
          
          # 如果有足够的SNP进行分析
          if (nrow(prepared_data$coloc_data) >= 5) {  # 至少需要5个SNP
            # 执行共定位分析
            coloc_result <- coloc.abf(prepared_data$coloc_input)
            
            # 保存结果
            region_name <- paste0("chr", chr, ":", overlap_start, "-", overlap_end)
            coloc_results[[region_name]] <- list(
              region = region_name,
              chr = chr,
              start = overlap_start,
              end = overlap_end,
              coloc_result = coloc_result,
              coloc_data = prepared_data$coloc_data,
              n_snps = nrow(prepared_data$coloc_data)
            )
            
            # 打印结果摘要
            cat(sprintf("区域 %s 的共定位分析结果:\n", region_name))
            print(coloc_result$summary)
            cat("\n")
          } else {
            cat(sprintf("区域 chr%s:%s-%s 的SNP数量不足，跳过分析\n", chr, overlap_start, overlap_end))
          }
        }
      }
    }
  }
} else {
  cat("未发现显著关联区域，使用全基因组数据进行共定位分析...\n")
  
  # 如果没有显著区域，使用全基因组数据进行分析
  # 为了计算效率，这里只使用前1000个SNP
  prepared_data <- prepare_coloc_data(
    exposure_data = cd28_treg_data %>% slice(1:1000),
    outcome_data = ibd_data %>% slice(1:1000),
    region_chr = 1,
    region_start = 0,
    region_end = 1e12
  )
  
  # 执行共定位分析
  coloc_result <- coloc.abf(prepared_data$coloc_input)
  
  # 保存结果
  coloc_results[["genome_wide"]] <- list(
    region = "genome_wide",
    chr = "all",
    start = 0,
    end = 1e12,
    coloc_result = coloc_result,
    coloc_data = prepared_data$coloc_data,
    n_snps = nrow(prepared_data$coloc_data)
  )
  
  # 打印结果摘要
  cat("全基因组共定位分析结果:\n")
  print(coloc_result$summary)
  cat("\n")
}

# ----------------------
# 5. 共定位结果可视化
# ----------------------

cat("开始共定位结果可视化...\n")

# 5.1 共定位结果森林图
if (length(coloc_results) > 0) {
  # 提取共定位结果
  coloc_summary <- lapply(coloc_results, function(result) {
    data.frame(
      region = result$region,
      n_snps = result$n_snps,
      H0 = result$coloc_result$summary["H0"],
      H1 = result$coloc_result$summary["H1"],
      H2 = result$coloc_result$summary["H2"],
      H3 = result$coloc_result$summary["H3"],
      H4 = result$coloc_result$summary["H4"]
    )
  }) %>% bind_rows()
  
  # 创建共定位结果图
  coloc_plot <- coloc_summary %>%
    pivot_longer(cols = H0:H4, names_to = "hypothesis", values_to = "probability") %>%
    ggplot(aes(x = probability, y = region, fill = hypothesis)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(
      values = c("#2c3e50", "#3498db", "#e74c3c", "#2ecc71", "#f39c12"),
      labels = c("H0: 无关联", "H1: 仅暴露关联", "H2: 仅结局关联", "H3: 独立关联", "H4: 共定位")
    ) +
    labseller::scale_x_continuous(labels = scales::percent_format()) +
    labseller::labs(
      title = "共定位分析结果",
      subtitle = "各假设的后验概率",
      x = "后验概率",
      y = "基因组区域",
      fill = "假设"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
  
  # 保存共定位结果图
  ggsave(
    filename = here(output_dir, "coloc_results_plot.png"),
    plot = coloc_plot,
    width = 12,
    height = 8,
    dpi = 300
  )
}

# 5.2 共定位区域的精细定位图
# 对每个共定位区域创建精细定位图
for (region_name in names(coloc_results)) {
  result <- coloc_results[[region_name]]
  
  # 如果区域有足够的SNP，创建精细定位图
  if (result$n_snps >= 5) {
    # 准备数据
    plot_data <- result$coloc_data %>%
      mutate(
        neg_log10_p_exposure = -log10(pval.exposure),
        neg_log10_p_outcome = -log10(pval.outcome)
      ) %>%
      # 添加SNP位置信息
      left_join(
        cd28_treg_data %>% select(SNP, POS),
        by = "SNP"
      ) %>%
      arrange(POS)
    
    # 创建精细定位图
    finemapping_plot <- ggplot(plot_data, aes(x = POS/1e6)) +
      # CD28+Treg细胞特征的p值
      geom_point(aes(y = neg_log10_p_exposure, color = "CD28+Treg细胞比例"), size = 2) +
      # IBD的p值
      geom_point(aes(y = neg_log10_p_outcome, color = "IBD"), size = 2) +
      # 添加显著线
      geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed") +
      geom_hline(yintercept = -log10(1e-5), color = "blue", linetype = "dashed") +
      # 设置颜色
      scale_color_manual(values = c("#3498db", "#e74c3c")) +
      # 添加标题和标签
      labs(
        title = paste("区域", region_name, "的精细定位图"),
        subtitle = sprintf("共定位概率 H4 = %.3f", result$coloc_result$summary["H4"]),
        x = "位置 (Mb)",
        y = "-log10(P值)",
        color = "表型"
      ) +
      # 设置主题
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
    
    # 保存精细定位图
    ggsave(
      filename = here(output_dir, paste0(gsub(":", "_", gsub("-", "_", region_name)), "_finemapping.png")),
      plot = finemapping_plot,
      width = 12,
      height = 6,
      dpi = 300
    )
  }
}

# ----------------------
# 6. 结果汇总和保存
# ----------------------

cat("\n开始结果汇总...\n")

# 6.1 汇总共定位结果
if (length(coloc_results) > 0) {
  coloc_summary_table <- lapply(coloc_results, function(result) {
    data.frame(
      region = result$region,
      chr = result$chr,
      start = result$start,
      end = result$end,
      n_snps = result$n_snps,
      H0 = result$coloc_result$summary["H0"],
      H1 = result$coloc_result$summary["H1"],
      H2 = result$coloc_result$summary["H2"],
      H3 = result$coloc_result$summary["H3"],
      H4 = result$coloc_result$summary["H4"],
      coloc_evidence = ifelse(result$coloc_result$summary["H4"] > 0.8, "强证据",
                             ifelse(result$coloc_result$summary["H4"] > 0.5, "中等证据",
                                    ifelse(result$coloc_result$summary["H4"] > 0.2, "弱证据", "无证据")))
    )
  }) %>% bind_rows()
  
  # 保存汇总结果
  write.csv(
    coloc_summary_table,
    here(output_dir, "coloc_summary_results.csv"),
    row.names = FALSE
  )
  
  # 打印汇总结果
  cat("共定位分析汇总结果:\n")
  print(coloc_summary_table %>% select(region, n_snps, H4, coloc_evidence))
}

# 6.2 保存显著共定位区域的详细信息
if (length(coloc_results) > 0) {
  # 找出有强证据或中等证据的共定位区域
  significant_regions <- coloc_summary_table %>%
    filter(coloc_evidence %in% c("强证据", "中等证据")) %>%
    pull(region)
  
  if (length(significant_regions) > 0) {
    for (region in significant_regions) {
      result <- coloc_results[[region]]
      
      # 获取该区域的详细SNP信息
      snp_details <- result$coloc_data %>%
        left_join(
          cd28_treg_data %>% select(SNP, CHR, POS, beta.exposure = beta, se.exposure = se, pval.exposure = pval),
          by = "SNP"
        ) %>%
        left_join(
          ibd_data %>% select(SNP, beta.outcome = beta, se.outcome = se, pval.outcome = pval),
          by = "SNP"
        ) %>%
        mutate(
          OR.exposure = exp(beta.exposure),
          OR.outcome = exp(beta.outcome)
        ) %>%
        select(SNP, CHR, POS, beta.exposure, se.exposure, pval.exposure, OR.exposure,
               beta.outcome, se.outcome, pval.outcome, OR.outcome) %>%
        arrange(pval.exposure)
      
      # 保存SNP详细信息
      write.csv(
        snp_details,
        here(output_dir, paste0(gsub(":", "_", gsub("-", "_", region)), "_snp_details.csv")),
        row.names = FALSE
      )
    }
  }
}

# ----------------------
# 7. 生成分析报告
# ----------------------

cat("\n生成分析报告...\n")

# 创建分析报告
report_content <- glue(
  "# CD28-Treg-IBD-Study: 共定位分析报告\n\n",
  
  "## 分析概述\n",
  "- **分析日期**: {Sys.Date()}\n",
  "- **暴露因素**: CD28+Treg细胞比例\n",
  "- **结局因素**: 炎症性肠病(IBD)\n",
  "- **分析方法**: 共定位分析 (coloc)\n\n",
  
  "## 数据概况\n",
  "- **暴露数据SNP数量**: {nrow(cd28_treg_data)}\n",
  "- **结局数据SNP数量**: {nrow(ibd_data)}\n",
  "- **重叠SNP数量**: {length(intersect(cd28_treg_data$SNP, ibd_data$SNP))}\n\n",
  
  "## 显著关联区域\n",
  "### CD28+Treg细胞特征的显著关联区域\n",
  "{if (nrow(cd28_significant) > 0) {\n",
  "  paste(sprintf('- %s: %d个SNP', cd28_significant$region, cd28_significant$n_snps), collapse = '\\n')\n",
  "} else {\n",
  "  '未发现显著关联区域 (p < 5e-8)'\n",
  "}}\n\n",
  
  "### IBD的显著关联区域\n",
  "{if (nrow(ibd_significant) > 0) {\n",
  "  paste(sprintf('- %s: %d个SNP', ibd_significant$region, ibd_significant$n_snps), collapse = '\\n')\n",
  "} else {\n",
  "  '未发现显著关联区域 (p < 5e-8)'\n",
  "}}\n\n",
  
  "## 共定位分析结果\n",
  "{if (length(coloc_results) > 0) {\n",
  "  paste(sprintf('- 分析了%d个基因组区域', length(coloc_results)), collapse = '\\n')\n",
  "} else {\n",
  "  '未进行共定位分析'\n",
  "}}\n\n",
  
  "### 主要发现\n",
  "{if (length(coloc_results) > 0) {\n",
  "  significant_regions <- coloc_summary_table %>% filter(coloc_evidence %in% c('强证据', '中等证据'))\n",
  "  if (nrow(significant_regions) > 0) {\n",
  "    paste(c('发现以下区域存在共定位证据:', \n",
  "            sprintf('- %s: %s (H4 = %.3f)', \n",
  "                    significant_regions$region, \n",
  "                    significant_regions$coloc_evidence, \n",
  "                    significant_regions$H4)), collapse = '\\n')\n",
  "  } else {\n",
  "    '未发现有强证据或中等证据支持的共定位区域'\n",
  "  }\n",
  "} else {\n",
  "  '无共定位分析结果'\n",
  "}}\n\n",
  
  "### 共定位假设解释\n",
  "- **H0**: 该区域与暴露和结局均无关联\n",
  "- **H1**: 该区域仅与暴露有关联\n",
  "- **H2**: 该区域仅与结局有关联\n",
  "- **H3**: 该区域与暴露和结局均有关联，但由不同的因果变异引起\n",
  "- **H4**: 该区域与暴露和结局均有关联，且由相同的因果变异引起（共定位）\n\n",
  
  "## 敏感性分析\n",
  "- 使用默认的先验概率设置（P1=1e-4, P2=1e-4, P3=1e-5, P4=1e-5）\n",
  "- 假设所有SNP的次要等位基因频率(MAF)为0.2（实际分析中应使用真实MAF）\n",
  "- 假设性状为数量性状（quant）\n\n",
  
  "## 局限性\n",
  "- 本分析基于示例数据，实际应用中应使用真实GWAS数据\n",
  "- 共定位分析假设每个区域只有一个因果变异\n",
  "- 未考虑连锁不平衡(LD)结构对结果的影响\n",
  "- 先验概率的选择可能影响结果解释\n\n",
  
  "## 建议\n",
  "- 使用真实GWAS数据和完整的LD信息进行分析\n",
  "- 对有共定位证据的区域进行功能注释，探究可能的机制\n",
  "- 使用其他方法（如SMR）验证共定位结果\n",
  "- 在实验模型中验证关键发现\n\n",
  
  "## 输出文件\n",
  "- `gwas_manhattan_plots.png`: GWAS曼哈顿图\n",
  "- `coloc_results_plot.png`: 共定位分析结果图\n",
  "- `coloc_summary_results.csv`: 共定位分析汇总结果\n",
  "- 各共定位区域的精细定位图\n",
  "- 显著共定位区域的SNP详细信息\n\n",
  
  "---\n",
  "*报告生成时间: {Sys.time()}*"
)

# 保存报告
writeLines(
  report_content,
  here(output_dir, "coloc_analysis_report.md")
)

# ----------------------
# 8. 完成分析
# ----------------------

cat("共定位分析完成！\n")
cat("分析结果已保存到:", output_dir, "\n")
cat("主要输出文件包括:\n")
cat("- gwas_manhattan_plots.png\n")
cat("- coloc_results_plot.png\n")
cat("- coloc_summary_results.csv\n")
cat("- coloc_analysis_report.md\n")

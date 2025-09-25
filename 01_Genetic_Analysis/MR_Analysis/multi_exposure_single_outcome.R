#' CD28-Treg-IBD-Study: 多暴露-单结局孟德尔随机化分析
#' 
#' 本脚本用于分析多个暴露因素（CD28+Treg细胞特征）对单一结局（IBD）的因果效应。
#' 采用多种MR方法，并进行异质性和多效性检验，以全面评估因果关系。

# 加载必要的包
library(TwoSampleMR)
library(MendelianRandomization)
library(mr.raps)
library(RadialMR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(here)
library(glue)

# 加载自定义工具函数
source(here("01_Genetic_Analysis", "MR_Analysis", "utility_functions.R"))

# 设置输出目录
output_dir <- here("results", "MR_Analysis", "multi_exposure_single_outcome")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------
# 1. 数据准备
# ----------------------

# 注意：在实际分析中，应从data目录加载真实数据
# 这里使用示例数据进行演示

# 创建示例暴露数据（CD28+Treg细胞特征）
set.seed(1234)
n_snps <- 50

# 暴露因素1：CD28+Treg细胞比例
cd28_treg_prop <- data.frame(
  SNP = paste0("rs", sample(100000:999999, n_snps)),
  beta = rnorm(n_snps, mean = 0.1, sd = 0.05),
  se = runif(n_snps, min = 0.01, max = 0.03),
  Effect_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  Other_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  pval = rbeta(n_snps, shape1 = 0.5, shape2 = 5)  # 生成小p值
) %>%
  filter(pval < 5e-8)  # 筛选显著SNP

# 暴露因素2：CD28+Treg细胞活化状态
cd28_treg_activation <- data.frame(
  SNP = paste0("rs", sample(100000:999999, n_snps)),
  beta = rnorm(n_snps, mean = 0.08, sd = 0.04),
  se = runif(n_snps, min = 0.01, max = 0.03),
  Effect_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  Other_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  pval = rbeta(n_snps, shape1 = 0.5, shape2 = 5)
) %>%
  filter(pval < 5e-8)

# 暴露因素3：CD28+Treg细胞功能评分
cd28_treg_function <- data.frame(
  SNP = paste0("rs", sample(100000:999999, n_snps)),
  beta = rnorm(n_snps, mean = 0.12, sd = 0.06),
  se = runif(n_snps, min = 0.01, max = 0.03),
  Effect_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  Other_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  pval = rbeta(n_snps, shape1 = 0.5, shape2 = 5)
) %>%
  filter(pval < 5e-8)

# 创建示例结局数据（IBD）
ibd_data <- data.frame(
  SNP = c(cd28_treg_prop$SNP, cd28_treg_activation$SNP, cd28_treg_function$SNP),
  beta = rnorm(length(c(cd28_treg_prop$SNP, cd28_treg_activation$SNP, cd28_treg_function$SNP)), 
               mean = 0.05, sd = 0.02),
  se = runif(length(c(cd28_treg_prop$SNP, cd28_treg_activation$SNP, cd28_treg_function$SNP)), 
             min = 0.005, max = 0.015),
  Effect_allele = sample(c("A", "T", "C", "G"), 
                         length(c(cd28_treg_prop$SNP, cd28_treg_activation$SNP, cd28_treg_function$SNP)), 
                         replace = TRUE),
  Other_allele = sample(c("A", "T", "C", "G"), 
                        length(c(cd28_treg_prop$SNP, cd28_treg_activation$SNP, cd28_treg_function$SNP)), 
                        replace = TRUE),
  pval = rbeta(length(c(cd28_treg_prop$SNP, cd28_treg_activation$SNP, cd28_treg_function$SNP)), 
               shape1 = 1, shape2 = 10)
)

# 创建暴露数据列表和名称列表
exposure_data_list <- list(
  cd28_treg_prop = cd28_treg_prop,
  cd28_treg_activation = cd28_treg_activation,
  cd28_treg_function = cd28_treg_function
)

exposure_names <- c(
  "CD28+Treg细胞比例",
  "CD28+Treg细胞活化状态", 
  "CD28+Treg细胞功能评分"
)

outcome_name <- "炎症性肠病(IBD)"

# ----------------------
# 2. 工具变量选择
# ----------------------

cat("开始工具变量选择...\n")

# 对每个暴露因素选择独立的工具变量
selected_instruments <- lapply(exposure_data_list, function(data) {
  select_instruments(data, p_threshold = 5e-8, clump_r2 = 0.001, clump_kb = 10000)
})

# 打印每个暴露因素的工具变量数量
cat("工具变量选择结果:\n")
for (i in seq_along(selected_instruments)) {
  cat(sprintf("%s: %d个SNP\n", exposure_names[i], nrow(selected_instruments[[i]])))
}

# ----------------------
# 3. 多暴露MR分析
# ----------------------

cat("\n开始多暴露MR分析...\n")

# 执行多暴露MR分析
mr_results <- perform_multiple_exposure_mr(
  exposure_data_list = selected_instruments,
  outcome_data = ibd_data,
  exposure_names = exposure_names,
  outcome_name = outcome_name
)

# 打印主要结果
cat("多暴露MR分析主要结果:\n")
print(mr_results$mr_results)

# ----------------------
# 4. 单暴露MR分析（作为补充）
# ----------------------

cat("\n开始单暴露MR分析...\n")

# 对每个暴露因素单独进行MR分析
single_exposure_results <- list()

for (i in seq_along(selected_instruments)) {
  cat(sprintf("分析 %s -> %s...\n", exposure_names[i], outcome_name))
  
  # 格式化数据
  formatted_data <- format_mr_data(
    exposure_data = selected_instruments[[i]],
    outcome_data = ibd_data,
    exposure_name = exposure_names[i],
    outcome_name = outcome_name
  )
  
  # 执行MR分析
  mr_result <- perform_mr_analysis(
    mr_data = formatted_data$mr_data,
    methods = c("mr_wald_ratio", "mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
  )
  
  # 保存结果
  single_exposure_results[[i]] <- list(
    formatted_data = formatted_data,
    mr_result = mr_result
  )
  
  # 打印报告
  print_mr_report(mr_result, exposure_names[i], outcome_name)
}

names(single_exposure_results) <- exposure_names

# ----------------------
# 5. 结果可视化
# ----------------------

cat("\n开始结果可视化...\n")

# 5.1 多暴露MR结果森林图
multi_exposure_forest <- mr_results$mr_results %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    exposure = exposure_names,
    ci_lower = b - 1.96 * se,
    ci_upper = b + 1.96 * se
  ) %>%
  ggplot(aes(x = b, y = exposure)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
  geom_text(aes(label = sprintf("%.3f (%.3f, %.3f)", b, ci_lower, ci_upper)), 
            hjust = -0.1, size = 3.5) +
  labs(
    title = "多暴露MR分析结果",
    subtitle = "CD28+Treg细胞特征对IBD的因果效应",
    x = "效应量 (95% CI)",
    y = "暴露因素"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.major.x = element_line(linetype = "dashed", alpha = 0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.y = element_blank()
  ) +
  expand_limits(x = c(min(multi_exposure_forest$data$ci_lower) * 1.2, max(multi_exposure_forest$data$ci_upper) * 1.2))

# 保存多暴露森林图
ggsave(
  filename = here(output_dir, "multi_exposure_forest_plot.png"),
  plot = multi_exposure_forest,
  width = 10,
  height = 6,
  dpi = 300
)

# 5.2 单暴露MR结果可视化
for (i in seq_along(single_exposure_results)) {
  exposure <- exposure_names[i]
  result <- single_exposure_results[[i]]
  
  # 创建综合可视化图
  p_forest <- plot_mr_forest(
    mr_results = result$mr_result,
    title = sprintf("%s -> %s", exposure, outcome_name),
    xlab = "效应量 (95% CI)"
  )
  
  p_scatter <- plot_mr_scatter(
    mr_data = result$formatted_data$mr_data,
    mr_results = result$mr_result,
    title = sprintf("%s -> %s", exposure, outcome_name)
  )
  
  p_funnel <- plot_mr_funnel(
    mr_data = result$formatted_data$mr_data,
    title = sprintf("%s -> %s", exposure, outcome_name)
  )
  
  # 组合图形
  combined_plot <- (p_forest + p_scatter) / p_funnel +
    plot_annotation(
      title = sprintf("%s对%s的MR分析结果", exposure, outcome_name),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    )
  
  # 保存图形
  ggsave(
    filename = here(output_dir, sprintf("%s_mr_analysis.png", gsub(" ", "_", gsub("[^a-zA-Z0-9_ ]", "", exposure)))),
    plot = combined_plot,
    width = 12,
    height = 10,
    dpi = 300
  )
}

# ----------------------
# 6. 结果汇总和保存
# ----------------------

cat("\n开始结果汇总...\n")

# 6.1 多暴露MR结果汇总
multi_exposure_summary <- mr_results$mr_results %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    exposure = exposure_names,
    OR = exp(b),
    OR_lower = exp(b - 1.96 * se),
    OR_upper = exp(b + 1.96 * se),
    outcome = outcome_name
  ) %>%
  select(exposure, outcome, method, nsnp, b, se, OR, OR_lower, OR_upper, pval)

# 保存多暴露结果
write.csv(
  multi_exposure_summary,
  here(output_dir, "multi_exposure_mr_results.csv"),
  row.names = FALSE
)

# 6.2 单暴露MR结果汇总
single_exposure_summary <- lapply(seq_along(single_exposure_results), function(i) {
  exposure <- exposure_names[i]
  result <- single_exposure_results[[i]]
  
  summarize_mr_results(
    mr_results_list = result$mr_result,
    exposure_name = exposure,
    outcome_name = outcome_name
  )
}) %>% bind_rows()

# 保存单暴露结果
write.csv(
  single_exposure_summary,
  here(output_dir, "single_exposure_mr_results.csv"),
  row.names = FALSE
)

# 6.3 保存工具变量信息
tool_vars_info <- lapply(seq_along(selected_instruments), function(i) {
  exposure <- exposure_names[i]
  data <- selected_instruments[[i]] %>%
    mutate(exposure = exposure) %>%
    select(exposure, everything())
  
  return(data)
}) %>% bind_rows()

write.csv(
  tool_vars_info,
  here(output_dir, "instrumental_variables_info.csv"),
  row.names = FALSE
)

# ----------------------
# 7. 生成分析报告
# ----------------------

cat("\n生成分析报告...\n")

# 创建分析报告
report_content <- glue(
  "# CD28-Treg-IBD-Study: 多暴露-单结局孟德尔随机化分析报告\n\n",
  
  "## 分析概述\n",
  "- **分析日期**: {Sys.Date()}\n",
  "- **暴露因素**: {paste(exposure_names, collapse = ', ')}\n",
  "- **结局因素**: {outcome_name}\n",
  "- **分析方法**: 多暴露MR分析，单暴露MR分析（IVW, MR Egger, 加权中位数, 加权众数）\n\n",
  
  "## 数据概况\n",
  "### 工具变量选择结果\n",
  "| 暴露因素 | 工具变量数量 |\n",
  "|----------|--------------|\n",
  "{paste(sprintf('| %s | %d |', exposure_names, sapply(selected_instruments, nrow)), collapse = '\\n')}\n\n",
  
  "## 主要发现\n",
  "### 多暴露MR分析结果\n",
  "| 暴露因素 | 效应量 | 标准误 | OR (95% CI) | P值 |\n",
  "|----------|--------|--------|-------------|-----|\n",
  "{paste(sprintf('| %s | %.3f | %.3f | %.3f (%.3f-%.3f) | %.3f |', 
         multi_exposure_summary$exposure,
         multi_exposure_summary$b,
         multi_exposure_summary$se,
         multi_exposure_summary$OR,
         multi_exposure_summary$OR_lower,
         multi_exposure_summary$OR_upper,
         multi_exposure_summary$pval), collapse = '\\n')}\n\n",
  
  "### 关键结论\n",
  "{if (any(multi_exposure_summary$pval < 0.05)) {\n",
  "  significant_exposures <- multi_exposure_summary %>% filter(pval < 0.05) %>% pull(exposure)\n",
  "  paste(sprintf('- %s与%s存在显著因果关联 (p < 0.05)', significant_exposures, outcome_name), collapse = '\\n')\n",
  "} else {\n",
  "  sprintf('- 未发现任何暴露因素与%s存在显著因果关联 (p ≥ 0.05)', outcome_name)\n",
  "}}\n\n",
  
  "## 敏感性分析\n",
  "### 异质性检验\n",
  "{paste(lapply(seq_along(single_exposure_results), function(i) {\n",
  "  exposure <- exposure_names[i]\n",
  "  hetero <- single_exposure_results[[i]]$mr_result$heterogeneity %>% filter(method == 'Inverse variance weighted')\n",
  "  if (hetero$Q_pval < 0.05) {\n",
  "    sprintf('- %s: 存在显著异质性 (Q = %.3f, p = %.3f, I² = %.3f)', exposure, hetero$Q, hetero$Q_pval, hetero$I2)\n",
  "  } else {\n",
  "    sprintf('- %s: 未检测到显著异质性 (Q = %.3f, p = %.3f, I² = %.3f)', exposure, hetero$Q, hetero$Q_pval, hetero$I2)\n",
  "  }\n",
  "}), collapse = '\\n')}\n\n",
  
  "### 多效性检验\n",
  "{paste(lapply(seq_along(single_exposure_results), function(i) {\n",
  "  exposure <- exposure_names[i]\n",
  "  pleio <- single_exposure_results[[i]]$mr_result$pleiotropy\n",
  "  if (pleio$egger_intercept_pval < 0.05) {\n",
  "    sprintf('- %s: MR Egger截距检验提示可能存在水平多效性 (intercept = %.3f, p = %.3f)', exposure, pleio$egger_intercept, pleio$egger_intercept_pval)\n",
  "  } else {\n",
  "    sprintf('- %s: MR Egger截距检验未提示水平多效性 (intercept = %.3f, p = %.3f)', exposure, pleio$egger_intercept, pleio$egger_intercept_pval)\n",
  "  }\n",
  "}), collapse = '\\n')}\n\n",
  
  "## 局限性\n",
  "- 本分析基于示例数据，实际应用中应使用真实GWAS数据\n",
  "- 工具变量的有效性假设（相关性、独立性、排他性）需要进一步验证\n",
  "- 未考虑潜在的基因-环境交互作用\n",
  "- 多暴露之间可能存在相关性，需要进一步的多变量分析\n\n",
  
  "## 建议\n",
  "- 使用更大样本量的GWAS数据验证结果\n",
  "- 进行多变量MR分析，控制暴露因素之间的相关性\n",
  "- 结合功能基因组学数据，探究因果关系的分子机制\n",
  "- 在实验模型中验证关键发现\n\n",
  
  "## 分析流程\n",
  "```mermaid\n",
  "graph TD\n",
  "    A[数据准备] --> B[工具变量选择]\n",
  "    B --> C[多暴露MR分析]\n",
  "    B --> D[单暴露MR分析]\n",
  "    C --> E[结果可视化]\n",
  "    D --> E\n",
  "    E --> F[结果汇总和报告]\n",
  "```\n\n",
  
  "## 输出文件\n",
  "- `multi_exposure_forest_plot.png`: 多暴露MR分析森林图\n",
  "- `multi_exposure_mr_results.csv`: 多暴露MR分析结果表\n",
  "- `single_exposure_mr_results.csv`: 单暴露MR分析结果表\n",
  "- `instrumental_variables_info.csv`: 工具变量信息表\n",
  "- 各暴露因素的MR分析综合图\n\n",
  
  "---\n",
  "*报告生成时间: {Sys.time()}*"
)

# 保存报告
writeLines(
  report_content,
  here(output_dir, "multi_exposure_single_outcome_mr_report.md")
)

# ----------------------
# 8. 完成分析
# ----------------------

cat("多暴露-单结局MR分析完成！\n")
cat("分析结果已保存到:", output_dir, "\n")
cat("主要输出文件包括:\n")
cat("- multi_exposure_forest_plot.png\n")
cat("- multi_exposure_mr_results.csv\n")
cat("- single_exposure_mr_results.csv\n")
cat("- instrumental_variables_info.csv\n")
cat("- multi_exposure_single_outcome_mr_report.md\n")

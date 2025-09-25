#' CD28-Treg-IBD-Study: 单暴露-多结局孟德尔随机化分析
#' 
#' 本脚本用于分析单一暴露因素（CD28+Treg细胞特征）对多个结局（IBD及其亚型）的因果效应。
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
output_dir <- here("results", "MR_Analysis", "single_exposure_multi_outcome")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------
# 1. 数据准备
# ----------------------

# 注意：在实际分析中，应从data目录加载真实数据
# 这里使用示例数据进行演示

# 创建示例暴露数据（CD28+Treg细胞特征）
set.seed(1234)
n_snps <- 50

# 暴露因素：CD28+Treg细胞比例
cd28_treg_prop <- data.frame(
  SNP = paste0("rs", sample(100000:999999, n_snps)),
  beta = rnorm(n_snps, mean = 0.1, sd = 0.05),
  se = runif(n_snps, min = 0.01, max = 0.03),
  Effect_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  Other_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  pval = rbeta(n_snps, shape1 = 0.5, shape2 = 5)  # 生成小p值
) %>%
  filter(pval < 5e-8)  # 筛选显著SNP

# 创建示例结局数据（IBD及其亚型）
# 结局1：总体IBD
ibd_data <- data.frame(
  SNP = cd28_treg_prop$SNP,
  beta = rnorm(n_snps, mean = 0.05, sd = 0.02),
  se = runif(n_snps, min = 0.005, max = 0.015),
  Effect_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  Other_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  pval = rbeta(n_snps, shape1 = 1, shape2 = 10)
)

# 结局2：克罗恩病(CD)
cd_data <- data.frame(
  SNP = cd28_treg_prop$SNP,
  beta = rnorm(n_snps, mean = 0.06, sd = 0.025),
  se = runif(n_snps, min = 0.005, max = 0.015),
  Effect_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  Other_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  pval = rbeta(n_snps, shape1 = 1, shape2 = 10)
)

# 结局3：溃疡性结肠炎(UC)
uc_data <- data.frame(
  SNP = cd28_treg_prop$SNP,
  beta = rnorm(n_snps, mean = 0.04, sd = 0.015),
  se = runif(n_snps, min = 0.005, max = 0.015),
  Effect_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  Other_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  pval = rbeta(n_snps, shape1 = 1, shape2 = 10)
)

# 结局4：IBD相关并发症
ibd_complication_data <- data.frame(
  SNP = cd28_treg_prop$SNP,
  beta = rnorm(n_snps, mean = 0.07, sd = 0.03),
  se = runif(n_snps, min = 0.005, max = 0.015),
  Effect_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  Other_allele = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  pval = rbeta(n_snps, shape1 = 1, shape2 = 10)
)

# 创建结局数据列表和名称列表
outcome_data_list <- list(
  ibd = ibd_data,
  cd = cd_data,
  uc = uc_data,
  ibd_complication = ibd_complication_data
)

outcome_names <- c(
  "总体炎症性肠病(IBD)",
  "克罗恩病(CD)",
  "溃疡性结肠炎(UC)",
  "IBD相关并发症"
)

exposure_name <- "CD28+Treg细胞比例"

# ----------------------
# 2. 工具变量选择
# ----------------------

cat("开始工具变量选择...\n")

# 选择独立的工具变量
selected_instruments <- select_instruments(
  cd28_treg_prop,
  p_threshold = 5e-8,
  clump_r2 = 0.001,
  clump_kb = 10000
)

# 打印工具变量数量
cat(sprintf("为%s选择了%d个工具变量\n", exposure_name, nrow(selected_instruments)))

# ----------------------
# 3. 单暴露-多结局MR分析
# ----------------------

cat("\n开始单暴露-多结局MR分析...\n")

# 对每个结局进行MR分析
mr_results_list <- list()

for (i in seq_along(outcome_data_list)) {
  cat(sprintf("分析 %s -> %s...\n", exposure_name, outcome_names[i]))
  
  # 格式化数据
  formatted_data <- format_mr_data(
    exposure_data = selected_instruments,
    outcome_data = outcome_data_list[[i]],
    exposure_name = exposure_name,
    outcome_name = outcome_names[i]
  )
  
  # 执行MR分析
  mr_result <- perform_mr_analysis(
    mr_data = formatted_data$mr_data,
    methods = c("mr_wald_ratio", "mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
  )
  
  # 执行径向MR分析（检测异常值）
  radial_result <- perform_radial_mr(formatted_data$mr_data)
  
  # 保存结果
  mr_results_list[[i]] <- list(
    formatted_data = formatted_data,
    mr_result = mr_result,
    radial_result = radial_result
  )
  
  # 打印报告
  print_mr_report(mr_result, exposure_name, outcome_names[i])
}

names(mr_results_list) <- outcome_names

# ----------------------
# 4. 结果可视化
# ----------------------

cat("\n开始结果可视化...\n")

# 4.1 多结局森林图
# 提取所有结局的IVW结果
ivw_results <- lapply(seq_along(mr_results_list), function(i) {
  outcome <- outcome_names[i]
  result <- mr_results_list[[i]]$mr_result$mr_results %>%
    filter(method == "Inverse variance weighted") %>%
    mutate(
      outcome = outcome,
      ci_lower = b - 1.96 * se,
      ci_upper = b + 1.96 * se
    )
  
  return(result)
}) %>% bind_rows()

# 创建多结局森林图
multi_outcome_forest <- ivw_results %>%
  ggplot(aes(x = b, y = outcome)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
  geom_text(aes(label = sprintf("%.3f (%.3f, %.3f)", b, ci_lower, ci_upper)), 
            hjust = -0.1, size = 3.5) +
  labs(
    title = "单暴露-多结局MR分析结果",
    subtitle = sprintf("%s对不同IBD表型的因果效应", exposure_name),
    x = "效应量 (95% CI)",
    y = "结局"
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
  expand_limits(x = c(min(ivw_results$ci_lower) * 1.2, max(ivw_results$ci_upper) * 1.2))

# 保存多结局森林图
ggsave(
  filename = here(output_dir, "multi_outcome_forest_plot.png"),
  plot = multi_outcome_forest,
  width = 10,
  height = 6,
  dpi = 300
)

# 4.2 OR值比较图
or_results <- ivw_results %>%
  mutate(
    OR = exp(b),
    OR_lower = exp(ci_lower),
    OR_upper = exp(ci_upper)
  )

or_comparison_plot <- or_results %>%
  ggplot(aes(x = OR, y = outcome)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper), height = 0.2) +
  geom_text(aes(label = sprintf("%.3f (%.3f, %.3f)", OR, OR_lower, OR_upper)), 
            hjust = -0.1, size = 3.5) +
  labs(
    title = "单暴露-多结局MR分析OR值比较",
    subtitle = sprintf("%s对不同IBD表型的因果效应", exposure_name),
    x = "比值比 (OR, 95% CI)",
    y = "结局"
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
  expand_limits(x = c(0.9, max(or_results$OR_upper) * 1.2))

# 保存OR值比较图
ggsave(
  filename = here(output_dir, "or_comparison_plot.png"),
  plot = or_comparison_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# 4.3 每个结局的详细分析图
for (i in seq_along(mr_results_list)) {
  outcome <- outcome_names[i]
  result <- mr_results_list[[i]]
  
  # 创建综合可视化图
  p_forest <- plot_mr_forest(
    mr_results = result$mr_result,
    title = sprintf("%s -> %s", exposure_name, outcome),
    xlab = "效应量 (95% CI)"
  )
  
  p_scatter <- plot_mr_scatter(
    mr_data = result$formatted_data$mr_data,
    mr_results = result$mr_result,
    title = sprintf("%s -> %s", exposure_name, outcome)
  )
  
  p_funnel <- plot_mr_funnel(
    mr_data = result$formatted_data$mr_data,
    title = sprintf("%s -> %s", exposure_name, outcome)
  )
  
  p_radial <- plot_radial_mr(
    radial_result = result$radial_result,
    title = sprintf("%s -> %s", exposure_name, outcome)
  )
  
  # 组合图形
  combined_plot <- (p_forest + p_scatter) / (p_funnel + p_radial) +
    plot_annotation(
      title = sprintf("%s对%s的MR分析结果", exposure_name, outcome),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    )
  
  # 保存图形
  ggsave(
    filename = here(output_dir, sprintf("%s_mr_analysis.png", gsub(" ", "_", gsub("[^a-zA-Z0-9_ ]", "", outcome)))),
    plot = combined_plot,
    width = 14,
    height = 10,
    dpi = 300
  )
}

# ----------------------
# 5. 结果汇总和保存
# ----------------------

cat("\n开始结果汇总...\n")

# 5.1 所有结局的MR结果汇总
all_results_summary <- lapply(seq_along(mr_results_list), function(i) {
  outcome <- outcome_names[i]
  result <- mr_results_list[[i]]
  
  summarize_mr_results(
    mr_results_list = result$mr_result,
    exposure_name = exposure_name,
    outcome_name = outcome
  )
}) %>% bind_rows()

# 保存汇总结果
write.csv(
  all_results_summary,
  here(output_dir, "all_outcomes_mr_results.csv"),
  row.names = FALSE
)

# 5.2 保存工具变量信息
write.csv(
  selected_instruments %>% mutate(exposure = exposure_name),
  here(output_dir, "instrumental_variables_info.csv"),
  row.names = FALSE
)

# 5.3 保存径向MR分析结果
radial_results_summary <- lapply(seq_along(mr_results_list), function(i) {
  outcome <- outcome_names[i]
  result <- mr_results_list[[i]]$radial_result
  
  data.frame(
    outcome = outcome,
    method = "Radial MR",
    b = result$coef,
    se = result$se,
    pval = result$pval,
    Q = result$Q,
    Q_df = result$Q_df,
    Q_pval = result$Q_pval,
    outliers_identified = ifelse(length(result$outliers) > 0, "Yes", "No"),
    n_outliers = length(result$outliers),
    outliers = ifelse(length(result$outliers) > 0, paste(result$outliers, collapse = ", "), NA)
  )
}) %>% bind_rows()

write.csv(
  radial_results_summary,
  here(output_dir, "radial_mr_results.csv"),
  row.names = FALSE
)

# ----------------------
# 6. 生成分析报告
# ----------------------

cat("\n生成分析报告...\n")

# 创建分析报告
report_content <- glue(
  "# CD28-Treg-IBD-Study: 单暴露-多结局孟德尔随机化分析报告\n\n",
  
  "## 分析概述\n",
  "- **分析日期**: {Sys.Date()}\n",
  "- **暴露因素**: {exposure_name}\n",
  "- **结局因素**: {paste(outcome_names, collapse = ', ')}\n",
  "- **分析方法**: 单暴露-多结局MR分析（IVW, MR Egger, 加权中位数, 加权众数, 径向MR）\n\n",
  
  "## 数据概况\n",
  "### 工具变量选择结果\n",
  "| 暴露因素 | 工具变量数量 |\n",
  "|----------|--------------|\n",
  "| {exposure_name} | {nrow(selected_instruments)} |\n\n",
  
  "## 主要发现\n",
  "### IVW方法主要结果\n",
  "| 结局 | 效应量 | 标准误 | OR (95% CI) | P值 |\n",
  "|------|--------|--------|-------------|-----|\n",
  "{paste(sprintf('| %s | %.3f | %.3f | %.3f (%.3f-%.3f) | %.3f |', 
         ivw_results$outcome,
         ivw_results$b,
         ivw_results$se,
         exp(ivw_results$b),
         exp(ivw_results$ci_lower),
         exp(ivw_results$ci_upper),
         ivw_results$pval), collapse = '\\n')}\n\n",
  
  "### 关键结论\n",
  "{if (any(ivw_results$pval < 0.05)) {\n",
  "  significant_outcomes <- ivw_results %>% filter(pval < 0.05) %>% pull(outcome)\n",
  "  paste(sprintf('- %s与%s存在显著因果关联 (p < 0.05)', exposure_name, significant_outcomes), collapse = '\\n')\n",
  "} else {\n",
  "  sprintf('- 未发现%s与任何结局存在显著因果关联 (p ≥ 0.05)', exposure_name)\n",
  "}}\n\n",
  
  "### 结局间比较\n",
  "{if (nrow(ivw_results) > 1) {\n",
  "  max_effect_outcome <- ivw_results %>% arrange(desc(abs(b))) %>% slice(1) %>% pull(outcome)\n",
  "  max_effect_size <- ivw_results %>% arrange(desc(abs(b))) %>% slice(1) %>% pull(b)\n",
  "  sprintf('- %s对%s的效应量最大 (β = %.3f)', exposure_name, max_effect_outcome, max_effect_size)\n",
  "}}\n\n",
  
  "## 敏感性分析\n",
  "### 异质性检验\n",
  "{paste(lapply(seq_along(mr_results_list), function(i) {\n",
  "  outcome <- outcome_names[i]\n",
  "  hetero <- mr_results_list[[i]]$mr_result$heterogeneity %>% filter(method == 'Inverse variance weighted')\n",
  "  if (hetero$Q_pval < 0.05) {\n",
  "    sprintf('- %s: 存在显著异质性 (Q = %.3f, p = %.3f, I² = %.3f)', outcome, hetero$Q, hetero$Q_pval, hetero$I2)\n",
  "  } else {\n",
  "    sprintf('- %s: 未检测到显著异质性 (Q = %.3f, p = %.3f, I² = %.3f)', outcome, hetero$Q, hetero$Q_pval, hetero$I2)\n",
  "  }\n",
  "}), collapse = '\\n')}\n\n",
  
  "### 多效性检验\n",
  "{paste(lapply(seq_along(mr_results_list), function(i) {\n",
  "  outcome <- outcome_names[i]\n",
  "  pleio <- mr_results_list[[i]]$mr_result$pleiotropy\n",
  "  if (pleio$egger_intercept_pval < 0.05) {\n",
  "    sprintf('- %s: MR Egger截距检验提示可能存在水平多效性 (intercept = %.3f, p = %.3f)', outcome, pleio$egger_intercept, pleio$egger_intercept_pval)\n",
  "  } else {\n",
  "    sprintf('- %s: MR Egger截距检验未提示水平多效性 (intercept = %.3f, p = %.3f)', outcome, pleio$egger_intercept, pleio$egger_intercept_pval)\n",
  "  }\n",
  "}), collapse = '\\n')}\n\n",
  
  "### 径向MR分析\n",
  "{paste(lapply(seq_along(radial_results_summary$outcome), function(i) {\n",
  "  outcome <- radial_results_summary$outcome[i]\n",
  "  n_outliers <- radial_results_summary$n_outliers[i]\n",
  "  if (n_outliers > 0) {\n",
  "    sprintf('- %s: 检测到%d个异常值，可能影响结果稳健性', outcome, n_outliers)\n",
  "  } else {\n",
  "    sprintf('- %s: 未检测到异常值，结果较为稳健', outcome)\n",
  "  }\n",
  "}), collapse = '\\n')}\n\n",
  
  "## 局限性\n",
  "- 本分析基于示例数据，实际应用中应使用真实GWAS数据\n",
  "- 工具变量的有效性假设（相关性、独立性、排他性）需要进一步验证\n",
  "- 未考虑潜在的基因-环境交互作用\n",
  "- 不同结局之间可能存在共享遗传背景，需要进一步分析\n\n",
  
  "## 建议\n",
  "- 使用更大样本量的GWAS数据验证结果\n",
  "- 进行多变量MR分析，控制潜在的混杂因素\n",
  "- 结合功能基因组学数据，探究因果关系的分子机制\n",
  "- 在实验模型中验证关键发现\n\n",
  
  "## 分析流程\n",
  "```mermaid\n",
  "graph TD\n",
  "    A[数据准备] --> B[工具变量选择]\n",
  "    B --> C[单暴露-多结局MR分析]\n",
  "    C --> D[敏感性分析（异质性、多效性、径向MR）]\n",
  "    D --> E[结果可视化]\n",
  "    E --> F[结果汇总和报告]\n",
  "```\n\n",
  
  "## 输出文件\n",
  "- `multi_outcome_forest_plot.png`: 多结局MR分析森林图\n",
  "- `or_comparison_plot.png`: OR值比较图\n",
  "- `all_outcomes_mr_results.csv`: 所有结局的MR分析结果表\n",
  "- `radial_mr_results.csv`: 径向MR分析结果表\n",
  "- `instrumental_variables_info.csv`: 工具变量信息表\n",
  "- 各结局的MR分析综合图\n\n",
  
  "---\n",
  "*报告生成时间: {Sys.time()}*"
)

# 保存报告
writeLines(
  report_content,
  here(output_dir, "single_exposure_multi_outcome_mr_report.md")
)

# ----------------------
# 7. 完成分析
# ----------------------

cat("单暴露-多结局MR分析完成！\n")
cat("分析结果已保存到:", output_dir, "\n")
cat("主要输出文件包括:\n")
cat("- multi_outcome_forest_plot.png\n")
cat("- or_comparison_plot.png\n")
cat("- all_outcomes_mr_results.csv\n")
cat("- radial_mr_results.csv\n")
cat("- instrumental_variables_info.csv\n")
cat("- single_exposure_multi_outcome_mr_report.md\n")

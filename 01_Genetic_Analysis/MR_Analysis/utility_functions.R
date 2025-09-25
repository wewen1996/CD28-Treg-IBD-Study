#' CD28-Treg-IBD-Study: 孟德尔随机化分析工具函数
#' 
#' 这是一个包含孟德尔随机化(MR)分析中常用工具函数的集合，
#' 用于辅助MR分析的各个步骤，包括数据预处理、工具变量选择、
#' MR分析和结果可视化。

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

#' 数据预处理函数：格式化暴露和结局数据
#' 
#' @param exposure_data 暴露数据框，包含SNP、beta、se、effect_allele等列
#' @param outcome_data 结局数据框，包含SNP、beta、se、Effect_allele等列
#' @param exposure_name 暴露因素名称
#' @param outcome_name 结局因素名称
#' 
#' @return 格式化后的MR数据列表
format_mr_data <- function(exposure_data, outcome_data, exposure_name, outcome_name) {
  # 格式化暴露数据
  exposure_mr <- format_data(
    exposure_data,
    type = "exposure",
    snps = exposure_data$SNP,
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "Effect_allele",
    other_allele_col = "Other_allele",
    pval_col = "pval",
    phenotype_col = exposure_name
  )
  
  # 格式化结局数据
  outcome_mr <- format_data(
    outcome_data,
    type = "outcome",
    snps = outcome_data$SNP,
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "Effect_allele",
    other_allele_col = "Other_allele",
    pval_col = "pval",
    phenotype_col = outcome_name
  )
  
  # 合并暴露和结局数据
  mr_data <- harmonise_data(exposure_mr, outcome_mr)
  
  # 返回格式化后的数据
  return(list(
    exposure_data = exposure_mr,
    outcome_data = outcome_mr,
    mr_data = mr_data
  ))
}

#' 工具变量选择函数：基于LD和P值筛选独立SNP
#' 
#' @param gwas_data GWAS数据框，包含SNP、beta、se、pval等列
#' @param p_threshold P值阈值，默认5e-8
#' @param clump_r2 LD阈值，默认0.001
#' @param clump_kb 聚类窗口大小，默认10000kb
#' @param population 参考人群，默认"EUR"
#' 
#' @return 筛选后的工具变量数据框
select_instruments <- function(gwas_data, p_threshold = 5e-8, clump_r2 = 0.001, clump_kb = 10000, population = "EUR") {
  # 筛选达到显著水平的SNP
  gwas_significant <- gwas_data %>% filter(pval < p_threshold)
  
  # 基于LD聚类，选择独立SNP
  instruments_snps <- clump_data(
    gwas_significant,
    clump_r2 = clump_r2,
    clump_kb = clump_kb,
    pop = population
  )
  
  # 返回筛选后的工具变量
  return(instrument_snps)
}

#' MR分析函数：执行多种MR方法并汇总结果
#' 
#' @param mr_data 经过harmonise_data处理的MR数据
#' @param methods 要使用的MR方法列表，默认包含常用方法
#' 
#' @return MR分析结果列表
perform_mr_analysis <- function(mr_data, methods = c("mr_wald_ratio", "mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")) {
  # 执行MR分析
  mr_results <- mr(
    mr_data,
    method_list = methods
  )
  
  # 计算异质性指标
  heterogeneity <- mr_heterogeneity(mr_data)
  
  # 执行水平多效性检验
  pleiotropy <- mr_pleiotropy_test(mr_data)
  
  # 返回结果列表
  return(list(
    mr_results = mr_results,
    heterogeneity = heterogeneity,
    pleiotropy = pleiotropy
  ))
}

#' MR结果可视化函数：森林图
#' 
#' @param mr_results MR分析结果
#' @param title 图表标题
#' @param xlab x轴标签
#' 
#' @return 森林图ggplot对象
plot_mr_forest <- function(mr_results, title = "MR分析森林图", xlab = "效应量 (95% CI)") {
  # 提取结果数据
  forest_data <- mr_results$mr_results %>%
    filter(method %in% c("Wald ratio", "Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode")) %>%
    mutate(
      method = factor(method, levels = c("Wald ratio", "Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode")),
      ci_lower = b - 1.96 * se,
      ci_upper = b + 1.96 * se
    )
  
  # 创建森林图
  p <- ggplot(forest_data, aes(x = b, y = method)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
    geom_text(aes(label = sprintf("%.3f (%.3f, %.3f)", b, ci_lower, ci_upper)), 
              hjust = -0.1, size = 3.5) +
    labseller::scale_y_discrete(labels = c("Wald ratio" = "Wald比", 
                                          "Inverse variance weighted" = "逆方差加权",
                                          "MR Egger" = "MR Egger",
                                          "Weighted median" = "加权中位数",
                                          "Weighted mode" = "加权众数")) +
    labs(title = title, x = xlab, y = "MR方法") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major.x = element_line(linetype = "dashed", alpha = 0.5),
      panel.grid.minor.x = element_blank(),
      panel.grid.y = element_blank()
    ) +
    expand_limits(x = c(min(forest_data$ci_lower) * 1.2, max(forest_data$ci_upper) * 1.2))
  
  return(p)
}

#' MR结果可视化函数：散点图
#' 
#' @param mr_data 经过harmonise_data处理的MR数据
#' @param mr_results MR分析结果
#' @param title 图表标题
#' 
#' @return 散点图ggplot对象
plot_mr_scatter <- function(mr_data, mr_results, title = "MR分析散点图") {
  # 提取IVW结果用于绘制回归线
  ivw_result <- mr_results$mr_results %>%
    filter(method == "Inverse variance weighted")
  
  # 创建散点图
  p <- ggplot(mr_data, aes(x = beta.exposure, y = beta.outcome)) +
    geom_point(aes(size = 1/se.outcome), alpha = 0.7) +
    geom_errorbar(aes(ymin = beta.outcome - 1.96*se.outcome, ymax = beta.outcome + 1.96*se.outcome), alpha = 0.5) +
    geom_errorbarh(aes(xmin = beta.exposure - 1.96*se.exposure, xmax = beta.exposure + 1.96*se.exposure), alpha = 0.5) +
    geom_abline(intercept = ivw_result$intercept, slope = ivw_result$b, color = "red", linewidth = 1) +
    geom_text(aes(label = SNP), size = 3, hjust = -0.1, vjust = 0) +
    labs(
      title = title,
      x = "暴露效应量",
      y = "结局效应量",
      size = "1/SE结局"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
  
  return(p)
}

#' MR结果可视化函数：漏斗图
#' 
#' @param mr_data 经过harmonise_data处理的MR数据
#' @param title 图表标题
#' 
#' @return 漏斗图ggplot对象
plot_mr_funnel <- function(mr_data, title = "MR分析漏斗图") {
  # 计算标准化效应量
  mr_data$standardized_effect <- mr_data$beta.outcome / mr_data$beta.exposure
  
  # 创建漏斗图
  p <- ggplot(mr_data, aes(x = standardized_effect, y = 1/se.outcome)) +
    geom_point(alpha = 0.7, size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(
      title = title,
      x = "标准化效应量",
      y = "1/SE结局"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  return(p)
}

#' 径向MR分析函数
#' 
#' @param mr_data 经过harmonise_data处理的MR数据
#' 
#' @return 径向MR分析结果
perform_radial_mr <- function(mr_data) {
  # 转换数据格式
  radial_data <- radialize(mr_data)
  
  # 执行径向MR分析
  radial_result <- radial.mr(radial_data)
  
  # 识别异常值
  outliers <- radial_result$outliers
  
  # 返回结果
  return(list(
    radial_data = radial_data,
    radial_result = radial_result,
    outliers = outliers
  ))
}

#' 径向MR结果可视化函数
#' 
#' @param radial_result 径向MR分析结果
#' @param title 图表标题
#' 
#' @return 径向MR图ggplot对象
plot_radial_mr <- function(radial_result, title = "径向MR分析图") {
  # 创建径向MR图
  p <- plot(radial_result) +
    ggtitle(title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  return(p)
}

#' MR结果汇总函数：生成综合报告
#' 
#' @param mr_results_list MR分析结果列表
#' @param exposure_name 暴露因素名称
#' @param outcome_name 结局因素名称
#' @param output_file 输出文件路径
#' 
#' @return 汇总结果数据框
summarize_mr_results <- function(mr_results_list, exposure_name, outcome_name, output_file = NULL) {
  # 提取主要结果
  main_results <- mr_results_list$mr_results %>%
    filter(method %in% c("Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode")) %>%
    select(method, b, se, pval, nsnp) %>%
    mutate(
      OR = exp(b),
      OR_lower = exp(b - 1.96 * se),
      OR_upper = exp(b + 1.96 * se),
      exposure = exposure_name,
      outcome = outcome_name
    ) %>%
    select(exposure, outcome, method, nsnp, b, se, OR, OR_lower, OR_upper, pval)
  
  # 提取异质性结果
  heterogeneity <- mr_results_list$heterogeneity %>%
    select(method, Q, Q_df, Q_pval, I2)
  
  # 合并结果
  summary_results <- left_join(main_results, heterogeneity, by = "method")
  
  # 如果指定了输出文件，则保存结果
  if (!is.null(output_file)) {
    write.csv(summary_results, output_file, row.names = FALSE)
    cat("MR分析结果已保存到:", output_file, "\n")
  }
  
  # 返回汇总结果
  return(summary_results)
}

#' 多暴露MR分析函数
#' 
#' @param exposure_data_list 暴露数据列表，每个元素是一个暴露的GWAS数据框
#' @param outcome_data 结局数据框
#' @param exposure_names 暴露因素名称列表
#' @param outcome_name 结局因素名称
#' 
#' @return 多暴露MR分析结果
perform_multiple_exposure_mr <- function(exposure_data_list, outcome_data, exposure_names, outcome_name) {
  # 确保暴露数据列表和名称列表长度一致
  if (length(exposure_data_list) != length(exposure_names)) {
    stop("暴露数据列表和名称列表长度不一致")
  }
  
  # 格式化所有暴露数据
  formatted_exposures <- lapply(seq_along(exposure_data_list), function(i) {
    format_data(
      exposure_data_list[[i]],
      type = "exposure",
      snps = exposure_data_list[[i]]$SNP,
      beta_col = "beta",
      se_col = "se",
      effect_allele_col = "Effect_allele",
      other_allele_col = "Other_allele",
      pval_col = "pval",
      phenotype_col = exposure_names[i]
    )
  })
  
  # 格式化结局数据
  formatted_outcome <- format_data(
    outcome_data,
    type = "outcome",
    snps = outcome_data$SNP,
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "Effect_allele",
    other_allele_col = "Other_allele",
    pval_col = "pval",
    phenotype_col = outcome_name
  )
  
  # 合并所有数据
  combined_data <- harmonise_data(exposures = formatted_exposures, outcomes = formatted_outcome)
  
  # 执行多暴露MR分析
  mr_results <- mr(combined_data)
  
  # 返回结果
  return(list(
    combined_data = combined_data,
    mr_results = mr_results
  ))
}

#' 打印MR分析报告
#' 
#' @param mr_results_list MR分析结果列表
#' @param exposure_name 暴露因素名称
#' @param outcome_name 结局因素名称
print_mr_report <- function(mr_results_list, exposure_name, outcome_name) {
  cat("="*60, "\n")
  cat(sprintf("孟德尔随机化分析报告: %s -> %s\n", exposure_name, outcome_name))
  cat("="*60, "\n\n")
  
  # 打印主要结果
  cat("主要MR分析结果:\n")
  main_results <- mr_results_list$mr_results %>%
    filter(method %in% c("Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode")) %>%
    select(method, b, se, pval, nsnp) %>%
    mutate(
      OR = exp(b),
      OR_lower = exp(b - 1.96 * se),
      OR_upper = exp(b + 1.96 * se)
    )
  
  print(main_results)
  
  # 打印异质性结果
  cat("\n异质性检验结果:\n")
  heterogeneity <- mr_results_list$heterogeneity %>%
    select(method, Q, Q_df, Q_pval, I2)
  print(heterogeneity)
  
  # 打印多效性检验结果
  cat("\n多效性检验结果 (MR Egger intercept):\n")
  pleiotropy <- mr_results_list$pleiotropy %>%
    select(egger_intercept, egger_intercept_se, egger_intercept_pval)
  print(pleiotropy)
  
  # 风险评估
  cat("\n因果关系推断风险评估:\n")
  ivw_result <- main_results %>% filter(method == "Inverse variance weighted")
  if (ivw_result$pval < 0.05) {
    cat(sprintf("基于IVW方法，观察到%s对%s可能存在因果效应 (OR = %.3f, 95%% CI: %.3f-%.3f, p = %.3f)\n",
                exposure_name, outcome_name, ivw_result$OR, ivw_result$OR_lower, ivw_result$OR_upper, ivw_result$pval))
    
    # 评估异质性
    ivw_hetero <- heterogeneity %>% filter(method == "Inverse variance weighted")
    if (ivw_hetero$Q_pval < 0.05) {
      cat(sprintf("存在显著异质性 (Q = %.3f, p = %.3f, I2 = %.3f)，建议使用加权中位数或加权众数方法进行验证\n",
                  ivw_hetero$Q, ivw_hetero$Q_pval, ivw_hetero$I2))
    } else {
      cat(sprintf("未检测到显著异质性 (Q = %.3f, p = %.3f, I2 = %.3f)\n",
                  ivw_hetero$Q, ivw_hetero$Q_pval, ivw_hetero$I2))
    }
    
    # 评估多效性
    if (pleiotropy$egger_intercept_pval < 0.05) {
      cat(sprintf("MR Egger截距检验提示可能存在水平多效性 (intercept = %.3f, p = %.3f)\n",
                  pleiotropy$egger_intercept, pleiotropy$egger_intercept_pval))
    } else {
      cat(sprintf("MR Egger截距检验未提示水平多效性 (intercept = %.3f, p = %.3f)\n",
                  pleiotropy$egger_intercept, pleiotropy$egger_intercept_pval))
    }
  } else {
    cat(sprintf("基于IVW方法，未观察到%s对%s的显著因果效应 (OR = %.3f, 95%% CI: %.3f-%.3f, p = %.3f)\n",
                exposure_name, outcome_name, ivw_result$OR, ivw_result$OR_lower, ivw_result$OR_upper, ivw_result$pval))
  }
  
  cat("\n"*2)
}

# 函数加载完成提示
cat("MR分析工具函数加载完成！\n")
cat("可用函数包括：\n")
cat("- format_mr_data: 格式化MR分析数据\n")
cat("- select_instruments: 选择工具变量\n")
cat("- perform_mr_analysis: 执行MR分析\n")
cat("- plot_mr_forest: 绘制森林图\n")
cat("- plot_mr_scatter: 绘制散点图\n")
cat("- plot_mr_funnel: 绘制漏斗图\n")
cat("- perform_radial_mr: 执行径向MR分析\n")
cat("- plot_radial_mr: 绘制径向MR图\n")
cat("- summarize_mr_results: 汇总MR结果\n")
cat("- perform_multiple_exposure_mr: 执行多暴露MR分析\n")
cat("- print_mr_report: 打印MR分析报告\n")

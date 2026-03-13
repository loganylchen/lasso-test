#!/usr/bin/env Rscript
# LASSO Biomarker Selection - 快速演示
# 用法: Rscript demo.R

cat("=== LASSO Biomarker筛选快速演示 ===\n\n")

# 1. 加载函数和包
cat("1. 加载必要的包和函数...\n")
suppressPackageStartupMessages({
  library(glmnet)
  library(ggplot2)
  library(pROC)
  library(pheatmap)
  library(caret)
})

source("R/lasso_biomarker.R")
source("R/visualization.R")
source("R/evaluation.R")
source("data/generate_data.R")

# 2. 生成数据
cat("\n2. 生成示例数据...\n")
set.seed(42)
data <- generate_example_data(
  n_samples = 20,
  n_features = 500,
  n_informative = 20
)
cat("   - 样本数:", length(data$y), "(10正常 + 10疾病)\n")
cat("   - 特征数:", nrow(data$X), "\n")
cat("   - 真实biomarker数:", length(data$true_features), "\n")

# 3. LASSO交叉验证
cat("\n3. 运行LASSO交叉验证...\n")
cv_fit <- run_lasso_cv(data$X, data$y, nfolds = 5)
cat("   - λ.min:", round(cv_fit$lambda.min, 6), "\n")
cat("   - λ.1se:", round(cv_fit$lambda.1se, 6), "\n")

# 4. 筛选Biomarker
cat("\n4. 筛选Biomarker...\n")
biomarkers <- select_biomarkers(cv_fit, data$X, lambda = "lambda.1se")
cat("   - 筛选出", biomarkers$n_features, "个biomarker\n")

# 计算召回率
true_positives <- sum(biomarkers$selected_features %in% data$true_features)
recall <- true_positives / length(data$true_features)
cat("   - 真实biomarker召回率:", round(recall * 100, 1), "%\n")

# 5. 模型评估
cat("\n5. 评估模型性能...\n")
eval_result <- evaluate_model(cv_fit, data$X, data$y, lambda = "lambda.1se")
cat("   - AUC:", round(eval_result$auc, 3), "\n")
cat("   - 准确率:", round(eval_result$accuracy, 3), "\n")
cat("   - 灵敏度:", round(eval_result$sensitivity, 3), "\n")
cat("   - 特异度:", round(eval_result$specificity, 3), "\n")

# 6. 可视化
cat("\n6. 生成可视化图表...\n")

# 创建输出目录
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# ROC曲线
pdf("output/figures/roc_curve.pdf", width = 8, height = 6)
print(plot_roc_curve(cv_fit, data$X, data$y))
dev.off()
cat("   ✓ ROC曲线: output/figures/roc_curve.pdf\n")

# 系数路径图
pdf("output/figures/coefficient_path.pdf", width = 10, height = 6)
print(plot_coefficient_path(cv_fit))
dev.off()
cat("   ✓ 系数路径图: output/figures/coefficient_path.pdf\n")

# 交叉验证曲线
pdf("output/figures/cv_curve.pdf", width = 8, height = 6)
print(plot_cv_curve(cv_fit))
dev.off()
cat("   ✓ 交叉验证曲线: output/figures/cv_curve.pdf\n")

# Top biomarker箱线图
pdf("output/figures/biomarker_boxplot.pdf", width = 12, height = 6)
print(plot_biomarker_boxplot(data$X, data$y, biomarkers, top_n = 10))
dev.off()
cat("   ✓ Biomarker箱线图: output/figures/biomarker_boxplot.pdf\n")

# 热图（PNG格式，因为pheatmap不能直接用dev.off）
plot_biomarker_heatmap(
  data$X, 
  data$y, 
  biomarkers, 
  output_file = "output/figures/biomarker_heatmap.pdf"
)
cat("   ✓ Biomarker热图: output/figures/biomarker_heatmap.pdf\n")

# 7. 保存结果
cat("\n7. 保存结果...\n")
dir.create("output/results", showWarnings = FALSE, recursive = TRUE)

# 保存biomarker重要性
importance <- calculate_feature_importance(biomarkers)
write.csv(
  importance,
  "output/results/biomarker_importance.csv",
  row.names = FALSE
)
cat("   ✓ Biomarker重要性: output/results/biomarker_importance.csv\n")

# 保存评估指标
eval_summary <- data.frame(
  Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity", "N_Biomarkers", "Recall"),
  Value = c(
    eval_result$auc,
    eval_result$accuracy,
    eval_result$sensitivity,
    eval_result$specificity,
    biomarkers$n_features,
    recall
  )
)
write.csv(
  eval_summary,
  "output/results/evaluation_summary.csv",
  row.names = FALSE
)
cat("   ✓ 评估指标: output/results/evaluation_summary.csv\n")

# 保存模型
saveRDS(cv_fit, "output/results/lasso_model.rds")
cat("   ✓ LASSO模型: output/results/lasso_model.rds\n")

# 8. 总结
cat("\n=== 演示完成 ===\n")
cat("\n结果文件位置:\n")
cat("  - 图表: output/figures/\n")
cat("  - 结果: output/results/\n")
cat("\n下一步:\n")
cat("  - 查看完整教程: notebooks/tutorial.Rmd\n")
cat("  - 运行测试: Rscript -e \"testthat::test_file('tests/testthat/test_lasso.R')\"\n")

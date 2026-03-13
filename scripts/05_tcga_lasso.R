#!/usr/bin/env Rscript
# TCGA CRC数据 - LASSO Biomarker筛选

cat("=== TCGA CRC数据 - LASSO Biomarker筛选 ===\n\n")

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(dplyr)
})

# 加载LASSO函数
source("R/lasso_biomarker.R")
source("R/visualization.R")
source("R/evaluation.R")

save_dir <- "~/.openclaw/workspace/lasso-biomarker-tutorial/data/tcga_crc"

cat("1. 加载TCGA DESeq2结果...\n")
dds <- readRDS(file.path(save_dir, "tcga_coad_dds_analyzed.rds"))
deseq_res <- readRDS(file.path(save_dir, "tcga_coad_deseq_results.rds"))

cat("   ✓ 数据加载完成\n")
cat("   - 样本数:", ncol(dds), "\n")
cat("   - 基因数:", nrow(dds), "\n\n")

cat("2. 准备LASSO输入数据...\n")

sig_genes_df <- read.csv(file.path(save_dir, "tcga_coad_sig_genes.csv"), row.names = 1)
cat("   总差异基因数:", nrow(sig_genes_df), "\n")

# 选择top 500差异基因（按padj排序）
n_features <- min(500, nrow(sig_genes_df))
sig_gene_names <- rownames(sig_genes_df)[1:n_features]

cat("   使用", length(sig_gene_names), "个差异基因作为候选特征\n\n")

cat("3. 标准化数据（VST变换）...\n")
vsd <- vst(dds, blind = FALSE)
expr_matrix <- assay(vsd)

X <- expr_matrix[sig_gene_names, ]
cat("   ✓ 特征矩阵准备完成\n")
cat("   - 特征数:", nrow(X), "\n")
cat("   - 样本数:", ncol(X), "\n\n")

y <- as.integer(dds$condition == "Tumor")
cat("   标签分布:\n")
print(table(y))
cat("\n")

cat("4. 运行LASSO交叉验证...\n")

n_samples <- length(y)
n_per_class <- min(sum(y == 0), sum(y == 1))
nfolds <- min(10, floor(n_per_class / 3))  # 真实数据用10折

cat("   使用", nfolds, "折交叉验证\n")

cv_fit <- run_lasso_cv(
  X = X,
  y = y,
  nfolds = nfolds,
  alpha = 1
)

cat("   ✓ LASSO交叉验证完成\n")
cat("   - λ.min:", round(cv_fit$lambda.min, 6), "\n")
cat("   - λ.1se:", round(cv_fit$lambda.1se, 6), "\n\n")

cat("5. 筛选Biomarker（尝试两种lambda）...\n")

biomarkers_1se <- select_biomarkers(cv_fit, X, lambda = "lambda.1se")
cat("   lambda.1se: ", biomarkers_1se$n_features, "个\n")

biomarkers_min <- select_biomarkers(cv_fit, X, lambda = "lambda.min")
cat("   lambda.min: ", biomarkers_min$n_features, "个\n\n")

# 选择有足够biomarker的版本
if (biomarkers_1se$n_features >= 5) {
  biomarkers <- biomarkers_1se
  lambda_choice <- "lambda.1se"
  cat("   使用lambda.1se（更保守）\n")
} else if (biomarkers_min$n_features >= 5) {
  biomarkers <- biomarkers_min
  lambda_choice <- "lambda.min"
  cat("   使用lambda.min（更宽松）\n")
} else {
  # 如果还是太少，选择top 15差异基因
  cat("   ⚠️ LASSO筛选biomarker过少，手动选择top 15差异基因\n")
  top_genes <- rownames(sig_genes_df)[1:15]
  biomarkers <- list(
    selected_features = top_genes,
    coefficients = rep(1, 15),
    n_features = 15,
    lambda = cv_fit$lambda.min
  )
  lambda_choice <- "lambda.min"
}

cat("\n   ✓ 最终筛选出", biomarkers$n_features, "个biomarker\n")
cat("\n   Biomarker列表:\n")
print(biomarkers$selected_features)
cat("\n")

cat("6. 模型评估...\n")

eval_result <- evaluate_model(
  cv_fit = cv_fit,
  X = X,
  y = y,
  lambda = lambda_choice
)

cat("   ✓ 模型性能:\n")
cat("   - AUC:", round(eval_result$auc, 3), "\n")
cat("   - 准确率:", round(eval_result$accuracy, 3), "\n")
cat("   - 灵敏度:", round(eval_result$sensitivity, 3), "\n")
cat("   - 特异度:", round(eval_result$specificity, 3), "\n\n")

cat("   混淆矩阵:\n")
print(eval_result$confusion_matrix)
cat("\n")

cat("7. 生成可视化...\n")

# ROC曲线
pdf(file.path(save_dir, "tcga_lasso_roc_curve.pdf"), width = 8, height = 6)
print(plot_roc_curve(cv_fit, X, y, lambda = lambda_choice))
dev.off()
cat("   ✓ ROC曲线\n")

# 系数路径
pdf(file.path(save_dir, "tcga_lasso_coefficient_path.pdf"), width = 12, height = 6)
print(plot_coefficient_path(cv_fit))
dev.off()
cat("   ✓ 系数路径\n")

# 交叉验证曲线
pdf(file.path(save_dir, "tcga_lasso_cv_curve.pdf"), width = 8, height = 6)
print(plot_cv_curve(cv_fit))
dev.off()
cat("   ✓ 交叉验证曲线\n")

# Biomarker热图
if (biomarkers$n_features >= 2) {
  tryCatch({
    plot_biomarker_heatmap(
      X = X,
      y = y,
      biomarkers = biomarkers,
      output_file = file.path(save_dir, "tcga_lasso_biomarker_heatmap.pdf")
    )
    cat("   ✓ Biomarker热图\n")
  }, error = function(e) {
    cat("   ⚠️ 热图生成失败\n")
  })
  
  # Biomarker箱线图
  if (biomarkers$n_features >= 3) {
    pdf(file.path(save_dir, "tcga_lasso_biomarker_boxplot.pdf"), width = 14, height = 6)
    print(plot_biomarker_boxplot(X, y, biomarkers, top_n = min(15, biomarkers$n_features)))
    dev.off()
    cat("   ✓ Biomarker箱线图\n")
  }
}

cat("\n8. 保存结果...\n")

# 保存biomarker重要性
importance <- calculate_feature_importance(biomarkers)
write.csv(
  importance,
  file.path(save_dir, "tcga_lasso_biomarker_importance.csv"),
  row.names = FALSE
)

# 保存评估指标
eval_summary <- data.frame(
  Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity", 
             "N_Biomarkers", "N_Samples", "N_Features", "N_Normal", "N_Tumor"),
  Value = c(
    eval_result$auc,
    eval_result$accuracy,
    eval_result$sensitivity,
    eval_result$specificity,
    biomarkers$n_features,
    length(y),
    nrow(X),
    sum(y == 0),
    sum(y == 1)
  )
)
write.csv(
  eval_summary,
  file.path(save_dir, "tcga_lasso_evaluation_summary.csv"),
  row.names = FALSE
)

# 保存LASSO模型
saveRDS(cv_fit, file.path(save_dir, "tcga_lasso_cv_fit.rds"))
saveRDS(biomarkers, file.path(save_dir, "tcga_lasso_biomarkers.rds"))

cat("   ✓ 所有结果已保存\n\n")

cat("=== TCGA LASSO分析完成 ===\n")
cat("\n📊 最终总结:\n")
cat("  - 数据来源: TCGA-COAD (真实结肠癌数据)\n")
cat("  - 样本数:", length(y), "(", sum(y == 0), "Normal +", sum(y == 1), "Tumor )\n")
cat("  - 输入特征数:", nrow(X), "\n")
cat("  - 筛选出biomarker数:", biomarkers$n_features, "\n")
cat("  - 模型AUC:", round(eval_result$auc, 3), "\n")
cat("  - 模型准确率:", round(eval_result$accuracy, 3), "\n")
cat("\n所有结果已保存到:", save_dir, "\n")

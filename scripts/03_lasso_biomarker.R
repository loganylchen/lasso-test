#!/usr/bin/env Rscript
# CRC真实数据分析 - Step 3: LASSO Biomarker筛选
# 
# 输入：DESeq2结果（crc_dds.rds, crc_sig_genes.csv）
# 输出：LASSO筛选的biomarker

cat("=== CRC数据 - LASSO Biomarker筛选 ===\n\n")

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(dplyr)
})

# 加载LASSO函数
source("R/lasso_biomarker.R")
source("R/visualization.R")
source("R/evaluation.R")

save_dir <- "~/.openclaw/workspace/lasso-biomarker-tutorial/data/crc_analysis"

cat("1. 加载DESeq2结果...\n")
dds <- readRDS(file.path(save_dir, "crc_dds.rds"))
deseq_res <- readRDS(file.path(save_dir, "crc_deseq_results.rds"))

cat("   ✓ 数据加载完成\n")
cat("   - 样本数:", ncol(dds), "\n")
cat("   - 基因数:", nrow(dds), "\n\n")

cat("2. 准备LASSO输入数据...\n")

# 读取显著差异基因
sig_genes_file <- file.path(save_dir, "crc_sig_genes.csv")
if (file.exists(sig_genes_file)) {
  sig_genes_df <- read.csv(sig_genes_file, row.names = 1)
  sig_gene_names <- rownames(sig_genes_df)
  
  # 如果差异基因太多，选择top 500（按padj排序）
  if (length(sig_gene_names) > 500) {
    sig_gene_names <- rownames(sig_genes_df)[1:500]
    cat("   差异基因过多，选择top 500\n")
  }
} else {
  # 备选：选择padj < 0.01的基因
  sig_genes <- subset(deseq_res, padj < 0.01 & abs(log2FoldChange) > 1)
  sig_gene_names <- rownames(sig_genes)
  
  if (length(sig_gene_names) > 500) {
    sig_gene_names <- rownames(sig_genes)[1:500]
  }
}

cat("   使用", length(sig_gene_names), "个差异基因作为候选特征\n\n")

# 提取标准化的count数据（VST变换）
cat("3. 标准化数据（VST变换）...\n")
vsd <- vst(dds, blind = FALSE)
expr_matrix <- assay(vsd)

# 提取差异基因的表达矩阵
X <- expr_matrix[sig_gene_names, ]
cat("   ✓ 特征矩阵准备完成\n")
cat("   - 特征数:", nrow(X), "\n")
cat("   - 样本数:", ncol(X), "\n\n")

# 准备标签（0=Normal, 1=Tumor）
y <- as.integer(dds$condition == "Tumor")
cat("   标签分布:\n")
print(table(y))
cat("\n")

cat("4. 运行LASSO交叉验证...\n")

# 计算合适的nfolds（确保每折至少有5个样本）
n_samples <- length(y)
n_per_class <- min(sum(y == 0), sum(y == 1))
nfolds <- min(5, floor(n_per_class / 2))

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

cat("5. 筛选Biomarker...\n")

biomarkers <- select_biomarkers(
  cv_fit = cv_fit,
  X = X,
  lambda = "lambda.1se"
)

cat("   ✓ 筛选出", biomarkers$n_features, "个biomarker\n")
cat("\n   Biomarker列表:\n")
print(biomarkers$selected_features)
cat("\n")

cat("6. 模型评估...\n")

eval_result <- evaluate_model(
  cv_fit = cv_fit,
  X = X,
  y = y,
  lambda = "lambda.1se"
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
pdf(file.path(save_dir, "lasso_roc_curve.pdf"), width = 8, height = 6)
print(plot_roc_curve(cv_fit, X, y, lambda = "lambda.1se"))
dev.off()
cat("   ✓ ROC曲线: lasso_roc_curve.pdf\n")

# 系数路径
pdf(file.path(save_dir, "lasso_coefficient_path.pdf"), width = 10, height = 6)
print(plot_coefficient_path(cv_fit))
dev.off()
cat("   ✓ 系数路径: lasso_coefficient_path.pdf\n")

# 交叉验证曲线
pdf(file.path(save_dir, "lasso_cv_curve.pdf"), width = 8, height = 6)
print(plot_cv_curve(cv_fit))
dev.off()
cat("   ✓ 交叉验证曲线: lasso_cv_curve.pdf\n")

# Biomarker热图
if (biomarkers$n_features > 0) {
  plot_biomarker_heatmap(
    X = X,
    y = y,
    biomarkers = biomarkers,
    output_file = file.path(save_dir, "lasso_biomarker_heatmap.pdf")
  )
  cat("   ✓ Biomarker热图: lasso_biomarker_heatmap.pdf\n")
  
  # Biomarker箱线图
  if (biomarkers$n_features >= 5) {
    pdf(file.path(save_dir, "lasso_biomarker_boxplot.pdf"), width = 12, height = 6)
    print(plot_biomarker_boxplot(X, y, biomarkers, top_n = min(10, biomarkers$n_features)))
    dev.off()
    cat("   ✓ Biomarker箱线图: lasso_biomarker_boxplot.pdf\n")
  }
}

cat("\n8. 保存结果...\n")

# 保存biomarker重要性
importance <- calculate_feature_importance(biomarkers)
write.csv(
  importance,
  file.path(save_dir, "lasso_biomarker_importance.csv"),
  row.names = FALSE
)
cat("   ✓ Biomarker重要性: lasso_biomarker_importance.csv\n")

# 保存评估指标
eval_summary <- data.frame(
  Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity", "N_Biomarkers", "N_Samples", "N_Features"),
  Value = c(
    eval_result$auc,
    eval_result$accuracy,
    eval_result$sensitivity,
    eval_result$specificity,
    biomarkers$n_features,
    length(y),
    nrow(X)
  )
)
write.csv(
  eval_summary,
  file.path(save_dir, "lasso_evaluation_summary.csv"),
  row.names = FALSE
)
cat("   ✓ 评估指标: lasso_evaluation_summary.csv\n")

# 保存LASSO模型
saveRDS(cv_fit, file.path(save_dir, "lasso_cv_fit.rds"))
saveRDS(biomarkers, file.path(save_dir, "lasso_biomarkers.rds"))
cat("   ✓ LASSO模型和biomarker已保存\n\n")

cat("=== LASSO分析完成 ===\n")
cat("\n总结:\n")
cat("  - 输入特征数:", nrow(X), "\n")
cat("  - 样本数:", length(y), "(", sum(y == 0), "Normal +", sum(y == 1), "Tumor )\n")
cat("  - 筛选出biomarker数:", biomarkers$n_features, "\n")
cat("  - 模型AUC:", round(eval_result$auc, 3), "\n")
cat("  - 模型准确率:", round(eval_result$accuracy, 3), "\n")
cat("\n所有结果已保存到:", save_dir, "\n")

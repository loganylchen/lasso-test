#!/usr/bin/env Rscript
# 从Python下载的TCGA数据创建DESeq2对象

cat("=== 从Python下载的TCGA数据创建DESeq2对象 ===\n\n")

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
})

python_data_dir <- "~/.openclaw/workspace/lasso-biomarker-tutorial/data/tcga_coad_python"
save_dir <- python_data_dir

cat("1. 加载Python下载的数据...\n")

# 读取count矩阵
count_matrix <- read.csv(file.path(python_data_dir, "count_matrix.csv"), 
                         row.names = 1, check.names = FALSE)

# 读取样本信息
sample_info <- read.csv(file.path(python_data_dir, "sample_info.csv"))

cat("   ✓ 数据加载完成\n")
cat("   - 基因数:", nrow(count_matrix), "\n")
cat("   - 样本数:", ncol(count_matrix), "\n\n")

cat("   样本类型:\n")
print(table(sample_info$condition))
cat("\n")

cat("2. 创建DESeq2对象...\n")

# 确保count矩阵是整数
count_matrix <- round(count_matrix)

# 创建colData
condition <- factor(sample_info$condition, levels = c("Normal", "Tumor"))
coldata_df <- data.frame(
  condition = condition,
  row.names = sample_info$sample_id
)

# 确保列名匹配
colnames(count_matrix) <- sample_info$sample_id

# 创建DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata_df,
  design = ~ condition
)

cat("   ✓ DESeqDataSet创建完成\n")
cat("   - Normal:", sum(condition == "Normal"), "\n")
cat("   - Tumor:", sum(condition == "Tumor"), "\n\n")

cat("3. 预过滤低表达基因...\n")
keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep, ]
cat("   保留", nrow(dds), "个基因\n\n")

cat("4. 保存DESeq2对象...\n")

saveRDS(dds, file.path(save_dir, "tcga_python_dds.rds"))

cat("   ✓ 已保存到:", file.path(save_dir, "tcga_python_dds.rds"), "\n\n")

cat("=== DESeq2对象创建完成 ===\n")
cat("\n下一步：运行DESeq2分析\n")
cat("Rscript scripts/08_tcga_python_deseq2.R\n")

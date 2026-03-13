#!/usr/bin/env Rscript
# 使用curatedTCGAData下载TCGA CRC数据（预处理好的）
# 这个包不需要httr，直接从ExperimentHub下载

cat("=== 使用curatedTCGAData获取TCGA CRC数据 ===\n\n")

cat("1. 安装curatedTCGAData...\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

required_packages <- c("curatedTCGAData", "MultiAssayExperiment", "DESeq2", "SummarizedExperiment")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("   安装", pkg, "...\n")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

cat("   ✓ 包准备完成\n\n")

suppressPackageStartupMessages({
  library(curatedTCGAData)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
  library(DESeq2)
})

save_dir <- "~/.openclaw/workspace/lasso-biomarker-tutorial/data/tcga_crc"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

cat("2. 下载TCGA-COAD RNA-seq数据...\n")

# 下载COAD（结肠癌）的RNASeq2GeneNorm数据
coad_data <- curatedTCGAData(
  diseaseCode = "COAD",
  assays = "RNASeq2GeneNorm",
  version = "2.0.1",
  dry.run = FALSE
)

cat("   ✓ 数据下载完成\n\n")

cat("3. 提取RNA-seq数据...\n")

# 提取RNA-seq assay
rnaseq <- experiments(coad_data)[[1]]

cat("   原始数据:\n")
cat("   - 样本数:", ncol(rnaseq), "\n")
cat("   - 基因数:", nrow(rnaseq), "\n\n")

# 提取样本信息
sample_info <- colData(rnaseq)

# 识别正常和肿瘤样本
# TCGA barcode: TCGA-XX-XXXX-01A (肿瘤) vs TCGA-XX-XXXX-11A (正常)
barcode_type <- substr(colnames(rnaseq), 14, 15)
sample_type <- ifelse(as.integer(barcode_type) >= 10, "Normal", "Tumor")

cat("   样本类型分布:\n")
print(table(sample_type))
cat("\n")

# 选择样本（每组最多30个）
cat("4. 选择样本...\n")

normal_idx <- which(sample_type == "Normal")
tumor_idx <- which(sample_type == "Tumor")

set.seed(42)
selected_normal <- sample(normal_idx, min(30, length(normal_idx)))
selected_tumor <- sample(tumor_idx, min(30, length(tumor_idx)))
selected_idx <- c(selected_normal, selected_tumor)

rnaseq_subset <- rnaseq[, selected_idx]
sample_type_subset <- sample_type[selected_idx]

cat("   选择", length(selected_normal), "个正常样本\n")
cat("   选择", length(selected_tumor), "个肿瘤样本\n\n")

cat("5. 准备count矩阵...\n")

# curatedTCGAData提供的是normalized数据，我们需要转换回count
# 对于DESeq2，我们需要raw counts
# 这里使用RSEM expected_count（如果有）或者四舍五入normalized值

expr_data <- assay(rnaseq_subset)

# 检查数据类型
cat("   数据范围:", range(expr_data, na.rm = TRUE), "\n")
cat("   数据均值:", mean(expr_data, na.rm = TRUE), "\n\n")

# 如果数据是log-transformed，先转回来
if (max(expr_data, na.rm = TRUE) < 50) {
  cat("   检测到log-transformed数据，转换回线性空间\n")
  expr_data <- 2^expr_data - 1  # 假设是log2(x+1)
}

# 四舍五入为整数（模拟count）
counts_matrix <- round(expr_data)
counts_matrix[counts_matrix < 0] <- 0

cat("   转换后count范围:", range(counts_matrix, na.rm = TRUE), "\n\n")

cat("6. 创建DESeq2对象...\n")

# 创建colData
condition <- factor(
  sample_type_subset,
  levels = c("Normal", "Tumor")
)

coldata_df <- data.frame(
  condition = condition,
  row.names = colnames(counts_matrix)
)

# 创建DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = coldata_df,
  design = ~ condition
)

cat("   ✓ DESeqDataSet创建完成\n")
cat("   - Normal:", sum(condition == "Normal"), "\n")
cat("   - Tumor:", sum(condition == "Tumor"), "\n\n")

cat("7. 保存数据...\n")

saveRDS(dds, file.path(save_dir, "tcga_coad_dds.rds"))
saveRDS(rnaseq_subset, file.path(save_dir, "tcga_coad_raw.rds"))

cat("   ✓ 数据已保存\n")
cat("   - tcga_coad_dds.rds (DESeq2对象)\n")
cat("   - tcga_coad_raw.rds (原始数据)\n\n")

cat("=== 数据准备完成 ===\n")
cat("\n下一步：运行DESeq2分析\n")
cat("Rscript scripts/04_tcga_deseq2.R\n")

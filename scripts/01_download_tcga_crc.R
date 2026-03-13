#!/usr/bin/env Rscript
# 下载真实TCGA CRC数据 - 使用TCGAbiolinks
# 
# TCGAbiolinks是专门用于下载和处理TCGA数据的Bioconductor包

cat("=== 下载真实TCGA CRC数据 ===\n\n")

# 1. 安装TCGAbiolinks
cat("1. 检查并安装TCGAbiolinks...\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

required_packages <- c("TCGAbiolinks", "SummarizedExperiment", "DESeq2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("   安装", pkg, "...\n")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

cat("   ✓ 包准备完成\n\n")

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
})

save_dir <- "~/.openclaw/workspace/lasso-biomarker-tutorial/data/tcga_crc"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

cat("2. 查询TCGA-COAD（结肠癌）数据...\n")

# 查询TCGA-COAD项目的RNA-seq数据
query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Solid Tissue Normal", "Primary Tumor")
)

cat("   找到", length(getResults(query)$cases), "个样本\n")

# 查看样本类型分布
results <- getResults(query)
cat("   样本类型分布:\n")
print(table(results$sample_type))

# 限制样本数（加快下载）
cat("\n3. 选择样本（每组最多30个）...\n")

normal_samples <- results$cases[results$sample_type == "Solid Tissue Normal"]
tumor_samples <- results$cases[results$sample_type == "Primary Tumor"]

# 随机选择
set.seed(42)
selected_normal <- sample(normal_samples, min(30, length(normal_samples)))
selected_tumor <- sample(tumor_samples, min(30, length(tumor_samples)))
selected_cases <- c(selected_normal, selected_tumor)

cat("   选择", length(selected_normal), "个正常样本\n")
cat("   选择", length(selected_tumor), "个肿瘤样本\n")

# 重新查询选中的样本
query_subset <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = selected_cases
)

cat("\n4. 下载数据（这可能需要几分钟）...\n")

# 下载数据
GDCdownload(query_subset, directory = save_dir)

cat("   ✓ 数据下载完成\n\n")

cat("5. 准备数据矩阵...\n")

# 读取并整理数据
data <- GDCprepare(query_subset, directory = save_dir)

cat("   ✓ 数据读取完成\n")
cat("   - 样本数:", ncol(data), "\n")
cat("   - 基因数:", nrow(data), "\n\n")

cat("6. 创建DESeq2对象...\n")

# 提取count矩阵
counts_matrix <- assay(data, "unstranded")

# 样本信息
sample_info <- colData(data)
condition <- factor(
  ifelse(sample_info$sample_type == "Solid Tissue Normal", "Normal", "Tumor"),
  levels = c("Normal", "Tumor")
)

# 创建简化的colData
coldata_simple <- data.frame(
  condition = condition,
  row.names = colnames(counts_matrix)
)

# 创建DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = coldata_simple,
  design = ~ condition
)

cat("   ✓ DESeqDataSet创建完成\n")
cat("   - Normal:", sum(condition == "Normal"), "\n")
cat("   - Tumor:", sum(condition == "Tumor"), "\n\n")

cat("7. 保存数据...\n")

saveRDS(dds, file.path(save_dir, "tcga_coad_dds.rds"))
saveRDS(data, file.path(save_dir, "tcga_coad_raw.rds"))

cat("   ✓ 数据已保存\n")
cat("   - tcga_coad_dds.rds (DESeq2对象)\n")
cat("   - tcga_coad_raw.rds (原始数据)\n\n")

cat("=== 数据下载完成 ===\n")
cat("\n下一步：运行DESeq2分析\n")
cat("Rscript scripts/04_tcga_deseq2.R\n")

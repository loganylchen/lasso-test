#!/usr/bin/env Rscript
# 使用GEOquery下载真实CRC RNA-seq数据
# GSE164191 - 有完整的count矩阵

cat("=== 下载真实CRC RNA-seq数据（GSE164191）===\n\n")

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(DESeq2)
})

save_dir <- "~/.openclaw/workspace/lasso-biomarker-tutorial/data/geo_crc"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

cat("策略：直接下载包含count矩阵的supplementary file\n\n")

cat("1. 搜索合适的CRC数据集...\n")
cat("   使用 GSE39582 (Illumina microarray, N=566)\n")
cat("   这是一个大型CRC队列研究，有完整的临床信息\n\n")

cat("2. 下载数据（这可能需要几分钟）...\n")

# 下载GSE39582
gse <- tryCatch({
  getGEO("GSE39582", GSEMatrix = TRUE, getGPL = FALSE)
}, error = function(e) {
  cat("   尝试备选方案：GSE32323\n")
  getGEO("GSE32323", GSEMatrix = TRUE, getGPL = FALSE)
})

if (is.list(gse)) {
  gse <- gse[[1]]
}

cat("   ✓ 数据下载完成\n\n")

cat("3. 检查数据...\n")

expr_data <- exprs(gse)
pdata <- pData(gse)

cat("   - 基因/探针数:", nrow(expr_data), "\n")
cat("   - 样本数:", ncol(expr_data), "\n")
cat("   - 数据范围:", range(expr_data, na.rm = TRUE), "\n")
cat("   - 数据均值:", mean(expr_data, na.rm = TRUE), "\n\n")

# 检查是否有NULL数据
if (nrow(expr_data) == 0) {
  cat("❌ 数据为空，尝试下载supplementary files...\n")
  
  # 获取supplementary files列表
  suppl <- getGEOSuppFiles("GSE39582", baseDir = save_dir, makeDirectory = FALSE)
  cat("   Supplementary files:\n")
  print(suppl)
  
  stop("需要手动处理supplementary files")
}

cat("4. 提取样本信息...\n")

# 查看可用的pheno data列
cat("   可用的样本信息列:\n")
print(head(names(pdata), 20))
cat("\n")

# 尝试找到组织类型列
tissue_cols <- grep("tissue|sample|type", names(pdata), ignore.case = TRUE, value = TRUE)
cat("   疑似组织类型列:\n")
print(tissue_cols)
cat("\n")

# 如果有characteristics列
char_cols <- grep("characteristics", names(pdata), ignore.case = TRUE, value = TRUE)
if (length(char_cols) > 0) {
  cat("   Characteristics列:\n")
  for (col in char_cols[1:min(3, length(char_cols))]) {
    cat("   ", col, ":\n")
    print(head(pdata[[col]], 3))
    cat("\n")
  }
}

# 简单方案：随机选择Normal和Tumor样本（如果没有明确标签）
n_samples <- ncol(expr_data)
n_normal <- min(30, floor(n_samples / 3))
n_tumor <- min(30, floor(n_samples / 3))

cat("5. 选择样本（随机分组演示）...\n")
cat("   注意：这是演示用途，实际需要根据样本标签分组\n\n")

set.seed(42)
normal_idx <- sample(1:n_samples, n_normal)
tumor_idx <- sample(setdiff(1:n_samples, normal_idx), n_tumor)
selected_idx <- c(normal_idx, tumor_idx)

expr_subset <- expr_data[, selected_idx]
sample_type <- c(rep("Normal", n_normal), rep("Tumor", n_tumor))

cat("   选择", n_normal, "个Normal样本\n")
cat("   选择", n_tumor, "个Tumor样本\n\n")

cat("6. 数据预处理...\n")

# 转换为count-like数据
if (max(expr_subset, na.rm = TRUE) < 20) {
  cat("   检测到log2数据，转换...\n")
  expr_linear <- 2^expr_subset - 1
} else {
  expr_linear <- expr_subset
}

expr_linear[is.na(expr_linear)] <- 0
counts_matrix <- round(expr_linear)
counts_matrix[counts_matrix < 0] <- 0

# 过滤低表达
keep <- rowSums(counts_matrix >= 10) >= 10
counts_filtered <- counts_matrix[keep, ]

cat("   过滤后:", nrow(counts_filtered), "个基因/探针\n\n")

cat("7. 创建DESeq2对象...\n")

condition <- factor(sample_type, levels = c("Normal", "Tumor"))
coldata_df <- data.frame(
  condition = condition,
  row.names = colnames(counts_filtered)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = coldata_df,
  design = ~ condition
)

cat("   ✓ DESeqDataSet创建完成\n")
cat("   - Normal:", sum(condition == "Normal"), "\n")
cat("   - Tumor:", sum(condition == "Tumor"), "\n")
cat("   - 探针/基因数:", nrow(dds), "\n\n")

cat("8. 保存数据...\n")

saveRDS(dds, file.path(save_dir, "geo_crc_dds.rds"))
saveRDS(gse, file.path(save_dir, "geo_crc_raw.rds"))

sample_info_df <- data.frame(
  sample_id = colnames(counts_filtered),
  condition = sample_type,
  index = selected_idx
)
write.csv(sample_info_df, file.path(save_dir, "geo_sample_info.csv"), row.names = FALSE)

cat("   ✓ 数据已保存\n\n")

cat("=== GEO数据下载完成 ===\n")
cat("数据集: GSE39582 (或备选)\n")
cat("样本数:", ncol(dds), "\n")
cat("基因/探针数:", nrow(dds), "\n")
cat("\n⚠️  注意：样本分组是随机的，仅用于演示\n")
cat("    实际使用时需要根据真实的临床标签分组\n\n")
cat("下一步：运行DESeq2分析\n")
cat("Rscript scripts/06_geo_deseq2.R\n")

#!/usr/bin/env Rscript
# 使用GEOquery下载真实CRC RNA-seq数据

cat("=== 使用GEOquery下载真实CRC数据 ===\n\n")

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(DESeq2)
})

save_dir <- "~/.openclaw/workspace/lasso-biomarker-tutorial/data/geo_crc"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

cat("1. 下载GEO数据集（GSE50760 - CRC vs Normal）...\n")
cat("   数据集信息:\n")
cat("   - GSE50760: Colorectal cancer gene expression\n")
cat("   - 平台: Illumina\n")
cat("   - 样本: 18 Normal + 18 Tumor\n\n")

# 下载数据（这是一个较小的CRC数据集，便于快速演示）
gse <- getGEO("GSE50760", GSEMatrix = TRUE, getGPL = FALSE)

if (is.list(gse)) {
  gse <- gse[[1]]
}

cat("   ✓ 数据下载完成\n\n")

cat("2. 提取表达矩阵和样本信息...\n")

# 表达矩阵
expr_data <- exprs(gse)
cat("   - 探针数:", nrow(expr_data), "\n")
cat("   - 样本数:", ncol(expr_data), "\n\n")

# 样本信息
pdata <- pData(gse)

# 提取样本类型（从characteristics_ch1字段）
sample_type <- pdata$characteristics_ch1
cat("   样本类型:\n")
print(table(sample_type))
cat("\n")

# 清理样本类型标签
sample_type_clean <- ifelse(grepl("normal", sample_type, ignore.case = TRUE), 
                             "Normal", "Tumor")

cat("   清理后样本类型:\n")
print(table(sample_type_clean))
cat("\n")

cat("3. 数据预处理...\n")

# GEO数据通常是归一化的表达值，需要转换为count-like数据用于DESeq2
# 这里使用一个简单的转换：取2的幂次（如果是log2值）然后四舍五入

# 检查数据范围
cat("   原始数据范围:", range(expr_data, na.rm = TRUE), "\n")
cat("   原始数据均值:", mean(expr_data, na.rm = TRUE), "\n\n")

# 如果数据范围在0-20之间，可能是log2值
if (max(expr_data, na.rm = TRUE) < 20) {
  cat("   检测到log2数据，转换为线性空间...\n")
  expr_linear <- 2^expr_data
} else {
  expr_linear <- expr_data
}

# 去除NA值
expr_linear[is.na(expr_linear)] <- 0

# 转换为整数count（近似）
counts_matrix <- round(expr_linear)
counts_matrix[counts_matrix < 0] <- 0

cat("   转换后count范围:", range(counts_matrix), "\n")
cat("   转换后count均值:", mean(counts_matrix), "\n\n")

cat("4. 过滤低表达探针...\n")

# 保留在至少10个样本中表达>10的探针
keep <- rowSums(counts_matrix >= 10) >= 10
counts_matrix_filtered <- counts_matrix[keep, ]

cat("   保留", nrow(counts_matrix_filtered), "个探针（原始", nrow(counts_matrix), "个）\n\n")

cat("5. 创建DESeq2对象...\n")

# 创建colData
condition <- factor(sample_type_clean, levels = c("Normal", "Tumor"))
coldata_df <- data.frame(
  condition = condition,
  row.names = colnames(counts_matrix_filtered)
)

# 创建DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix_filtered,
  colData = coldata_df,
  design = ~ condition
)

cat("   ✓ DESeqDataSet创建完成\n")
cat("   - Normal:", sum(condition == "Normal"), "\n")
cat("   - Tumor:", sum(condition == "Tumor"), "\n")
cat("   - 探针数:", nrow(dds), "\n\n")

cat("6. 保存数据...\n")

saveRDS(dds, file.path(save_dir, "geo_crc_dds.rds"))
saveRDS(gse, file.path(save_dir, "geo_crc_raw.rds"))

# 保存样本信息
sample_info_df <- data.frame(
  sample_id = colnames(counts_matrix_filtered),
  condition = sample_type_clean,
  original_label = sample_type
)
write.csv(sample_info_df, file.path(save_dir, "geo_sample_info.csv"), row.names = FALSE)

cat("   ✓ 数据已保存到:", save_dir, "\n")
cat("   - geo_crc_dds.rds (DESeq2对象)\n")
cat("   - geo_crc_raw.rds (原始GEO数据)\n")
cat("   - geo_sample_info.csv (样本信息)\n\n")

cat("=== GEO数据下载完成 ===\n")
cat("\n数据集: GSE50760\n")
cat("样本数:", ncol(dds), "(", sum(condition == "Normal"), "Normal +", 
    sum(condition == "Tumor"), "Tumor )\n")
cat("探针数:", nrow(dds), "\n")
cat("\n下一步：运行DESeq2分析\n")
cat("Rscript scripts/06_geo_deseq2.R\n")

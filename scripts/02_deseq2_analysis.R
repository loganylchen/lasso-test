#!/usr/bin/env Rscript
# CRC真实数据分析 - Step 2: DESeq2差异分析
# 
# 输入：crc_rse_raw.rds
# 输出：差异基因列表

cat("=== CRC数据 - DESeq2差异分析 ===\n\n")

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
})

save_dir <- "~/.openclaw/workspace/lasso-biomarker-tutorial/data/crc_analysis"

cat("1. 加载数据...\n")
rse <- readRDS(file.path(save_dir, "crc_rse_raw.rds"))
cat("   - 样本数:", ncol(rse), "\n")
cat("   - 基因数:", nrow(rse), "\n\n")

cat("2. 准备分组信息...\n")

# 检查样本信息
coldata <- colData(rse)

# 自动识别Normal vs Tumor分组
# 通常在TCGA数据中，样本类型在tcga.tcga_barcode或其他列
if ("tcga.gdc_cases.samples.sample_type" %in% colnames(coldata)) {
  sample_type_col <- "tcga.gdc_cases.samples.sample_type"
} else if ("tcga.cgc_sample_sample_type" %in% colnames(coldata)) {
  sample_type_col <- "tcga.cgc_sample_sample_type"
} else {
  # 手动指定（需要根据实际数据调整）
  cat("   ⚠️ 未找到标准TCGA样本类型列\n")
  cat("   可用列:", paste(head(colnames(coldata), 10), collapse = ", "), "...\n")
  cat("   请手动检查colData并指定分组列\n")
  
  # 尝试使用sample id模式识别
  # TCGA: 正常样本通常是 10A-19Z，肿瘤是 01A-09Z
  barcode_col <- colnames(coldata)[grepl("barcode", colnames(coldata), ignore.case = TRUE)][1]
  
  if (!is.na(barcode_col)) {
    barcodes <- coldata[[barcode_col]]
    # TCGA样本编码：TCGA-XX-XXXX-01A-11R-XXXX-XX
    #              位置13-14: 01-09=肿瘤，10-19=正常
    sample_types <- ifelse(grepl("-1[0-9][A-Z]-", barcodes), "Normal", "Tumor")
    coldata$condition <- factor(sample_types)
    sample_type_col <- "condition"
    cat("   根据TCGA barcode推断分组\n")
  } else {
    stop("无法自动识别样本分组，请手动指定")
  }
} else {
  coldata$condition <- factor(coldata[[sample_type_col]])
  sample_type_col <- "condition"
}

cat("   分组列:", sample_type_col, "\n")
cat("   分组分布:\n")
print(table(coldata[[sample_type_col]]))

# 筛选正常和肿瘤样本
normal_samples <- which(grepl("normal|solid tissue normal", 
                               coldata[[sample_type_col]], ignore.case = TRUE))
tumor_samples <- which(grepl("primary|tumor|metastatic", 
                              coldata[[sample_type_col]], ignore.case = TRUE))

cat("\n   正常样本:", length(normal_samples), "\n")
cat("   肿瘤样本:", length(tumor_samples), "\n\n")

if (length(normal_samples) < 3 || length(tumor_samples) < 3) {
  stop("样本数不足（每组至少需要3个）")
}

# 限制样本数（加快分析）
max_samples_per_group <- 30
if (length(normal_samples) > max_samples_per_group) {
  normal_samples <- sample(normal_samples, max_samples_per_group)
}
if (length(tumor_samples) > max_samples_per_group) {
  tumor_samples <- sample(tumor_samples, max_samples_per_group)
}

selected_samples <- c(normal_samples, tumor_samples)
rse_subset <- rse[, selected_samples]

cat("3. 准备DESeq2对象...\n")

# 创建简化的分组信息
condition <- factor(
  c(rep("Normal", length(normal_samples)), 
    rep("Tumor", length(tumor_samples))),
  levels = c("Normal", "Tumor")
)

# 创建DESeqDataSet
colData(rse_subset)$condition <- condition
dds <- DESeqDataSet(rse_subset, design = ~ condition)

cat("   ✓ DESeqDataSet创建完成\n")
cat("   - 样本数:", ncol(dds), "(", sum(condition == "Normal"), "Normal +", 
    sum(condition == "Tumor"), "Tumor )\n")
cat("   - 基因数:", nrow(dds), "\n\n")

cat("4. 预过滤低表达基因...\n")
# 保留至少10个样本中count >= 10的基因
keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep, ]
cat("   保留", nrow(dds), "个基因\n\n")

cat("5. 运行DESeq2（这可能需要几分钟）...\n")
dds <- DESeq(dds)
cat("   ✓ DESeq2分析完成\n\n")

cat("6. 提取差异基因...\n")
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
res <- res[order(res$padj), ]

# 统计显著差异基因
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
cat("   显著差异基因 (padj < 0.05, |log2FC| > 1):", nrow(sig_genes), "\n")
cat("   - 上调:", sum(sig_genes$log2FoldChange > 0), "\n")
cat("   - 下调:", sum(sig_genes$log2FoldChange < 0), "\n\n")

# 保存结果
saveRDS(dds, file.path(save_dir, "crc_dds.rds"))
saveRDS(res, file.path(save_dir, "crc_deseq_results.rds"))
write.csv(as.data.frame(sig_genes), 
          file.path(save_dir, "crc_sig_genes.csv"))

cat("7. 生成可视化...\n")

# MA plot
pdf(file.path(save_dir, "MA_plot.pdf"), width = 8, height = 6)
plotMA(res, ylim = c(-5, 5), main = "CRC DESeq2 MA Plot")
dev.off()

# Volcano plot
res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
                              "Significant", "Not Significant")

p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c("gray70", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "CRC Volcano Plot", x = "log2(Fold Change)", y = "-log10(padj)") +
  theme_minimal()

ggsave(file.path(save_dir, "volcano_plot.pdf"), p, width = 8, height = 6)

cat("   ✓ 可视化已保存\n\n")

cat("=== DESeq2分析完成 ===\n")
cat("\n输出文件:\n")
cat("  - crc_dds.rds (DESeq2对象)\n")
cat("  - crc_deseq_results.rds (完整结果)\n")
cat("  - crc_sig_genes.csv (显著差异基因)\n")
cat("  - MA_plot.pdf\n")
cat("  - volcano_plot.pdf\n")
cat("\n下一步：用LASSO从差异基因中筛选biomarker\n")

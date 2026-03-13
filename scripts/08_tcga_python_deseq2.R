#!/usr/bin/env Rscript
# TCGA Python数据 - DESeq2分析

cat("=== TCGA Python数据 - DESeq2分析 ===\n\n")

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
})

save_dir <- "~/.openclaw/workspace/lasso-biomarker-tutorial/data/tcga_coad_python"

cat("1. 加载DESeq2对象...\n")
dds <- readRDS(file.path(save_dir, "tcga_python_dds.rds"))

cat("   ✓ 数据加载完成\n")
cat("   - 样本数:", ncol(dds), "\n")
cat("   - 基因数:", nrow(dds), "\n\n")

cat("2. 运行DESeq2（真实TCGA数据，可能需要几分钟）...\n")
dds <- DESeq(dds)
cat("   ✓ DESeq2分析完成\n\n")

cat("3. 提取差异基因...\n")
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
res <- res[order(res$padj), ]

# 统计显著差异基因
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
cat("   显著差异基因 (padj < 0.05, |log2FC| > 1):", nrow(sig_genes), "\n")
cat("   - 上调:", sum(sig_genes$log2FoldChange > 0), "\n")
cat("   - 下调:", sum(sig_genes$log2FoldChange < 0), "\n\n")

# 显示top 20差异基因
cat("   Top 20差异基因:\n")
print(head(as.data.frame(sig_genes)[, c("log2FoldChange", "padj")], 20))
cat("\n")

# 保存结果
saveRDS(dds, file.path(save_dir, "tcga_python_dds_analyzed.rds"))
saveRDS(res, file.path(save_dir, "tcga_python_deseq_results.rds"))
write.csv(as.data.frame(sig_genes), 
          file.path(save_dir, "tcga_python_sig_genes.csv"))

cat("4. 生成可视化...\n")

# MA plot
pdf(file.path(save_dir, "tcga_python_MA_plot.pdf"), width = 8, height = 6)
plotMA(res, ylim = c(-10, 10), main = "TCGA-COAD (Python) DESeq2 MA Plot")
dev.off()

# Volcano plot
res_df <- as.data.frame(res)
res_df$significant <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
                              "Significant", "Not Significant")

p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.4, size = 0.5) +
  scale_color_manual(values = c("gray70", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  xlim(-10, 10) +
  labs(title = "TCGA-COAD (Python) Volcano Plot", 
       x = "log2(Fold Change)", 
       y = "-log10(adjusted p-value)") +
  theme_minimal()

ggsave(file.path(save_dir, "tcga_python_volcano_plot.pdf"), p, width = 10, height = 6)

# PCA plot
pdf(file.path(save_dir, "tcga_python_pca_plot.pdf"), width = 8, height = 6)
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition") +
  ggtitle("TCGA-COAD (Python) PCA Plot") +
  theme_minimal()
dev.off()

cat("   ✓ 可视化已保存\n\n")

cat("=== DESeq2分析完成 ===\n")
cat("\n输出文件:\n")
cat("  - tcga_python_dds_analyzed.rds\n")
cat("  - tcga_python_deseq_results.rds\n")
cat("  - tcga_python_sig_genes.csv (", nrow(sig_genes), "个)\n")
cat("  - tcga_python_MA_plot.pdf\n")
cat("  - tcga_python_volcano_plot.pdf\n")
cat("  - tcga_python_pca_plot.pdf\n")
cat("\n下一步：运行LASSO分析\n")
cat("Rscript scripts/09_tcga_python_lasso.R\n")

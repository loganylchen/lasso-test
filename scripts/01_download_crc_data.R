#!/usr/bin/env Rscript
# CRC真实数据分析 - Step 1: 下载recount3数据并DESeq2差异分析
# 
# 流程：
# 1. 用recount3下载CRC (结直肠癌) 数据
# 2. DESeq2找差异基因
# 3. 用LASSO筛选biomarker

cat("=== CRC数据分析 - 真实数据LASSO测试 ===\n\n")

# 检查并安装必要的包
cat("1. 检查并安装必要的R包...\n")

required_bioc <- c("recount3", "DESeq2", "SummarizedExperiment")
required_cran <- c("dplyr", "ggplot2")

for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("   安装Bioconductor包:", pkg, "\n")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, ask = FALSE)
  }
}

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("   安装CRAN包:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

cat("   ✓ 所有包已准备\n\n")

# 加载包
suppressPackageStartupMessages({
  library(recount3)
  library(DESeq2)
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
})

cat("2. 搜索recount3中的CRC数据集...\n")

# 搜索可用的CRC相关数据
human_projects <- available_projects()

# 筛选CRC相关研究
crc_projects <- human_projects %>%
  filter(grepl("colorectal|colon|rectal|CRC", project_title, ignore.case = TRUE) |
         grepl("colorectal|colon|rectal|CRC", project, ignore.case = TRUE))

if (nrow(crc_projects) == 0) {
  cat("   未找到CRC相关数据集，搜索TCGA-COAD...\n")
  crc_projects <- human_projects %>%
    filter(grepl("TCGA", project) & grepl("COAD|READ", project))
}

cat("   找到", nrow(crc_projects), "个CRC相关数据集\n")
print(head(crc_projects[, c("project", "organism", "file_source")], 3))

# 选择第一个数据集
project_info <- crc_projects[1, ]
cat("\n   选择数据集:", project_info$project, "\n")

cat("\n3. 下载数据（可能需要几分钟）...\n")

# 下载数据
rse <- create_rse(project_info)

cat("   ✓ 数据下载完成\n")
cat("   - 样本数:", ncol(rse), "\n")
cat("   - 基因数:", nrow(rse), "\n\n")

# 查看样本信息
cat("4. 检查样本分组信息...\n")

# 查看可用的表型列
coldata_cols <- colnames(colData(rse))
cat("   可用的表型列:", length(coldata_cols), "个\n")

# 尝试识别疾病状态列
possible_group_cols <- coldata_cols[grepl("disease|condition|tissue|type|sample", 
                                           coldata_cols, ignore.case = TRUE)]
cat("   可能的分组列:\n")
print(head(possible_group_cols, 10))

# 保存中间结果
save_dir <- "~/.openclaw/workspace/lasso-biomarker-tutorial/data/crc_analysis"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(rse, file.path(save_dir, "crc_rse_raw.rds"))
cat("\n   ✓ 原始数据已保存:", file.path(save_dir, "crc_rse_raw.rds"), "\n")

cat("\n=== 数据下载完成 ===\n")
cat("\n下一步：\n")
cat("1. 检查 colData(rse) 确定正常/疾病分组列\n")
cat("2. 创建 DESeq2 对象\n")
cat("3. 运行差异分析\n")
cat("4. 用LASSO筛选biomarker\n")

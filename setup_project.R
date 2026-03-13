#!/usr/bin/env Rscript
# 项目初始化脚本

cat("=== 创建LASSO Biomarker Tutorial项目结构 ===\n")

# 创建目录
dirs <- c(
  "data",
  "R",
  "notebooks",
  "tests/testthat",
  "output/figures",
  "output/results"
)

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat(paste0("✓ 创建目录: ", dir, "\n"))
  }
}

cat("\n=== 检查R包依赖 ===\n")

# 检查并安装必要的包
required_packages <- c(
  "glmnet",      # LASSO回归
  "ggplot2",     # 可视化
  "pROC",        # ROC曲线
  "pheatmap",    # 热图
  "testthat",    # 单元测试
  "caret",       # 混淆矩阵
  "tidyr",       # 数据整理
  "dplyr"        # 数据操作
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("安装包: ", pkg, "\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  } else {
    cat(paste0("✓ ", pkg, "\n"))
  }
}

cat("\n=== 项目初始化完成 ===\n")
cat("下一步: 运行 Rscript tests/testthat/test_lasso.R 开始TDD开发\n")

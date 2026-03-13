# LASSO Biomarker Selection Tutorial

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

使用LASSO回归从RNA-seq位点数据中筛选区分正常和疾病样本的生物标志物（biomarker）完整教程。

## 📋 项目简介

本项目提供了一个完整的、可复现的LASSO生物标志物筛选工作流，包括：

- ✅ 完整的R函数库（TDD开发）
- ✅ 交互式RMarkdown教程
- ✅ 单元测试
- ✅ 丰富的可视化
- ✅ 详细文档

## 🎯 适用场景

- RNA-seq差异表达分析
- RNA修饰位点筛选
- 蛋白组学biomarker发现
- 代谢组学特征选择
- 任何高维组学数据的特征筛选

## 📁 项目结构

```
lasso-biomarker-tutorial/
├── README.md                    # 本文件
├── data/
│   ├── generate_data.R          # 示例数据生成脚本
│   └── example_data.rds         # 示例数据（运行后生成）
├── R/
│   ├── lasso_biomarker.R        # LASSO核心函数
│   ├── visualization.R          # 可视化函数
│   └── evaluation.R             # 模型评估函数
├── notebooks/
│   └── tutorial.Rmd             # 完整教程（推荐从这里开始）
├── tests/
│   └── testthat/
│       └── test_lasso.R         # 单元测试
└── output/                      # 输出目录
    ├── figures/                 # 图表
    └── results/                 # 结果文件
```

## 🚀 快速开始

### 1. 安装依赖

```r
# 安装必要的R包
install.packages(c(
  "glmnet",      # LASSO回归
  "ggplot2",     # 可视化
  "pROC",        # ROC曲线
  "pheatmap",    # 热图
  "caret",       # 混淆矩阵
  "testthat"     # 单元测试
))
```

### 2. 克隆项目

```bash
git clone https://github.com/yourusername/lasso-biomarker-tutorial.git
cd lasso-biomarker-tutorial
```

### 3. 运行教程

**方法1：RMarkdown（推荐）**
```r
# 打开RStudio
# 打开 notebooks/tutorial.Rmd
# 点击 "Knit" 生成HTML报告
```

**方法2：R脚本**
```r
# 加载函数
source("R/lasso_biomarker.R")
source("R/visualization.R")
source("R/evaluation.R")
source("data/generate_data.R")

# 生成示例数据
data <- generate_example_data(n_samples = 20, n_features = 500)

# 运行LASSO
cv_fit <- run_lasso_cv(data$X, data$y)

# 筛选biomarker
biomarkers <- select_biomarkers(cv_fit, data$X)

# 评估模型
eval_result <- evaluate_model(cv_fit, data$X, data$y)

# 可视化
plot_roc_curve(cv_fit, data$X, data$y)
plot_biomarker_heatmap(data$X, data$y, biomarkers)
```

## 📊 核心功能

### 1. 数据生成
```r
generate_example_data(
  n_samples = 20,       # 样本数
  n_features = 500,     # 特征数
  n_informative = 20,   # 真实biomarker数
  seed = 42
)
```

### 2. LASSO训练
```r
run_lasso_cv(
  X = feature_matrix,   # 特征矩阵（features × samples）
  y = labels,           # 标签向量（0/1）
  nfolds = 5,           # 交叉验证折数
  alpha = 1             # LASSO参数
)
```

### 3. Biomarker筛选
```r
select_biomarkers(
  cv_fit = cv_fit,
  X = feature_matrix,
  lambda = "lambda.1se"  # 或 "lambda.min"
)
```

### 4. 模型评估
```r
evaluate_model(
  cv_fit = cv_fit,
  X = feature_matrix,
  y = labels,
  lambda = "lambda.1se"
)
# 返回: AUC, 准确率, 灵敏度, 特异度, 混淆矩阵
```

### 5. 可视化
```r
# ROC曲线
plot_roc_curve(cv_fit, X, y)

# 系数路径图
plot_coefficient_path(cv_fit)

# Biomarker热图
plot_biomarker_heatmap(X, y, biomarkers)

# Top biomarker箱线图
plot_biomarker_boxplot(X, y, biomarkers, top_n = 10)

# 交叉验证曲线
plot_cv_curve(cv_fit)
```

## 🧪 运行测试

```r
library(testthat)
test_file("tests/testthat/test_lasso.R")
```

所有测试应该通过 ✅

## 📈 示例结果

使用20个样本、500个特征（其中20个真实biomarker）的模拟数据：

- **AUC**: 0.95+
- **准确率**: 90%+
- **Biomarker召回率**: 70-90%

## 💡 使用建议

### Lambda选择
- **lambda.min**: AUC最大，特征较多，适合预测
- **lambda.1se**: 更保守，特征较少，适合生物学解释（**推荐**）

### 样本量要求
- 最少：20个样本（10 vs 10）
- 推荐：50+个样本（25 vs 25）
- 理想：100+个样本（50 vs 50）

### 数据预处理
1. 标准化（z-score）
2. 去除低方差特征
3. 批次效应校正（如有多批次）
4. 缺失值处理

## 📚 相关资源

### 理论基础
- Tibshirani, R. (1996). "Regression Shrinkage and Selection via the Lasso". *JRSS-B*.
- Friedman, J., Hastie, T., & Tibshirani, R. (2010). "Regularization Paths for Generalized Linear Models via Coordinate Descent". *JSS*.

### 工具文档
- [glmnet官方文档](https://glmnet.stanford.edu/)
- [Statistical Learning (Hastie & Tibshirani)](https://www.statlearning.com/)

## 🤝 贡献

欢迎提交Issue和Pull Request！

## 📄 License

MIT License

## ✉️ 联系

如有问题，请提交Issue或联系：[your.email@example.com]

---

**⭐ 如果这个项目对你有帮助，请给个星标！**

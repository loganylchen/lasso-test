# LASSO Biomarker Selection Tutorial

## 项目目标
使用LASSO回归从RNA-seq位点数据中筛选区分正常和疾病样本的生物标志物。

## 技术栈
- R (>= 4.0)
- glmnet - LASSO回归
- ggplot2 - 可视化
- pROC - ROC曲线
- pheatmap - 热图
- caret - 交叉验证

## 数据规格
- 样本数：20（假设10正常 + 10疾病）
- 特征：RNA修饰位点（数量待定，假设100-1000个位点）
- 数据格式：位点 × 样本矩阵

## 输出内容
1. RMarkdown教程（tutorial.Rmd）
2. R脚本（lasso_biomarker.R）
3. 示例数据生成脚本
4. 可视化函数集
5. 完整测试用例
6. README文档

## 项目结构
```
lasso-biomarker-tutorial/
├── README.md
├── renv.lock              # R环境管理
├── data/
│   ├── example_data.rds   # 示例数据
│   └── generate_data.R    # 数据生成
├── R/
│   ├── lasso_biomarker.R  # 核心函数
│   ├── visualization.R    # 可视化
│   └── evaluation.R       # 模型评估
├── notebooks/
│   └── tutorial.Rmd       # 主教程
├── tests/
│   └── test_lasso.R       # 单元测试
└── output/                # 输出图表
    ├── figures/
    └── results/
```

## 分析流程
1. 数据加载与预处理
2. 特征标准化
3. LASSO回归训练（交叉验证选lambda）
4. Biomarker筛选（非零系数）
5. 模型评估（ROC、AUC、混淆矩阵）
6. 可视化（系数路径图、ROC曲线、热图、箱线图）
7. 结果解释

## 关键参数
- alpha = 1（LASSO）
- cv.folds = 5（5折交叉验证）
- family = "binomial"（二分类）

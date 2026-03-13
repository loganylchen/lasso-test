# LASSO Biomarker Selection - 测试用例
# TDD: 先写测试，再实现功能

library(testthat)
library(glmnet)

# 测试数据生成
test_that("生成示例数据函数工作正常", {
  source("../../data/generate_data.R")
  
  data <- generate_example_data(n_samples = 20, n_features = 100, seed = 42)
  
  expect_equal(nrow(data$X), 100)  # 100个特征
  expect_equal(ncol(data$X), 20)   # 20个样本
  expect_equal(length(data$y), 20) # 20个标签
  expect_true(all(data$y %in% c(0, 1)))  # 标签为0/1
  expect_equal(sum(data$y), 10)    # 10个疾病样本
})

# 测试LASSO交叉验证
test_that("LASSO交叉验证函数正确运行", {
  source("../../R/lasso_biomarker.R")
  source("../../data/generate_data.R")
  
  data <- generate_example_data(n_samples = 20, n_features = 100, seed = 42)
  
  result <- run_lasso_cv(data$X, data$y, nfolds = 5, alpha = 1)
  
  expect_s3_class(result, "cv.glmnet")
  expect_true(!is.null(result$lambda.min))
  expect_true(!is.null(result$lambda.1se))
})

# 测试Biomarker筛选
test_that("Biomarker筛选函数返回正确结果", {
  source("../../R/lasso_biomarker.R")
  source("../../data/generate_data.R")
  
  data <- generate_example_data(n_samples = 20, n_features = 100, seed = 42)
  cv_fit <- run_lasso_cv(data$X, data$y)
  
  biomarkers <- select_biomarkers(cv_fit, data$X, lambda = "lambda.1se")
  
  expect_type(biomarkers, "list")
  expect_true("selected_features" %in% names(biomarkers))
  expect_true("coefficients" %in% names(biomarkers))
  expect_true(length(biomarkers$selected_features) > 0)
})

# 测试模型评估
test_that("模型评估函数计算正确指标", {
  source("../../R/evaluation.R")
  source("../../R/lasso_biomarker.R")
  source("../../data/generate_data.R")
  
  data <- generate_example_data(n_samples = 20, n_features = 100, seed = 42)
  cv_fit <- run_lasso_cv(data$X, data$y)
  
  eval_result <- evaluate_model(cv_fit, data$X, data$y, lambda = "lambda.1se")
  
  expect_type(eval_result, "list")
  expect_true("auc" %in% names(eval_result))
  expect_true("accuracy" %in% names(eval_result))
  expect_true("confusion_matrix" %in% names(eval_result))
  expect_true(eval_result$auc >= 0 && eval_result$auc <= 1)
})

# 测试可视化函数
test_that("可视化函数生成图表", {
  source("../../R/visualization.R")
  source("../../R/lasso_biomarker.R")
  source("../../data/generate_data.R")
  
  data <- generate_example_data(n_samples = 20, n_features = 100, seed = 42)
  cv_fit <- run_lasso_cv(data$X, data$y)
  
  # 系数路径图
  p1 <- plot_coefficient_path(cv_fit)
  expect_s3_class(p1, "ggplot")
  
  # ROC曲线
  p2 <- plot_roc_curve(cv_fit, data$X, data$y)
  expect_s3_class(p2, "ggplot")
})

cat("✅ 所有测试用例定义完成\n")

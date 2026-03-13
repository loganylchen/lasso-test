# LASSO Biomarker筛选核心函数
library(glmnet)

#' LASSO回归交叉验证
#' 
#' @param X 特征矩阵（features × samples）
#' @param y 标签向量（0/1）
#' @param nfolds 交叉验证折数（默认5）
#' @param alpha LASSO参数（1=LASSO, 0=Ridge, 0-1=Elastic Net）
#' @param family 模型类型（"binomial"用于二分类）
#' @return cv.glmnet对象
#' @export
run_lasso_cv <- function(X, y, nfolds = 5, alpha = 1, family = "binomial") {
  # 转置矩阵：glmnet需要样本×特征格式
  X_t <- t(X)
  
  # 运行交叉验证
  cv_fit <- cv.glmnet(
    x = X_t,
    y = y,
    family = family,
    alpha = alpha,
    nfolds = nfolds,
    type.measure = "auc"  # 使用AUC作为评估指标
  )
  
  return(cv_fit)
}


#' 从LASSO模型中筛选Biomarker
#' 
#' @param cv_fit cv.glmnet对象
#' @param X 原始特征矩阵
#' @param lambda lambda选择（"lambda.min" 或 "lambda.1se"）
#' @return list包含selected_features（特征名）和coefficients（系数）
#' @export
select_biomarkers <- function(cv_fit, X, lambda = "lambda.1se") {
  # 获取lambda值
  lambda_value <- if (lambda == "lambda.min") {
    cv_fit$lambda.min
  } else {
    cv_fit$lambda.1se
  }
  
  # 提取系数
  coef_matrix <- coef(cv_fit, s = lambda_value)
  coef_vector <- as.vector(coef_matrix)
  names(coef_vector) <- rownames(coef_matrix)
  
  # 筛选非零系数（排除截距）
  nonzero_idx <- which(coef_vector[-1] != 0)  # 去掉第一个（截距）
  selected_features <- rownames(X)[nonzero_idx]
  selected_coefs <- coef_vector[-1][nonzero_idx]
  
  return(list(
    selected_features = selected_features,
    coefficients = selected_coefs,
    n_features = length(selected_features),
    lambda = lambda_value
  ))
}


#' 预测新样本
#' 
#' @param cv_fit cv.glmnet对象
#' @param X_new 新样本特征矩阵
#' @param lambda lambda选择
#' @param type 预测类型（"response"返回概率，"class"返回类别）
#' @return 预测结果
#' @export
predict_biomarker <- function(cv_fit, X_new, lambda = "lambda.1se", type = "response") {
  lambda_value <- if (lambda == "lambda.min") {
    cv_fit$lambda.min
  } else {
    cv_fit$lambda.1se
  }
  
  X_new_t <- t(X_new)
  predictions <- predict(cv_fit, newx = X_new_t, s = lambda_value, type = type)
  
  return(as.vector(predictions))
}


#' 提取LASSO路径上的所有系数
#' 
#' @param cv_fit cv.glmnet对象
#' @return data.frame包含lambda、特征、系数
#' @export
extract_coefficient_path <- function(cv_fit) {
  coef_matrix <- coef(cv_fit$glmnet.fit)
  
  # 转换为data.frame
  path_df <- data.frame(
    lambda = cv_fit$glmnet.fit$lambda,
    n_nonzero = cv_fit$glmnet.fit$df
  )
  
  return(path_df)
}

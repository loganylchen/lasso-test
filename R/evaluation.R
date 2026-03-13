# 模型评估函数
library(pROC)
library(caret)

#' 评估LASSO模型性能
#' 
#' @param cv_fit cv.glmnet对象
#' @param X 特征矩阵
#' @param y 真实标签
#' @param lambda lambda选择
#' @return list包含AUC、准确率、混淆矩阵等指标
#' @export
evaluate_model <- function(cv_fit, X, y, lambda = "lambda.1se") {
  # 转置矩阵
  X_t <- t(X)
  
  # 获取lambda值
  lambda_value <- if (lambda == "lambda.min") {
    cv_fit$lambda.min
  } else {
    cv_fit$lambda.1se
  }
  
  # 预测概率
  pred_prob <- predict(cv_fit, newx = X_t, s = lambda_value, type = "response")
  pred_prob <- as.vector(pred_prob)
  
  # 预测类别（阈值0.5）
  pred_class <- ifelse(pred_prob > 0.5, 1, 0)
  
  # 计算ROC和AUC
  roc_obj <- roc(y, pred_prob, quiet = TRUE)
  auc_value <- as.numeric(auc(roc_obj))
  
  # 混淆矩阵
  cm <- confusionMatrix(
    factor(pred_class, levels = c(0, 1)),
    factor(y, levels = c(0, 1))
  )
  
  # 计算其他指标
  accuracy <- cm$overall["Accuracy"]
  sensitivity <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]
  
  return(list(
    auc = auc_value,
    accuracy = as.numeric(accuracy),
    sensitivity = as.numeric(sensitivity),
    specificity = as.numeric(specificity),
    confusion_matrix = cm$table,
    roc_object = roc_obj,
    predictions = data.frame(
      true_label = y,
      pred_prob = pred_prob,
      pred_class = pred_class
    )
  ))
}


#' 计算特征重要性得分
#' 
#' @param biomarkers select_biomarkers()的返回值
#' @return data.frame包含特征名和重要性得分
#' @export
calculate_feature_importance <- function(biomarkers) {
  importance_df <- data.frame(
    feature = biomarkers$selected_features,
    coefficient = biomarkers$coefficients,
    abs_coefficient = abs(biomarkers$coefficients)
  )
  
  # 按绝对系数排序
  importance_df <- importance_df[order(-importance_df$abs_coefficient), ]
  rownames(importance_df) <- NULL
  
  return(importance_df)
}


#' 交叉验证性能评估
#' 
#' @param X 特征矩阵
#' @param y 标签向量
#' @param nfolds 折数
#' @return data.frame包含每折的性能指标
#' @export
cross_validate_performance <- function(X, y, nfolds = 5) {
  n_samples <- length(y)
  fold_indices <- sample(rep(1:nfolds, length.out = n_samples))
  
  results <- data.frame(
    fold = integer(),
    auc = numeric(),
    accuracy = numeric(),
    n_biomarkers = integer()
  )
  
  for (fold in 1:nfolds) {
    # 分割数据
    test_idx <- which(fold_indices == fold)
    train_idx <- which(fold_indices != fold)
    
    X_train <- X[, train_idx]
    y_train <- y[train_idx]
    X_test <- X[, test_idx]
    y_test <- y[test_idx]
    
    # 训练模型
    cv_fit <- run_lasso_cv(X_train, y_train, nfolds = nfolds - 1)
    
    # 评估
    eval_result <- evaluate_model(cv_fit, X_test, y_test)
    biomarkers <- select_biomarkers(cv_fit, X_train)
    
    results <- rbind(results, data.frame(
      fold = fold,
      auc = eval_result$auc,
      accuracy = eval_result$accuracy,
      n_biomarkers = biomarkers$n_features
    ))
  }
  
  return(results)
}

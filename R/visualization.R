# 可视化函数
library(ggplot2)
library(pheatmap)
library(pROC)

#' 绘制LASSO系数路径图
#' 
#' @param cv_fit cv.glmnet对象
#' @return ggplot对象
#' @export
plot_coefficient_path <- function(cv_fit) {
  # 提取系数矩阵
  coef_matrix <- as.matrix(coef(cv_fit$glmnet.fit))[-1, ]  # 去掉截距
  lambda_seq <- cv_fit$glmnet.fit$lambda
  
  # 转换为long格式
  df_list <- list()
  for (i in 1:nrow(coef_matrix)) {
    df_list[[i]] <- data.frame(
      lambda = lambda_seq,
      coefficient = coef_matrix[i, ],
      feature = rownames(coef_matrix)[i]
    )
  }
  plot_df <- do.call(rbind, df_list)
  
  # 只绘制非零系数的特征
  nonzero_features <- unique(plot_df$feature[plot_df$coefficient != 0])
  plot_df <- plot_df[plot_df$feature %in% nonzero_features, ]
  
  # 绘图
  p <- ggplot(plot_df, aes(x = log(lambda), y = coefficient, color = feature)) +
    geom_line(alpha = 0.7) +
    geom_vline(xintercept = log(cv_fit$lambda.min), linetype = "dashed", color = "red") +
    geom_vline(xintercept = log(cv_fit$lambda.1se), linetype = "dashed", color = "blue") +
    labs(
      title = "LASSO Coefficient Path",
      x = "log(λ)",
      y = "Coefficient",
      caption = "红线=lambda.min, 蓝线=lambda.1se"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(p)
}


#' 绘制ROC曲线
#' 
#' @param cv_fit cv.glmnet对象
#' @param X 特征矩阵
#' @param y 真实标签
#' @param lambda lambda选择
#' @return ggplot对象
#' @export
plot_roc_curve <- function(cv_fit, X, y, lambda = "lambda.1se") {
  # 获取评估结果
  eval_result <- evaluate_model(cv_fit, X, y, lambda)
  roc_obj <- eval_result$roc_object
  
  # 提取ROC数据
  roc_df <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
  
  # 绘图
  p <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line(color = "blue", size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    labs(
      title = paste0("ROC Curve (AUC = ", round(eval_result$auc, 3), ")"),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme_minimal() +
    coord_fixed()
  
  return(p)
}


#' 绘制Biomarker热图
#' 
#' @param X 特征矩阵
#' @param y 标签向量
#' @param biomarkers select_biomarkers()的返回值
#' @param output_file 输出文件路径（可选）
#' @export
plot_biomarker_heatmap <- function(X, y, biomarkers, output_file = NULL) {
  # 提取选中的biomarker数据
  selected_X <- X[biomarkers$selected_features, , drop = FALSE]
  
  # 标准化（z-score）
  selected_X_scaled <- t(scale(t(selected_X)))
  
  # 创建列注释（样本分组）
  annotation_col <- data.frame(
    Group = factor(y, levels = c(0, 1), labels = c("Normal", "Disease"))
  )
  rownames(annotation_col) <- colnames(selected_X_scaled)
  
  # 颜色设置
  annotation_colors <- list(
    Group = c(Normal = "#4DAF4A", Disease = "#E41A1C")
  )
  
  # 绘制热图
  pheatmap(
    selected_X_scaled,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = paste0("Biomarker Heatmap (", biomarkers$n_features, " features)"),
    filename = output_file
  )
}


#' 绘制Biomarker箱线图
#' 
#' @param X 特征矩阵
#' @param y 标签向量
#' @param biomarkers select_biomarkers()的返回值
#' @param top_n 展示前N个重要的biomarker
#' @return ggplot对象
#' @export
plot_biomarker_boxplot <- function(X, y, biomarkers, top_n = 10) {
  # 选择前N个重要特征
  importance <- calculate_feature_importance(biomarkers)
  top_features <- head(importance$feature, top_n)
  
  # 准备数据
  plot_data <- data.frame()
  for (feature in top_features) {
    feature_values <- X[feature, ]
    plot_data <- rbind(plot_data, data.frame(
      Feature = feature,
      Value = feature_values,
      Group = factor(y, levels = c(0, 1), labels = c("Normal", "Disease"))
    ))
  }
  
  # 绘图
  p <- ggplot(plot_data, aes(x = Feature, y = Value, fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = c("Normal" = "#4DAF4A", "Disease" = "#E41A1C")) +
    labs(
      title = paste0("Top ", top_n, " Biomarkers"),
      x = "Biomarker",
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}


#' 绘制交叉验证曲线
#' 
#' @param cv_fit cv.glmnet对象
#' @return ggplot对象
#' @export
plot_cv_curve <- function(cv_fit) {
  plot_df <- data.frame(
    lambda = cv_fit$lambda,
    cvm = cv_fit$cvm,
    cvsd = cv_fit$cvsd
  )
  
  p <- ggplot(plot_df, aes(x = log(lambda), y = cvm)) +
    geom_line(color = "blue") +
    geom_ribbon(aes(ymin = cvm - cvsd, ymax = cvm + cvsd), alpha = 0.2, fill = "blue") +
    geom_vline(xintercept = log(cv_fit$lambda.min), linetype = "dashed", color = "red") +
    geom_vline(xintercept = log(cv_fit$lambda.1se), linetype = "dashed", color = "darkgreen") +
    labs(
      title = "Cross-Validation Curve",
      x = "log(λ)",
      y = "AUC",
      caption = "红线=lambda.min, 绿线=lambda.1se"
    ) +
    theme_minimal()
  
  return(p)
}

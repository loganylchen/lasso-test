# 生成LASSO Biomarker筛选的示例数据
# 模拟RNA-seq位点数据，用于区分正常和疾病样本

#' 生成示例数据
#' 
#' @param n_samples 样本数量（默认20）
#' @param n_features 特征数量（位点数，默认500）
#' @param n_informative 真正有信息的特征数（默认20）
#' @param seed 随机种子
#' @return list包含X（特征矩阵）、y（标签向量）、true_features（真实biomarker位点）
#' @export
generate_example_data <- function(n_samples = 20, 
                                   n_features = 500, 
                                   n_informative = 20,
                                   seed = 42) {
  set.seed(seed)
  
  # 样本分组：前一半正常，后一半疾病
  n_normal <- n_samples %/% 2
  n_disease <- n_samples - n_normal
  
  # 生成标签（0=正常，1=疾病）
  y <- c(rep(0, n_normal), rep(1, n_disease))
  
  # 初始化特征矩阵（特征 × 样本）
  X <- matrix(NA, nrow = n_features, ncol = n_samples)
  rownames(X) <- paste0("site_", sprintf("%03d", 1:n_features))
  colnames(X) <- paste0("sample_", sprintf("%02d", 1:n_samples))
  
  # 随机选择真正的biomarker位点
  true_biomarker_idx <- sort(sample(1:n_features, n_informative))
  
  # 生成数据
  for (i in 1:n_features) {
    if (i %in% true_biomarker_idx) {
      # 真实biomarker：正常和疾病样本有显著差异
      # 正常样本：均值5，标准差1
      X[i, y == 0] <- rnorm(n_normal, mean = 5, sd = 1)
      # 疾病样本：均值8，标准差1.5（差异+噪音）
      X[i, y == 1] <- rnorm(n_disease, mean = 8, sd = 1.5)
    } else {
      # 噪音特征：两组无差异
      X[i, ] <- rnorm(n_samples, mean = 6, sd = 2)
    }
  }
  
  return(list(
    X = X,
    y = y,
    true_features = rownames(X)[true_biomarker_idx],
    feature_names = rownames(X),
    sample_names = colnames(X)
  ))
}


#' 保存示例数据到RDS
#' 
#' @param output_path 输出路径
#' @export
save_example_data <- function(output_path = "data/example_data.rds") {
  data <- generate_example_data(
    n_samples = 20,
    n_features = 500,
    n_informative = 20,
    seed = 42
  )
  
  saveRDS(data, output_path)
  cat("✅ 示例数据已保存到:", output_path, "\n")
  cat("   - 样本数:", length(data$y), "\n")
  cat("   - 特征数:", nrow(data$X), "\n")
  cat("   - 真实biomarker数:", length(data$true_features), "\n")
  
  invisible(data)
}


# 如果直接运行此脚本，生成并保存数据
if (sys.nframe() == 0) {
  save_example_data()
}

#' grouping similarity
#' @description \code{grouping_similarity} calculates the similarity of two (or more) groupings (e.g. respondents). 
#' This function can be used, for example, to check how similar different clustering solutions are. 
#' @param gr1 a vector or a data frame (each columns should indicate different groupings) 
#' @param gr2 a vector to compare \code{gr1}, Default: NULL
#' @param groupwise.similarity Should we also show the similarity by each group?, Default: FALSE
#' @return \describe{
#' \item{wholelevel.similarity}{A matrix showing the similarity of two groupings (ranging from 0 to 1).}
#' \item{groupwise.similarity}{Similarity of each group in gr1 (if the input is the vector) or the first column 
#' (if the input is a data frame) and the rest of the groupings.}} 
#' @details Based on the groupings in \code{gr1} and \code{gr2},
#' \code{grouping_similarity} first construct the matrices in that the entry 
#' indicates 1 if i and j belong to the same group and 0 otherwise. 
#' It then calculate the standardized hamming distance among matrices
#'  (i.e., the realized hamming distance divided by the maximum hamming distance).
#' @examples
#' ## EXAMPLE1: input two vectors
#' gr1 <- c(rep(1, 3), rep(2, 3), rep(3, 4))
#' gr2 <- c(rep(1, 2), rep(2, 4), rep(3, 4))
#' grouping_similarity(gr1 = gr1, gr2 = gr2)
#' 
#' ## EXAMPLE2: input data.frame
#' gr3 <- data.frame(gr1, gr2)
#' grouping_similarity(gr1 = gr3, groupwise.similarity = TRUE)
#'
#' @export
grouping_similarity <- function(
    gr1, gr2 = NULL, groupwise.similarity = FALSE){
  # 内部で使用する関数
  .same_or_diff <- function(vec){
    n <- length(vec)
    mat <- matrix(0, n, n)
    for(i in 1:n){for(j in 1:n){
      if(vec[i] == vec[j]){
        mat[i,j] <- 1
      }
    }}
    diag(mat) <- 0
    mat
  }
  
  # grを一つに統合する
  gr <- gr1
  if(is.null(gr2)==F){
    gr <- cbind(gr, gr2)  
  }
  # colnameをストックしておく
  gr_colnames <- colnames(gr)
  n <- nrow(gr)
  
  # 同じグループに所属するものを1としたマトリクスをベクトル化
  same_mat_vectorized <- apply(X = gr, MARGIN = 2, FUN = .same_or_diff)
  # マトリクスに戻す
  same_mat_array <- array(same_mat_vectorized, dim = c(n, n, ncol(same_mat_vectorized)))
  
  # ハミング距離
  m <- dim(same_mat_array)[3]
  ham_dist <- matrix(0, m, m)
  for(i in 1:m){for(j in 1:m){
    ham_dist[i,j] <- sum(abs(same_mat_array[,,i]-same_mat_array[,,j]))
  }}
  whole_diff_ratio <- (n*(n-1)-ham_dist)/(n*(n-1))
  dimnames(whole_diff_ratio) <- list(gr_colnames, gr_colnames)
  
  if(groupwise.similarity==FALSE){
    out <- whole_diff_ratio
  }else{
    # グループ単位の類似性---------
    # グループ1のクラス
    gr1_class <- unique(gr[,1])
    # グループ1の各クラスと他のグルーピングとの類似性
    groupwise.ratio <- matrix(0, nrow = length(gr1_class), ncol = m)
    dimnames(groupwise.ratio) <- list(gr1_class, gr_colnames)
    for(i in 1:length(gr1_class)){for(j in 1:m){
      focal_class <- gr[,1]==gr1_class[i]
      focal_class_n <- sum(focal_class)
      focal_ham_dist <- sum(abs(same_mat_array[focal_class, focal_class, 1] - same_mat_array[focal_class, focal_class, j]))
      focal_diff_ratio <- (focal_class_n*(focal_class_n-1)-focal_ham_dist)/(focal_class_n*(focal_class_n-1))
      groupwise.ratio[i,j] <- focal_diff_ratio
    }}
    out <- list(wholelevel.similarity = whole_diff_ratio,
                groupwise.similarity = groupwise.ratio)
  }
  out
}


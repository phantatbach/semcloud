# semcloud: Post-processing of token-level clouds.
# Copyright (C) 2021 Mariana Montes
#
# See full notice in README.md file.

#' Log-transform distance matrices
#'
#' @param mat Distance matrix to transform
#' @param asDist Boolean: whether to force the matrix into a symmetric distance matrix.
#'
#' @return Double logged transformation of ranked rows of the input matrix,
#'  as \code{matrix}.
#' @export
transformMats <- function(mat, asDist = TRUE){
  ranked <- log(1+ log(t(apply(mat, 1, rank))))
  if (asDist) {
    return(as.matrix(stats::dist(ranked, diag=T, upper = T)))
  } else {
    return(ranked)
  }
}

#' Turn tibble into matrix
#'
#' @param df Tibble, for example, output of \code{\link{loadCloud}} with distances
#' @param rownames Name of the column that will become the row names. Default is
#'   "X_model", as it would work with the distances tibble.
#'
#' @return Matrix with right row names
#' @export
matricizeCloud <- function(df, rownames = "X_model"){
  df %>% data.frame(row.names = rownames) %>% as.matrix()
}

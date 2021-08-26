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
#' @param df Tibble
#' @param rownames Name of the column that will become the row names. Default is
#'   "X_model", as it would work with the distances tibble.
#'
#' @return Matrix with right row names
#' @export
matricizeCloud <- function(df, rownames = "X_model"){
  df %>% data.frame(row.names = rownames) %>% as.matrix()
}

#' Get Distance Matrix
#'
#' Store distances as `tsv` and return `dist` object.
#'
#' @param wwmx Distance matrix as outputted from \code{\link{focdistsFromCsv}}.
#' @param source_dir Directory to store the distance matrix in
#' @param lemma Name of the lemma, for the filename
#' @param suffix Suffix to add after `[lemma].models.dist` in the filename
#'
#' @return Distance matrix as a \code{dist} object.
#' @export
getDistMtx <- function(wwmx, source_dir, lemma, suffix = "") {
  # Save distance matrix as tsv
  suffix <- if (nchar(suffix) == 0) "tsv" else paste0(suffix, ".tsv")
  filename <- file.path(source_dir, paste(lemma, "models.dist", suffix, sep= "."))
  as.data.frame(wwmx) %>%
    dplyr::mutate(`_model` = row.names(wwmx)) %>%
    readr::write_tsv(filename)
  print(sprintf("Distance matrix saved in %s", filename))
  return(stats::as.dist(wwmx))
}


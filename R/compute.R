# semcloud: Post-processing of token-level clouds.
# Copyright (C) 2021 Mariana Montes
#
# See full notice in README.md file.

#' Dimensionality reduction for visualization
#'
#' Wrapper around \code{\link[vegan]{metaMDS}} and \code{\link[Rtsne]{Rtsne}}
#' for dimensionality reduction with certain set parameters.
#'
#' \code{\link[vegan]{metaMDS}} is run with \code{trymax=20} and \code{trace=FALSE}
#'  by default, so that the search for best solution is not printed.
#'  \code{\link[Rtsne]{Rtsne}} uses as default parameters
#'  \code{theta=0.0, check.duplicates = FALSE, max_iter = 1000}.
#'
#' @param d distance matrix as a \code{matrix}
#' @param dim number of dimensions, by default 2. As argument \code{k} for
#'   \code{\link[vegan]{metaMDS}} and \code{dims} for \code{\link[Rtsne]{Rtsne}}
#' @param technique either "mds" to run \code{\link[vegan]{metaMDS}},
#'   "tsne" to run \code{\link[Rtsne]{Rtsne}} or "umap" to run \code{\link[umap]{umap}}
#' @param perp perplexity value for \code{\link[Rtsne]{Rtsne}}, default is 30.
#'   This value is ignored when \code{\link[vegan]{metaMDS}} is run.
#' @param seed seed to keep randomness at check
#'
#' @return output of either \code{\link[vegan]{metaMDS}} or \code{\link[Rtsne]{Rtsne}}.
#' @export
getFit <- function(d, dim = 2, technique, perp = 30, seed = 8541){
  set.seed(seed)
  dst <- stats::as.dist(d)
  if (technique == 'mds') {
    vegan::metaMDS(dst, k=dim, trymax=20, trace=FALSE)
  } else if (technique == "tsne") {
    Rtsne::Rtsne(dst, dims=dim, perplexity=perp,
                        theta=0.0, check.duplicates = FALSE,
                        max_iter = 1000, is_distance = TRUE)
  } else if (technique == "umap") {
    umap::umap(d, input="dist")
  } else {
    stop("`technique` must be 'mds', 'tsne' or 'umap'.")
  }
}

#' Extract coordinates from fit
#'
#' @param fit Dimensionality reduction result; output from \code{\link{getFit}}
#' @param modelname Name of the model, for the names of the columns
#' @param rownames List of names for the rows
#' @param d Output from previous run or, if it's the first run, empty string.
#' @param source Technique used in dimensionality reduction: either "mds"
#'    for output from \code{\link[vegan]{metaMDS}}, "tsne" for \code{\link[Rtsne]{Rtsne}},
#'    "umap" for \code{\link[umap]{umap}}.
#'    (it matches `technique` from \code{\link{getFit}}).
#'
#' @return a [tibble][tibble::tibble-package] with the coordinates of each element.
#' @export
#'
#' @importFrom rlang .data
getCoords<-function(fit, modelname, rownames, d = "", source = "tsne"){
  if (source == "mds") { #output of metaMDS
    x = fit$points[,"MDS1"]
    y = fit$points[,"MDS2"]
  } else if (source == "tsne") { #output of Rtsne
    x = fit$Y[,1]
    y = fit$Y[,2]
  } else if (source == "umap") {
    x = fit$layout[,1]
    y = fit$layout[,2]
  } else {
    stop("`source` must be 'mds', 'tsne' or 'umap'.")
  }
  df <- tibble::tibble(
    `_id` = rownames,
    !!rlang::sym(paste0(modelname, ".x")) := x,
    !!rlang::sym(paste0(modelname, ".y")) := y)
  if (!is.character(d)) {
    df  <- df %>%
      dplyr::full_join(dplyr::select(d, !dplyr::starts_with(modelname)), by="_id")
    df <- replace(df, is.na(df), 0.0)
  }
  return(df)
}


#' Procrustes between two matrices
#'
#' @param mat1 First matrix
#' @param mat2 Second matrix
#' @param mat1_name Name of the first matrix, for error log (default: "first matrix")
#' @param mat2_name Name of the second matrix, for error log (default: "second matrix")
#' @param transformed Whether the matrices have been transformed, just for error log
#'
#' @return Output from \code{\link[vegan]{procrustes}} between the two matrices
#'
#' @family distances
#' @export
procMats <- function(mat1, mat2,
                     mat1_name = "first matrix", mat2_name = "second matrix",
                     transformed = TRUE){
  tokenlist <- row.names(mat1)[row.names(mat1) %in% row.names(mat2)]
  mat1b <- mat1[tokenlist, tokenlist]
  mat2b <- mat2[tokenlist, tokenlist]
  result <- tryCatch(vegan::procrustes(mat1b, mat2b, symmetric = T)$ss, error = return)
  if (inherits(result, "error")) {
    message <- paste(mat1_name, mat2_name, transformed, Sys.time())
    write(message, "procrustes.log", append = T)
    result <- vegan::procrustes(mat2b, mat1b, symmetric = T)$ss
  }
  return(result)
}

#' Mantel statistic between two matrices
#'
#' @param mat1 First matrix
#' @param mat2 Second matrix
#'
#' @return Statistic from \code{\link[vegan]{mantel}} between the two matrices.
#'
#' @family distances
#' @export
mantelMats <- function(mat1, mat2){
  # mantel equivalent of procMats
  tokenlist <- row.names(mat1)[row.names(mat1) %in% row.names(mat2)]
  res <- vegan::mantel(
    mat1[tokenlist, tokenlist],
    mat2[tokenlist, tokenlist]
    )$statistic
  return(res)
}

#' Mean euclidean distance
#'
#' Computes the euclidean distance between the individual vectors
#' and takes their mean.
#'
#' @param mat1 First matrix
#' @param mat2 Second matrix
#'
#' @return Distance between matrix
#'
#' @family distances
#' @export
eucliMats <- function(mat1, mat2){
  tokenlist <- row.names(mat1)[row.names(mat1) %in% row.names(mat2)]
  (mat1[tokenlist, tokenlist]-mat2[tokenlist, tokenlist])^2 %>%
    rowSums() %>%
    sqrt() %>%
    mean()
}

#' F-score
#'
#' Compute F-score from precision and recall.
#'
#' @param precision Precision
#' @param recall Recall
#' @param b Weight
#'
#' @return F-score
#' @export
fscore <- function(precision, recall, b = 1) {
  (1 + b^2) * (precision * recall) / (b^2 * precision + recall)
}

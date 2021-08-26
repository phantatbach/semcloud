# semcloud: Post-processing of token-level clouds.
# Copyright (C) 2021 Mariana Montes
#
# See full notice in README.md file.

#' Open token distance matrix
#'
#' Open a distance matrix stored as a \code{.pac} compressed file
#'   with a \code{.npy} file with distances and a \code{.meta} file
#'   with dimension names. Might only work for token-to-token distance
#'   matrices.
#'
#' @param input_file Name of the file
#'
#' @return Distance matrix as matrix with ids as rownames and column names.
#' @export
tokensFromPac <- function(input_file){
  temp <- utils::unzip(input_file, unzip="internal")
  tokvecs <- RcppCNPy::npyLoad(temp[2])
  metadata <- rjson::fromJSON(file = temp[1])
  dimid2item <- metadata$'row_items'
  dimnames(tokvecs) <- list(dimid2item, dimid2item)
  file.remove(temp[1], temp[2])
  return(tokvecs)
}

#' Open type level distance matrices
#' Open a distance matrix stored as a \code{.csv} file. Replaces
#' \code{\link{tokensFromPac}} for type level matrices.
#'
#' @inheritParams tokensFromPac
#'
#' @return Distance matrix as matrix with ids as rownames and column names.
#' @export
focdistsFromCsv <- function(input_file){
  focdists <- suppressWarnings(
    readr::read_tsv(input_file, col_types = readr::cols())
  ) %>%
    data.frame(row.names = "X1") %>%
    as.matrix()
  dimnames(focdists)[2] <- dimnames(focdists)[1]
  return(focdists)
}


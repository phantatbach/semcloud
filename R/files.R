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
#' @param input_directory Directory where distance matrices
#'   (as output from python module) are stored
#' @param filename Not full name of a distance matrix stored
#'   in \code{input_directory}.
#'
#' @return Distance matrix as matrix with ids as rownames and column names.
#' @export
tokens_from_pac <- function(input_directory, filename){
  input_file <- file.path(input_directory, filename)
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
#' \code{\link{tokens_from_pac}} for type level matrices.
#'
#' @inheritParams tokens_from_pac
#'
#' @return Distance matrix as matrix with ids as rownames and column names.
#' @export
focdists_from_csv <- function(input_directory, filename){
  input_file <- file.path(input_directory, filename)
  focdists <- suppressWarnings(
    readr::read_tsv(input_file, col_types = readr::cols())
  ) %>%
    data.frame(row.names = "X1") %>%
    as.matrix()
  dimnames(focdists)[2] <- dimnames(focdists)[1]
  return(focdists)
}


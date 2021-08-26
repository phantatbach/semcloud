# semcloud: Post-processing of token-level clouds.
# Copyright (C) 2021 Mariana Montes
#
# See full notice in README.md file.

#' Extract HDBSCAN info
#'
#' Run \code{\link[dbscan]{hdbscan}} on a distance matrix and gather some information.
#'
#' @param dstmtx Distance matrix
#' @param minPts Minimum points for \code{\link[dbscan]{hdbscan}}
#' @param includePlot Whether too include the plot (requires \code{cowplot}.)
#'
#' @return List: the \code{df} element is a tibble with information per token:
#'     \itemize{
#'         \item{**_id**: }{comes from the rownames of \code{dstmtx}}
#'         \item{**clusters**: }{gives the clustering of the elements}
#'         \item{**membprob**: }{indicates membership probabilities}
#'         \item{**eps**: }{returns the epsilon value}
#'     }
#'     If \code{includePlot} is \code{TRUE}, a \code{grob} of the plot is included under \code{plot}.
#' @export
#'
#' @importFrom rlang .data
extractHDBSCAN <- function(dstmtx, minPts = 8, includePlot = FALSE) {
  h <- dbscan::hdbscan(x = stats::as.dist(dstmtx), minPts = minPts)
  clusters <- stats::setNames(as.integer(h$cluster), row.names(dstmtx))
  membprob <- stats::setNames(h$membership_prob, row.names(dstmtx))
  eps <- purrr::map(
    attr(h, "hdbscan"), function(clus) {
      return(stats::setNames(clus$eps, row.names(dstmtx)[clus$contains]))
    }) %>% purrr::flatten_dbl()
  res <- list(
    df = tibble::tibble(
      `_id` = row.names(dstmtx),
      clusters = factor(clusters[.data$`_id`]),
      membprob = membprob[.data$`_id`],
      eps = eps[.data$`_id`]
    )
  )
  if (includePlot) {
    if (requireNamespace("cowplot", quietly = TRUE)) {
      res$hplot <- cowplot::as_grob(~plot(h, show_flat = TRUE))
    } else {
      warning("The plot could not be stored because `cowplot` is not installed.")
    }

  }

  return(res)
}

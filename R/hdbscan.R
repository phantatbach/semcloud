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
#' @return List: the \code{df} element is a [tibble][tibble::tibble-package] with information per token:
#'     \itemize{
#'         \item{**_id**: }{comes from the rownames of \code{dstmtx}}
#'         \item{**cluster**: }{gives the clustering of the elements}
#'         \item{**membprob**: }{indicates membership probabilities}
#'         \item{**eps**: }{returns the epsilon value}
#'     }
#'     If \code{includePlot} is \code{TRUE}, a \code{grob} of the plot is included under \code{hplot}.
#' @export
#'
#' @importFrom rlang .data
extractHDBSCAN <- function(dstmtx, minPts = 8, includePlot = FALSE) {
  requireHere('dbscan')
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
      cluster = factor(clusters[.data$`_id`]),
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

#' Map context words and HDBSCAN clusters
#'
#' The function expects a dataframe where at least you have token-id's (e.g. `_id`),
#' a column with character vectors of context words (e.g. `cws`)
#' and a column with names of clusters (e.g. `cluster`).
#' The example below shows how to also turn `;`-separated values into character vectors
#' within a [tibble][tibble::tibble-package] dataframe.
#'
#' @param variables Dataframe with IDs, clusters and lists of context words
#' @param cws_column Character string: Name of the column with the character vectors (one per row) of context words
#' @param cluster_column Character string: Name of the column with the name of the clusters (must be a factor)
#' @param b Weight for computing \code{\link{fscore}}
#'
#' @return a [tibble][tibble::tibble-package] with one row per context word per cluster, with frequency information.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' variables <- dplyr::mutate(variables, cws = stringr::str_split(cws, ";"))
#' cwsForClusters(variables, "cws", "cluster")
#'
#' }
#'
#' @importFrom rlang .data
cwsForClusters <- function(variables, cws_column, cluster_column, b=1){

  all_cws <- variables %>% dplyr::pull(cws_column) %>% purrr::flatten_chr() %>% unique()
  clusters <- variables %>% dplyr::pull(cluster_column) %>% unique() %>% sort()

  clus_cws_matrix <- matrix(nrow = length(clusters), ncol = length(all_cws),
                            dimnames = list(clusters, all_cws))
  clus_cws_matrix[is.na(clus_cws_matrix)] <- 0

  for (i in 1:nrow(variables)) {
    dat <- dplyr::slice(variables, i)
    clus <- dat[[cluster_column]]
    cws <- dat[[cws_column]][[1]]
    for (cw in unique(cws)) {
      attempt <- try(
        {clus_cws_matrix[clus,  cw] <- clus_cws_matrix[clus,  cw] + 1},
        silent = T)
      if (inherits(attempt, 'try-error')) {
        message(sprintf("cw is '%s' and cluster is '%s'", cw, clus))
        }
    }
  }
  purrr::map_df(as.character(clusters), function(clus) {
    clus_size <- nrow(dplyr::filter(variables, !!sym(cluster_column) == clus))
    this_row <- clus_cws_matrix[clus,]
    present <- names(this_row[this_row > 0])

    tibble::tibble(
      cw = present, TP = this_row[.data$cw],
      recall = .data$TP/clus_size,
      precision = purrr::map_dbl(.data$cw, ~this_row[[.x]]/sum(clus_cws_matrix[,.x])),
      Fscore = fscore(.data$precision, .data$recall, b),
      cluster = clus
    ) %>%
      dplyr::arrange(.data$cluster, dplyr::desc(.data$Fscore))
  })
}

#' Summarize HDBSCAN data for a model
#'
#' @param lemma Name of the lemma, for filenames
#' @param modelname Name of the model, for coordinates and filename
#' @param input_dir Directory where the distance matrix is stored
#' @param output_dir Directory where coordinates are stored. This directory must contain:
#'     \itemize{
#'         \item A file with the coordinates of the tokens, with a name combining \code{lemma} and \code{coords_name} and ending in \code{.tsv}.
#'         \item A file with coordinates for the context words, with a name combining \code{lemma} and \code{coords_name} and ending in \code{.cws.tsv}.
#'         \item A file with semicolon-separated lists of context words, with a name combining \code{lemma} and \code{.variables.tsv}
#'     }
#' @param coords_name The code in the coordinate files indicating the type of dimensionality reduction performed, for filenames
#' @inheritParams extractHDBSCAN
#'
#' @return list with at least two items:
#'    \itemize{
#'        \item{**coords**: }{a [tibble][tibble::tibble-package] with one row per token, the coordinates in the pertinent file, and information
#'        from \code{\link{extractHDBSCAN}}} as well as the \code{variables} file.
#'        \item{**cws**: }{a [tibble][tibble::tibble-package] with one row per context word and cluster, output
#'        from \code{\link{cwsForClusters}}}, combined with coordinates from the relevant file.
#'        \item{**hplot**: }{If \code{includePlot}, the HDBSCAN plot.}
#'    }
#' @export
#'
#' @importFrom rlang .data
summarizeHDBSCAN <- function(lemma, modelname, input_dir, output_dir, minPts = 8,
                             includePlot = FALSE, coords_name = ".tsne.30") {
  ttmx <- paste0(modelname, ".ttmx.dist.pac")
  coords_file <- file.path(output_dir, paste0(lemma, coords_name, ".tsv"))
  cw_coords_file <- file.path(output_dir, paste0(lemma, coords_name, ".cws.tsv"))
  variables_file <- file.path(output_dir, paste0(lemma, ".variables.tsv"))

  dstmtx <- tokensFromPac(file.path(input_dir, ttmx)) %>%
    transformMats(TRUE)
  variables_original <- readr::read_tsv(variables_file, show_col_types = FALSE,
                                        quote = F)
  variables <- variables_original %>%
    dplyr::select(.data$`_id`, # keep id column
           cws = stringr::str_replace(modelname, "(.+).LENGTH.*", "_cws.\\1"), # context words column
           !dplyr::starts_with("_") # columns without prefix, e.g. 'sense', if they exist
    ) %>%
    dplyr::mutate(cws = stringr::str_split(.data$cws, ";"))

  ctxt_name <- stringr::str_replace(modelname, "(.+).LENGTH.*", "_ctxt.\\1")
  if (ctxt_name %in% colnames(variables_original)) variables[['ctxt']] <- variables_original[[ctxt_name]]
  tokens <- intersect(row.names(dstmtx), variables$`_id`)

  coords <- readr::read_tsv(coords_file, show_col_types = FALSE) %>%
    dplyr::select(.data$`_id`, model.x = paste0(modelname, ".x"), model.y = paste0(modelname, ".y")) %>%
    dplyr::filter(.data$`_id` %in% tokens)
  dstmtx <- dstmtx[tokens, tokens]
  # Compute HDBSCAN and extract clusters and rest
  h <- extractHDBSCAN(dstmtx, minPts, includePlot = includePlot)

  # Merge into one dataframe
  coords <- coords %>%
    dplyr::left_join(variables, by = '_id') %>%
    dplyr::left_join(h$df, by = '_id')
  cws_per_cluster <- cwsForClusters(coords, "cws", "cluster")

  if (file.exists(cw_coords_file)) {
    cw_coords <- readr::read_tsv(file.path(cw_coords_file), show_col_types = FALSE)
    if (paste0(modelname, ".x") %in% colnames(cw_coords)) {
      cw_coords <- cw_coords %>%
        dplyr::select(cw = .data$`_id`, x = paste0(modelname, ".x"), y = paste0(modelname, ".y"))
      cws_per_cluster <- dplyr::left_join(
        cws_per_cluster,
        cw_coords,
        by = "cw"
      )
    } else {
      cws_per_cluster <- cws_per_cluster %>%
        dplyr::mutate(x = 0, y = 0)
    }
  }
  res <- list(coords = coords, cws = cws_per_cluster)
  if (includePlot) { res$hplot <- h$hplot }

  return(res)
}


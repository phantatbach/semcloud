# semcloud: Post-processing of token-level clouds.
# Copyright (C) 2021 Mariana Montes
#
# See full notice in README.md file.

#' Compute distances per cluster
#'
#' @param clustering Named vector with token IDs as names and (HDBSCAN) clusters
#'    as values. We assume that each cluster has at least 8 items.
#' @param dists Long format table with one row per pair of tokens (that are not
#'    the same), the distance between them, the cluster that the first token belongs
#'    to and whether they both belong to the same cluster.
#'
#' @return A tibble with one row per cluster and various distance-derived values:
#'    \describe{
#'       \item{min_, mean_ and max_identicals}{Minimum, mean and maximum number
#'       of identical tokens per token in the cluster.}
#'       \item{min_, mean_ and max_k8}{Minimum, mean and maximum distance from
#'       each token in the cluster and its 8th nearest neighbour.}
#'       \item{min_, mean_, max_ and sd_inner_dist}{Minimum, mean, and maximum
#'       distance, as well as their standard deviation, between each token of the
#'       cluster and all other tokens in the same cluster.}
#'       \item{min_, mean_, max_ and sd_outer_dist}{Minimum, mean, and maximum
#'       distance, as well as their standard deviation, between each token of the
#'       cluster and all other tokens in other clusters.}
#'    }
#' @export
#' @importFrom rlang .data
clusterDistance <- function(clustering, dists){
  identicals <- tibble::enframe(clustering, 'token_A', 'cluster_A') %>%
    dplyr::mutate(
      identicals = purrr::map_dbl(
        .data$token_A,
        ~nrow(dplyr::filter(dists, .data$token_A == .x, .data$same_cluster, .data$distance < 0.000001))),
      max_k8 = purrr::map_dbl(
        .data$token_A,
        ~sort(dplyr::filter(dists, .data$token_A == .x, .data$same_cluster)$distance)[8]
      )
    ) %>%
    dplyr::group_by(.data$cluster_A) %>%
    dplyr::summarize(
      'min_identicals' = min(.data$identicals),
      'mean_identicals' = mean(.data$identicals),
      'max_identicals' = max(.data$identicals),
      'min_k8' = min(.data$max_k8),
      'mean_k8' = mean(.data$max_k8),
      'max_k8' = max(.data$max_k8)
    )
  cluster_distances <- dists %>%
    dplyr::group_by(.data$cluster_A, .data$same_cluster) %>%
    dplyr::summarize(
      'min' = min(.data$distance), 'mean' = mean(.data$distance), 'max' = max(.data$distance), 'sd' = stats::sd(.data$distance),
      .groups = 'keep') %>%
    dplyr::mutate(same_cluster = dplyr::if_else(.data$same_cluster, 'inner_dist', 'outer_dist')) %>%
    tidyr::pivot_wider(names_from = 'same_cluster', values_from = c('min', 'mean', 'max', 'sd'))
  dplyr::full_join(identicals, cluster_distances, by = 'cluster_A') %>%
    dplyr::rename('cluster' = .data$cluster_A)
}

#' Summarize HDBSCAN data per cluster
#'
#' @param m Tibble with one token per row and HDBSCAN information. The \code{coords}
#'    element of a model resulting from \code{\link{summarizeHDBSCAN}}.
#'
#' @return Tibble with one row per cluster and various HDBSCAN-derived values:
#'    \describe{
#'       \item{min_, mean_, max_ and sd_cws}{Minimum, mean and maximum, as well
#'       as standard deviation, of the number of first-order context words per
#'       token in that cluster.}
#'       \item{min_, mean_, max_ and sd_eps}{Minimum, mean and maximum, as well
#'       as standard deviation, of the \eqn{\epsilon} value of the tokens in
#'       that cluster.}
#'       \item{size, rel_size}{Absolute number of tokens in the cluster and proportion
#'       of \emph{modelled} tokens covered by the cluster.}
#'       \item{deeper_than_noise}{Proportion of tokens in that cluster with an
#'       \eqn{\epsilon} value lower than the minimum \eqn{\epsilon} of noise
#'       tokens in that model.}
#'       \item{cw_tokens, _types, _ttratio}{Union of first-order context words
#'       of tokens in that cluster: number of types and of tokens and type-token
#'       ratio.}
#'    }
#'
#' @export
#'
#' @importFrom rlang .data
clusterHDBSCAN <- function(m) {
  min_eps <- m %>%
    dplyr::filter(!is.na(.data$eps)) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::summarize(eps = min(.data$eps)) %>%
    tibble::deframe()

  m %>%
    dplyr::mutate(deeper_than_noise = if ('0' %in% names(min_eps)) .data$eps < min_eps[['0']] else TRUE) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::mutate(
      cws_n = purrr::map_dbl(.data$cws, length),
      cws = purrr::map_chr(.data$cws, paste, collapse = ';')
    ) %>%
    dplyr::summarize(
      min_cws = min(.data$cws_n), mean_cws = mean(.data$cws_n),
      max_cws = max(.data$cws_n), sd_cws = stats::sd(.data$cws_n),
      min_eps = min(.data$eps), max_eps = max(.data$eps),
      mean_eps = mean(.data$eps), sd_eps = stats::sd(.data$eps),
      size = dplyr::n(), deeper_than_noise = sum(.data$deeper_than_noise)/.data$size,
      cwlist = paste(.data$cws, collapse = ';')
    ) %>%
    dplyr::mutate(
      cwlist = purrr::map(.data$cwlist, stringr::str_split, ';') %>% purrr::map(1),
      cw_tokens = purrr::map_dbl(.data$cwlist, length),
      cwlist = purrr::map(.data$cwlist, unique),
      cw_types = purrr::map_dbl(.data$cwlist, length),
      cw_ttratio = .data$cw_types/.data$cw_tokens,
      rel_size = .data$size/sum(.data$size)
    ) %>%
    dplyr::select(-'cwlist')
}

#' Compute Semvar values on clusters
#'
#' @inheritParams clusterHDBSCAN
#'
#' @return A tibble with one row per cluster and output from
#'   \code{\link{separationkNN}} and \code{\link[cluster]{silhouette}}
#'   for each class, based on the \emph{coordinates} in the input and both including
#'   and excluding noise tokens.
#'
#' @export
#'
#' @importFrom rlang .data
clusterSeparation <- function(m) {
  if (!requireNamespace('cluster', quietly = TRUE)) {
    stop("Package `cluster` or `semvar` needed.")
  }

  sil_func <- function(dists, classes){
    sils <- summary(cluster::silhouette(as.numeric(classes), dists))$clus.avg.widths
    clus_names <- names(sils)
    dim(sils) <- NULL
    names(sils) <- clus_names
    sils
  }
  knn_func <- function(dists, classes){
    separationkNN(dists, classes, k = 8)$classqual
  }

  classes <- as.character(m$cluster)
  if (length(unique(classes)) == 1) return(tibble::tibble(cluster = unique(classes)))

  if (min(table(classes)) == 1) { # cover for the occasional single noise point :)
    removed <- names(table(classes)[table(classes) == 1])
    m <- m %>% dplyr::filter(.data$cluster != removed)
    classes <- m$cluster
  }

  full_dists <- stats::dist(m[c("model.x", "model.y")], diag = TRUE, upper = TRUE)
  no_noise <- stats::dist(
    dplyr::filter(m, .data$cluster != '0')[c("model.x", "model.y")],
    diag = TRUE, upper = TRUE
    )

  tibble::tibble(
    cluster = unique(classes),
    kNN_full = knn_func(full_dists, classes)[.data$cluster],
    SIL_full = sil_func(full_dists, classes)[.data$cluster],
    kNN_no_noise = knn_func(no_noise, classes[classes != '0'])[.data$cluster],
    SIL_no_noise = sil_func(no_noise, classes[classes != '0'])[.data$cluster]
  )
}

#' Classify clouds in a model
#'
#' @param mdata Data of a model, as given by \code{\link{summarizeHDBSCAN}}.
#' @param mname Name of a model; the name of \code{mdata}.
#' @param ttmx_dir Directory where the token-level distance matrices are stored.
#' @param suffix Suffix of the file names for the token-level distance matrices;
#'   the function assumes that the name of the file is the name of the medoid
#'   plus the suffix.
#'
#' @return A table with one row per cluster in the model, the columns created by
#'    \code{\link{clusterSeparation}}, \code{\link{clusterHDBSCAN}} and
#'    \code{\link{clusterDistance}} and the classification of each cluster based
#'    on the Nephological Shapes from \insertCite{montes_2021;textual}{semcloud}
#'    (see \href{https://cloudspotting.marianamontes.me/shapes.html}{Chapter 5}
#'    for a full description and examples).
#'
#' @references
#'   \insertAllCited{}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' purrr::imap_dfr(models$medoidCoords, classifyModel, ttmx_dir = 'path/to/dir')
#' }
#' @importFrom Rdpack reprompt
classifyModel <- function(mdata, mname, ttmx_dir, suffix = '.ttmx.dist.pac'){
  m <- mdata$coords

  clustering <- tibble::deframe(dplyr::select(m, '_id', 'cluster'))

  dists <- tokensFromPac(file.path(ttmx_dir, paste0(mname, suffix))) %>%
    tibble::as_tibble(rownames = 'token_A') %>%
    tidyr::pivot_longer(-'token_A', names_to = 'token_B', values_to = 'distance') %>%
    dplyr::filter(.data$token_A != .data$token_B) %>%
    dplyr::mutate(cluster_A = clustering[.data$token_A],
                  same_cluster = clustering[.data$token_B] == .data$cluster_A)

  model_data <- dplyr::full_join(
    clusterSeparation(m), clusterHDBSCAN(m),
    by = 'cluster'
  ) %>%
    dplyr::full_join(clusterDistance(clustering, dists), by = 'cluster')

  if (nrow(model_data) == 1) {
    model_class <- dplyr::mutate(
      model_data,
      cloud_type = "Cirrostratus",
      Hail = .data$max_identicals >= 8
    )
    clouds <- 'Cirrostratus'
  } else {
    model_class <- dplyr::mutate(
      model_data,
      cloud_type = dplyr::case_when(
        .data$cluster == '0' ~ "Cirrostratus",
        .data$rel_size >= 0.5 ~ "Cumulonimbus",
        .data$kNN_full >= 0.75 & .data$SIL_full >= 0.5 &
          .data$mean_inner_dist <= 0.5 &
          (.data$deeper_than_noise >= 0.9 | sum(.data$rel_size) >= 0.9) ~ "Cumulus",
        (sum(.data$rel_size) <= 0.25 | .data$SIL_full >= 0.5 | .data$mean_inner_dist <= 0.2) &
          (.data$deeper_than_noise > 0.5 | sum(.data$rel_size) > 0.9) ~ "Stratocumulus",
        TRUE ~ "Cirrus"
      ),
      Hail = .data$max_identicals >= 8
    )
    clouds <- setdiff(sort(unique(model_class$cloud_type)), 'Cirrostratus')
  }

  if (sum(model_class$Hail) > 1) clouds <- c(clouds, 'Hail')
  clouds <- paste(clouds, collapse = '-')

  model_class %>%
    dplyr::mutate(model = mname, maincat = clouds) %>%
    dplyr::select('model', 'maincat', 'cluster', 'cloud_type', 'Hail', dplyr::everything()) %>%
    dplyr::arrange(.data$cluster)
}

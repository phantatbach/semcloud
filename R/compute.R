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
#' @param d distance matrix, which will be forced into a \code{matrix} object
#'   and then a \code{dist} object
#' @param dim number of dimensions, by default 2. As argument \code{k} for
#'   \code{\link[vegan]{metaMDS}} and \code{dims} for \code{\link[Rtsne]{Rtsne}}
#' @param technique either "mds" to run \code{\link[vegan]{metaMDS}} or
#'   "tsne" to run \code{\link[Rtsne]{Rtsne}}
#' @param perp perplexity value for \code{\link[Rtsne]{Rtsne}}, default is 30.
#'   This value is ignored when \code{\link[vegan]{metaMDS}} is run.
#' @param seed seed to keep randomness at check
#'
#' @return output of either \code{\link[vegan]{metaMDS}} or \code{\link[Rtsne]{Rtsne}}.
#' @export
getFit <- function(d, dim = 2, technique, perp = 30, seed = 8541){
  set.seed(seed)
  dst <- stats::as.dist(as.matrix(d))
  if (technique == 'mds') {
    fit <- vegan::metaMDS(dst, k=dim, trymax=20, trace=FALSE)
  } else if (technique == "tsne") {
    fit <- Rtsne::Rtsne(dst, dims=dim, perplexity=perp,
                        theta=0.0, check.duplicates = FALSE,
                        max_iter = 1000, is_distance = TRUE)
  } else {
    stop("`technique` must be 'mds' or 'tsne'.")
  }
  return(fit)
}

#' Extract coordinates from fit
#'
#' @param fit Dimensionality reduction result; output from \code{\link{getFit}}
#' @param modelname Name of the model, for the names of the columns
#' @param rownames List of names for the rows
#' @param d Output from previous run or, if it's the first run, empty string.
#' @param source Technique used in dimensionality reduction: either "mds"
#'    for output from \code{\link[vegan]{metaMDS}} or "tsne" for \code{\link[Rtsne]{Rtsne}}
#'    (it matches `technique` from \code{\link{getFit}}).
#'
#' @return Tibble with the coordinates of each element.
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
  } else {
    stop("`source` must be 'mds' or 'tsne'.")
  }
  df <- tibble::tibble(
    `_id` = rownames,
    model.x = x,
    model.y = y) %>%
    dplyr::rename_all(list(~stringr::str_replace(.data, "model", modelname)))
  if (!is.character(d)) {
    df  <- df %>%
      dplyr::full_join(dplyr::select(d, !dplyr::starts_with(modelname)), by="_id") %>%
      replace(.data, is.na(.data), 0.0)
  }
  return(df)
}

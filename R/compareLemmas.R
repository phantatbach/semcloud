# semcloud: Post-processing of token-level clouds.
# Copyright (C) 2021 Mariana Montes
#
# See full notice in README.md file.


#' Compute distances between models
#'
#' Three functions to compute distances between (token) matrices are supported:
#'    - "euclidean" runs \code{\link{eucliMats}}, which computes tokenwise euclidean distances and averages over them;
#'    - "procrustes" runs \code{\link[vegan]{procrustes}} via \code{\link{procMats}};
#'    - "mantel" runs \code{\link[vegan]{mantel}} via \code{\link{mantelMats}}.
#'
#' @param mnames List of model names.
#' @param input_dir Directory where the models are stored.
#' @param transformed Whether the distance matrices
#' @param fun Function to calculate the distances.
#'
#' @return A distance matrix (\code{matrix} object) with models as rows and columns.
#' @export
customDist <- function(mnames, input_dir, transformed = TRUE,
                        fun = c("euclidean", "procrustes", "mantel")) {
  res <- matrix(
    nrow = length(mnames),
    ncol = length(mnames),
    dimnames = list(mnames, mnames)) # create a matrix with the models as rows and columns
  names(mnames) <- mnames
  tokvecs <- purrr::map(mnames, function(m) {
    mat <- tokensFromPac(file.path(input_dir, paste0(m, "ttmx.dist.pac", sep = ".")))
    if (transformed == TRUE) return(transformMats(mat, asDist = fun != "euclidean")) else return(mat)
  })

  pb <- utils::txtProgressBar(min = 0, max = length(mnames), style = 3)
  for (row in mnames) {
    mat1 <- tokvecs[[row]]
    for (col in mnames) {
      if (row == col) {
        res[row, col] <- 0
      } else {
        mat2 <- tokvecs[[col]]
        if (!is.na(res[col, row])) {
          res[row, col] <- res[col, row]
        } else {
          if (fun == "procrustes") {
            res[row, col] <- procMats(mat1, mat2, row, col, transformed = transformed)
          } else if (fun == "mantel") {
            res[row, col] <- 1-mantelMats(mat1, mat2) #for distance instead of similarity
          } else if (fun == "euclidean") {
            res[row, col] <- eucliMats(mat1, mat2) #for distance instead of similarity
          } else {
            stop("Function non valid, mut be either 'euclidean', 'procrustes' or 'mantel'")
          }
        }
      }
    }
    utils::setTxtProgressBar(pb, which(mnames == row))
  }
  close(pb)

  print('Matrix created.')
  return(res)
}

#' Compare models of a lemma
#'
#' Compute pairwise distances between the lemmas, store the distance matrix,
#' reduce to two dimensions with \code{\link[vegan]{metaMDS}}, store and return summary..
#'
#' @param lemma Name of the lemma, for the filenames
#' @param output_dir Directory where the model information is and will be stored
#' @inheritParams customDist
#'
#' @return a [tibble][tibble::tibble-package] with minimal information
#' @export
#'
#' @importFrom rlang .data
compLemma <- function(lemma, input_dir, output_dir, transformed = TRUE,
                      fun = "euclidean"){

  # load data on the models
  models <- readr::read_tsv(file.path(output_dir, paste0(lemma, '.models.tsv')), col_types = readr::cols())
  if ('model.x' %in% colnames(models)) { models <- dplyr::select(models, -.data$model.x, -.data$model.y) }

  models.names <- models$`_model`

  distfile <- file.path(output_dir, paste0(lemma, ".models.dist.tsv"))
  if (file.exists(distfile)) {
    distmtx <- readr::read_tsv(distfile, col_types = readr::cols()) %>%
      matricizeCloud() %>% stats::as.dist()
  } else {
    distances <- customDist(models.names, input_dir, transformed = transformed, fun = fun)
    suffix <- if (fun == "euclidean") "" else fun

    distmtx <- getDistMtx(distances, output_dir, lemma, suffix)
  }

  dst.MDS <- vegan::metaMDS(distmtx, k = 2, trymax = 30)
  stress <- dst.MDS$stress
  print(stress)

  #### Add coordinates to models' files
  models.w.coords <- dst.MDS$points %>% as.data.frame %>%
    tibble::rownames_to_column %>% tibble::as_tibble %>%
    stats::setNames(.data, c('_model', 'model.x', 'model.y')) %>%
    dplyr::left_join(models, by = '_model')

  # Save new models file

  suffix <- if (fun == "euclidean") ".models" else paste0(".models.", fun)
  models.w.coords %>% readr::write_tsv(file.path(output_dir, paste0(lemma, suffix, '.tsv')))
  print("Models saved")

  reg <- tibble::tibble(
    type = lemma,
    models = nrow(models),
    stress = round(stress, 4),
    date = Sys.Date()
  )
  return(reg)

}

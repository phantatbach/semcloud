# semcloud: Post-processing of token-level clouds.
# Copyright (C) 2021 Mariana Montes
#
# See full notice in README.md file.

#' Get Clouds from distance matrix
#'
#' Read distance matrices from different models,
#' run dimensional reduction for visualization based on different techniques
#' and store the coordinates corresponding to each model in a dataframe per technique.
#' The names of the models will be found in the `models` file and their paths will be
#' searched for in `input_dir`: if a file is not found, a warning will be issued.
#'
#' @param input_dir Directory where the token distance matrices are stored.
#' @param output_dir Directory where the data will be stored.
#' @param files_list Liste of filenames within `input_dir`.
#' @param lemma Name of the lemma, for filenames
#' @param solutions Named list of techniques to run for visualization possible `technique` values in \code{\link{getFit}}.
#' @param logrank Whether to transform the matrices with \code{\link{transformMats}}.
#' @param type Whether to open the files with \code{\link{tokensFromPac}} (for "token") or \code{\link{focdistsFromCsv}} (otherwise).
#'
#' @return List of stresses (emtpy if "mds" is not given.)
#' @export
getClouds <- function(input_dir, output_dir, files_list, lemma, solutions, logrank = TRUE, type = "token"){

  d <- purrr::map(solutions, function(solution){ # set up main file
    suffix <- if (type == "token") ".tsv" else ".cws.tsv"

    filename.full <- file.path(output_dir, paste0(lemma, solution, suffix))
    df <- if (file.exists(filename.full)) readr::read_tsv(filename.full, show_col_types = FALSE, lazy = FALSE) else "" # in case we have stored it
    return(list(filename = filename.full, df = df))
  })

  stresses <- list()

  pb <- utils::txtProgressBar(min = 0, max = length(files_list), style = 3)
  for (file in files_list){
    fname <- file.path(input_dir, file)
    if (!file.exists(fname)) {
      warning(paste0(fname, " does not exist; model skipped."))
      next
    }
    # for each of the models
    modelname <- textreuse::filenames(textreuse::filenames(textreuse::filenames(file)))
    # obtain the distance matrix
    dists <- if (type == "token") tokensFromPac(fname) else focdistsFromCsv(fname)

    if (logrank) { dists <- transformMats(dists, TRUE) } # log-transform the matrix

    for (sol in names(solutions)) { # for each of the "solutions"

      technique <- stringr::str_extract(solutions[[sol]], "[a-z]+")
      perp <- readr::parse_integer(stringr::str_extract(solutions[[sol]], "[0-9]+"))

      # dimensionality reduction
      if (technique == "mds" | nrow(dists)/perp > 3) {
        # run algorithm
        fit <- getFit(dists, dim = 2, technique = technique, perp = perp)

        # extract coordinates from algorithm
        d[[sol]]$df <- getCoords(fit, modelname, row.names(dists), d[[sol]]$df, technique)

        if (technique == 'mds') {
          stresses[modelname] <- fit$stress
        }
      } else {
        print(paste0(nrow(dists), " - ", perp))
      }

    }
    utils::setTxtProgressBar(pb, which(files_list == file))
  }
  close(pb)

  for (sol in d) {
    readr::write_tsv(sol$df, sol$filename) #store the coordinates for each solution
  }

  return(stresses)
}

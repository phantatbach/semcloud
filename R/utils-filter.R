#' Extract window filter
#'
#' @inheritParams filterFoc
#'
#' @return Integer vector of length 2 (or list) with left and right window span lengths.
windowFilter <- function(foc_param) {
  windows <- foc_param %>%
    stringr::str_extract('\\d+-\\d+') %>%
    stringr::str_split('-')
  readr::parse_integer(windows[[1]])
}

#' Extract part-of-speech filter
#'
#' @inheritParams filterFoc
#'
#' @return Character vector: empty if no \code{pos} filter should be applied,
#'   values to be selected otherwise.
posFilter <- function(foc_param) {
  if (stringr::str_ends(foc_param, 'lex')) {
    c('noun', 'adj', 'adv', 'verb')
  } else {
    c()
  }
}


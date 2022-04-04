#' Correct Distance values
#'
#' Remove rows that should not count in the distance computation and recalculate distances.
#'
#' @param cws Dataframe as outputted by \code{listContextwords}
#' (\href{https://github.com/montesmariana/semasioFlow}{semasioFlow} Python module).
#' @param distance_corrector_func Function to be applied to the \code{word} variable to remove irrelevant lines.
#'
#' @return Tibble with the context word ids and the corrected distance as \code{conc_distance}.
#'
#' @importFrom rlang .data
correctDistance <- function(cws, distance_corrector_func){
  cws %>%
    dplyr::group_by(.data$token_id, .data$side) %>%
    dplyr::arrange(.data$distance) %>%
    dplyr::filter(distance_corrector_func(.data$word)) %>%
    dplyr::mutate(conc_distance = seq(dplyr::n())) %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$cw_id, .data$conc_distance)
}

#' Reorder PPMI columns
#'
#' The values of the \code{ppmi} will be used for weight selection and can be
#' included as superindices.
#'
#'
#' @param ppmi Dataframe with a row per context word and columns for frequency
#'    and association measures/
#' @param pmi_columnname Name of the column with the weighting values of interest.
#'
#' @return Dataframe with a \code{cw} column for the context word and a
#'    \code{my_weight} column for the weighting.
#'
#' @importFrom rlang .data
reorderPpmi <- function(ppmi, pmi_columnname){
  if (pmi_columnname %in% colnames(ppmi)) {
    ppmi %>%
      dplyr::select(.data$cw, my_weight = !!rlang::sym(pmi_columnname)) %>%
      dplyr::filter(!is.na(.data$my_weight))
  } else if (length(colnames(dplyr::select(ppmi, dplyr::starts_with(pmi_columnname)))) > 1) {
    ppmi %>% dplyr::select(.data$cw, dplyr::starts_with(pmi_columnname)) %>%
      dplyr::rename_all(stringr::str_remove, paste0(pmi_columnname, '_')) %>%
      tidyr::pivot_longer(-.data$cw, names_to = 'target_lemma', values_to = 'my_weight') %>%
      dplyr::filter(!is.na(.data$my_weight))
  } else {
    stop("You need to indicate a column present in the ppmi file.")
  }
}

#' Set up Concordancer
#'
#' Prepare dataframes for \code{\link{getContext}}.
#'
#' @param lemma Name of the lemma: for default filenames
#' @param input_dir Directory where the files are stored
#' @param cws_detail_path Path to a dataframe with one row per context word per token
#'   and context words with information from the token. Created by
#'   \code{listContextWords} in the \href{https://github.com/montesmariana/semasioFlow}{semasioFlow} Python module.
#' @param ppmi_path Path to a dataframe with one context word per row and frequency information
#' @param pmi_columnname Name (or prefix) of the column in the dataframe found in
#'   \code{ppmi_path} where weighting values are stored.
#' @param distance_corrector_func Function to filter the rows of the dataframe in
#'   \code{cws_detail_path} based on the values of the \code{word} column, to recalculate
#'   distances between words.
#' @param lemma_from_tid_fun Function to extract the target lemma from the tokenID.
#'
#' @return Enriched dataframe with one row per context word per token, weight values,
#'   corrected distances and a column indicating the right target lemma
#'   (in case you have more than one).
#' @export
setupConcordancer <- function(
  lemma = "", input_dir = "",
  cws_detail_path = file.path(input_dir, paste0(lemma, '.cws.detail.tsv')),
  ppmi_path = file.path(input_dir, paste0(lemma, '.ppmi.tsv')),
  pmi_columnname = 'pmi_4',
  distance_corrector_func = function(word) !stringr::str_starts(word, '<'),
  lemma_from_tid_fun = function(tid) paste(stringr::str_split(tid, '/')[[1]][-c(3, 4)], collapse = '/')
  ){
  files_exist <- {
    file.exists(cws_detail_path) & file.exists(ppmi_path)
  }
  if (!files_exist) {
    stop("You must give valid paths for `variables_path`, `cws_detail_path` and `ppmi_path`.")
  }
  cws <- readr::read_tsv(cws_detail_path, show_col_types = F)
  conc_distance <- correctDistance(cws, distance_corrector_func)
  cws <- cws %>% dplyr::left_join(conc_distance, by = 'cw_id')

  ppmi <- readr::read_tsv(ppmi_path, show_col_types = F) %>%
    reorderPpmi(pmi_columnname)

  if ('target_lemma' %in% colnames(ppmi)) {
    cws <- dplyr::left_join(cws, ppmi, by = c('cw', 'target_lemma'))
  } else {
    cws <- dplyr::left_join(cws, ppmi, by = 'cw')
  }
  cws
}


#' Clean Word column
#'
#' Function to convert word forms coming from the \emph{QLVLNewsCorpus} towards
#' something more friendly for HTML rendering.
#'
#' @param txt word form
#'
#' @return Converted string.
cleanWord <- function(txt) { # cleaning needed for QLVLNewsCorpus
  stringr::str_replace(txt, "and(.+);", "&\\1;") %>%
    stringr::str_replace("`", "'") %>% stringr::str_replace('</sentence>', '<br>')
}

#' Filter by First Order Parameters
#'
#' @param foc_param Character string coding the relevant first-order parameters.
#' @param tid_data Subsection of a context-word-by-token dataframe, as outputted
#'   by \code{\link{setupConcordancer}}, with information for one token.
#' @param cw_selection Vector of context words selected by the model for that token.
#' @param is_dep_fun Function that takes \code{foc_param} as input and returns
#'   \code{TRUE} if dependency information should be collected and \code{FALSE}
#'   if the model is based on bag-of-words instead.
#' @param max_steps_fun Function that takes \code{foc_param} as input and returns,
#'   for dependency-based models, the maximum number of steps in the dependency path
#'   to accept as viable context words.
#' @param window_filter_fun Function that takes \code{foc_param} as input and returns
#'   a vector or list with two elements: the left and right window sizes (for
#'   bag-of-words models).
#' @param pos_filter_fun Function that takes \code{foc_param} as input and returns
#'   a vector. If the vector is empty, no \code{pos} filter is implemented, while if
#'   it has values, the rows with \code{pos} included in that vector will be selected.
#' @param bound_filter_fun Function that takes \code{foc_param} as input and returns
#'   \code{TRUE} if words \strong{outside} the sentence are modelled and
#'   \code{FALSE} if they are not.
#'
#' @return Enriched dataframe including columns with filtering information.
#' @export
#'
#' @importFrom rlang .data
filterFoc <- function(
  foc_param, tid_data, cw_selection,
  is_dep_fun = function(foc_param) stringr::str_starts(foc_param, 'LEMMA'),
  max_steps_fun = function(foc_param) if (foc_param == 'LEMMAPATH2') 2 else 3,
  window_filter_fun = windowFilter, pos_filter_fun = posFilter,
  bound_filter_fun = function(foc_param) stringr::str_starts(foc_param, 'nobound')
){

  if (is_dep_fun(foc_param)) {
    # this is not super reliable, but reliable enough
    max_steps <- max_steps_fun(foc_param)
    tid_data <- tid_data %>%
      dplyr::mutate(
        focsel = .data$same_sentence & .data$steps <= max_steps & .data$cw %in% cw_selection
        )
  } else {
    window_sel <- window_filter_fun(foc_param)
    # modify pos_sel if you have other parameters :)
    pos_sel <- pos_filter_fun(foc_param)
    bound_sel <- bound_filter_fun(foc_param)
    tid_data <- tid_data %>%
      dplyr::mutate(
        by_bound = bound_sel | .data$same_sentence,
        by_win = !is.na(.data$conc_distance) & dplyr::if_else(
          .data$side == 'L',
          .data$conc_distance <= window_sel[[1]],
          .data$conc_distance <= window_sel[[2]]),
        by_pos = length(pos_sel) == 0 | .data$pos %in% pos_sel,
        focsel = .data$by_bound & .data$by_win & .data$by_pos & .data$cw %in% cw_selection
        )
  }
  tid_data
}

#' Filter by Weighting Information
#'
#' @param weight_param Character string coding the relevant weighting information.
#' @param tid_data Subsection of a context-word-by-token dataframe enriched by
#'   \code{\link{filterFoc}}.
#' @param weight_filter_fun Function that takes \code{weight_param} as input and
#'   returns \code{TRUE} if weighting should be included and \code{FALSE} if it
#'   should be ignored.
#' @param weight_as_sup Should weight values be included in superscript?
#' @param threshold Weighting threshold for a context word to be included
#'
#' @return Enriched dataframe including columns with filtering information.
#' @export
#'
#' @importFrom rlang .data
filterWeight <- function(
  weight_param, tid_data,
  weight_filter_fun = function(weightparam) stringr::str_ends(weightparam, 'no', negate = TRUE),
  weight_as_sup = FALSE,
  threshold = 0
  ) {
  tid_data <- tid_data %>%
    dplyr::mutate(
      weighted = weight_filter_fun(weight_param) | (!is.na(.data$my_weight) & .data$my_weight > threshold),
      bolden = .data$focsel & .data$weighted,
      sup = if (weight_as_sup) sprintf('<sup>%s</sup>', round(.data$my_weight, 3)) else '',
      word = dplyr::if_else(.data$bolden, sprintf('<strong>%s</strong>%s', .data$word, .data$sup), .data$word)
    )
  tid_data
}

#' Context from concordance dataframe
#'
#' From a dataframe with one row per context word of a token and columns for
#' elements (in particular, \code{word} with the appropriate enrichment),
#' create a concordance line.
#'
#' @param tid Token ID.
#' @param cws Dataframe with one row per context word per token, as outputted
#'   by \code{\link{setupConcordancer}}.
#' @param foc_param Character string with information on first order filters.
#' @param weight_param Character string with information on weighting filters.
#' @param cw_selection List or character string of semicolon-separated values
#'   with the context words captured by the model.
#' @param clean_word_fun Function to clean the \code{word} column.
#' @param to_remove Vector of \code{word} values to remove altogether.
#' @param ... Arguments to be passed to \code{\link{filterFoc}} and \code{\link{filterWeight}}.
#'
#' @return Character string with (weighted) concordance line.
#' @export
#'
#' @examples
#' \dontrun{
#' cws <- setupConcordancer(lemma, input_dir)
#' getContext(tokenID, cws, 'bound5-5lex', 'PPMIweight', c('a/det', 'horse/noun', 'fast/adj'))
#' }
#'
#' @importFrom rlang .data
getContext <- function(
  tid, cws, foc_param = NA, weight_param = NA,
  cw_selection = NA,
  clean_word_fun = cleanWord,
  to_remove = c('<sentence>'),
  ...) {
  extra_args <- list(...)

  tid_data <- cws %>%
    dplyr::filter(.data$token_id == tid) %>%
    dplyr::mutate(word = purrr::map_chr(.data$word, clean_word_fun)) %>%
    dplyr::filter(!.data$word %in% to_remove)
  target <- tid_data[tid_data$side == 'target', 'word']

  if (!is.na(foc_param) & !is.na(cw_selection)) {
    if (length(cw_selection) == 1) {
      cw_selection <- stringr::str_split(cw_selection, ';')[[1]]
    }
    tid_data <- R.utils::doCall(
      filterFoc,
      foc_param = foc_param, tid_data = tid_data, cw_selection = cw_selection,
      args = extra_args)
    if (!is.na(weight_param)) {
      tid_data <- R.utils::doCall(
        filterWeight,
        weight_param = weight_param, tid_data = tid_data,
        args = extra_args
        )
      }
  }

  left <- tid_data %>% dplyr::filter(.data$side == 'L') %>% dplyr::arrange(dplyr::desc(.data$distance))
  right <- tid_data %>% dplyr::filter(.data$side == 'R') %>% dplyr::arrange(.data$distance)
  stringr::str_glue('{paste(left$word, collapse = " ")} <span class="target">{target}</span> {paste(right$word, collapse = " ")}')
  # stringr::str_glue("{paste(left$word, collapse = ' ')} <span class='target'>{target}</span> {paste(right$word, collapse = ' ')}")

}

#' Tailor weighting function to a lemma
#'
#' @param cws Dataframe with one row per token per context word, output
#'   of \code{\link{setupConcordancer}}.
#' @param foc_param_fun Function that takes the name of a model and returns a
#'   character string with first-order filters (to be used as \code{foc_param}
#'   in \code{\link{getContext}} and thus \code{\link{filterFoc}}).
#' @param weight_param_fun Function that takes the name of a model and returns a
#'   character string with weighting filters (to be used as \code{weight_param}
#'   in \code{\link{getContext}} and thus \code{\link{filterWeight}}).
#' @param sup_weight_fun Function that takes the output of \code{weight_param_fun}
#'   as input and returns \code{TRUE} if weighting values should be added as
#'   superindices and \code{FALSE} if they should not. To be fed to \code{\link{filterWeight}}.
#' @param ... Arguments to be passed to \code{\link{getContext}}.
#'
#' @return Function that takes a token ID, a list of context words and a model name
#'   as input and calls \code{\link{getContext}}.
weightLemma <- function(
  cws,
  foc_param_fun = function(m) stringr::str_split(m, '\\.')[[1]][[1]],
  weight_param_fun = function(m) stringr::str_split(m, '\\.')[[1]][[2]],
  sup_weight_fun = function(weightparam) stringr::str_ends(weightparam, 'weight'),
  ...
  ){
  extra_args <- list(...)
  function(tid, cw_selection, m){
    extra_args$weight_as_sup <- sup_weight_fun(weight_param_fun(m))
    R.utils::doCall(
      getContext,
      tid = tid, cws = cws, cw_selection = cw_selection,
      foc_param = foc_param_fun(m),
      weight_param = weight_param_fun(m),
      args = extra_args,
      .ignoreUnusedArgs = FALSE
      )
  }
}

#' Create weighted concordance lines
#'
#' This code assumes that the lists of context words are stored in columns
#' starting with \code{_cws}, followed by the name of the lemma and a model name
#' from which first-order parameter settings can be extracted (all separated
#' by periods). Then it creates
#' the weighted concordance lines and stores them in columns following the same
#' name pattern, but starting with \code{_ctxt}. In addition, it creates a non
#' weighted column \code{_ctxt.raw}.
#'
#' @param variables Dataframe with one row per token ID and at least columns
#'   prefixed by \code{_cws.} with semicolon-separated context words.
#' @param cws Dataframe with one row per token ID per context word as outputted
#'   by \code{\link{setupConcordancer}}.
#' @param lemma Name of the lemma, to process model names
#' @param ... Arguments to be passed to \code{\link{weightLemma}} and
#'   \code{\link{getContext}}, in order to adapt to different ways of coding
#'   parameter settings. See \code{vignette('weightConcordance')}.
#'
#' @return Enriched \code{variables} dataframe with columns containing weighted
#'   concordance lines.
#' @export
#'
#' @examples
#' \dontrun{
#' cws <- setupConcordancer(lemma, input_dir)
#' variables <- readr::read_tsv('path/to/file', lazy = F)
#' ctxts <- weightConcordance(variables, cws, lemma)
#' }
#'
#' @importFrom rlang .data
weightConcordance <- function(variables, cws, lemma, ...) {
  extra_args <- list(...)
  models <- variables %>% dplyr::select(dplyr::starts_with('_cws.')) %>%
    dplyr::rename_all(stringr::str_remove, paste0('_cws.', lemma, '.')) %>%
    colnames()

  # add raw context
  rawModel <- function(tid) {R.utils::doCall(
    getContext,
    tid = tid, cws = cws, args = extra_args, .ignoreUnusedArgs = FALSE
  )}
  variables <- variables %>%
    dplyr::mutate(`_ctxt.raw` = purrr::map_chr(.data$`_id`, rawModel))


  weightModel <- R.utils::doCall(
    weightLemma,
    cws = cws,
    args = extra_args,
    .ignoreUnusedArgs = FALSE)

  for (m in models) {
    if (!paste('_ctxt', lemma, m, sep = '.') %in% colnames(variables)) {
      variables <- variables %>%
        dplyr::mutate(!!paste('_ctxt', lemma, m, sep='.') := purrr::map2_chr(
          .data$`_id`, !!rlang::sym(paste('_cws', lemma, m, sep='.')),
          weightModel, m = m
        )
        )
    }

  }

  variables
}

variables <- readr::read_tsv('blik.variables.tsv', show_col_types = F, lazy = F)
cws <- setupConcordancer(cws_detail_path = 'blik.cws.detail.tsv',
                         ppmi_path = 'blik.ppmi.tsv')

test_that("file is downloaded when paths are given", {
  expect_equal(nrow(cws), 160)
  expect_equal(length(colnames(cws)), 19)
  expect_true('my_weight' %in% colnames(cws))
})

test_that("file is downloaded when lemma and path are given", {
  cws2 <- setupConcordancer(lemma = "blik", input_dir = './')
  expect_equal(nrow(cws2), 160)
  expect_equal(length(colnames(cws2)), 19)
  expect_true('my_weight' %in% colnames(cws2))

})

ctxts <- weightConcordance(variables, cws, 'blik')

test_that("raw text is correct", {
  expect_match(ctxts$`_ctxt.raw`[[1]], "grote gebaar , maar dus niet van allebei , hè meneer . <br> Werp één <span class=\"target\">blik</span> in de vitrine van de servieswinkel en je beseft het : hier staat een Fransman")
})

test_that("weighted text is correct", {
  expect_match(ctxts$`_ctxt.blik.LEMMAREL2.PPMIweight`[[1]], "grote gebaar , maar dus niet van allebei , hè meneer . <br> <strong>Werp</strong><sup>5.85</sup> <strong>één</strong><sup>1.233</sup> <span class=\"target\">blik</span> in de vitrine van de servieswinkel en je beseft het : hier staat een Fransman")
})

test_that("unweighted text is correct", {
  expect_match(ctxts$`_ctxt.blik.bound5-5lex.PPMIselection`[[1]], "grote gebaar , maar dus niet van allebei , hè meneer . <br> <strong>Werp</strong> één <span class=\"target\">blik</span> in de <strong>vitrine</strong> van de servieswinkel en je beseft het : hier staat een Fransman")
})

nosup <- weightConcordance(variables, cws, 'blik', sup_weight_fun = function(weightparam) return(FALSE))

test_that("superindices can be overriden", {
  expect_match(nosup$`_ctxt.blik.LEMMAREL2.PPMIweight`[[1]], "grote gebaar , maar dus niet van allebei , hè meneer . <br> <strong>Werp</strong> <strong>één</strong> <span class=\"target\">blik</span> in de vitrine van de servieswinkel en je beseft het : hier staat een Fransman")
})

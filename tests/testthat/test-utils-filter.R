test_that("window is properly filtered", {
  both10 <- windowFilter("10-10")
  one5one20 <- windowFilter("whatever5-20kmqdfq")
  expect_equal(both10[[1]], 10)
  expect_equal(both10[[2]], 10)
  expect_equal(one5one20[[1]], 5)
  expect_equal(one5one20[[2]], 20)
})

test_that("pos filter works", {
  expect_equal(length(posFilter("10all")), 0)
  expect_equal(length(posFilter("10lex")), 4)
  expect_equal(length(posFilter("lex10")), 0)
})

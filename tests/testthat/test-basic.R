
test_that("DRBayes package loads correctly", {
  expect_true(require(DRBayes, quietly = TRUE))
})

test_that("Basic functions are available", {
  expect_true(exists("DRBayes.PC"))
  expect_true(exists("B.LM"))
  expect_true(exists("B.Logit"))
})

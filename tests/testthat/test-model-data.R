# Regression tests for the design-matrix defects fixed in the posterior coupling
# estimator. Each of these fails against the pre-fix code.

sim_data <- function(n = 120, seed = 42) {
  set.seed(seed)
  d <- data.frame(X1 = rnorm(n), X2 = rnorm(n), X3 = rnorm(n))
  d$A <- rbinom(n, 1, plogis(0.3 + 0.6 * d$X1 - 0.4 * d$X2))
  d$Y <- 1 + 2 * d$A + 0.8 * d$X1 + 0.5 * d$X2 + rnorm(n)
  d
}

test_that("row names are treated as labels, not positions", {
  d <- sim_data()

  # Row names that are a permutation of 1:n: indexing by position would silently
  # pick different rows while keeping the row count right, so nothing downstream
  # would notice.
  permuted <- d[sample(nrow(d)), ]
  reset    <- permuted
  rownames(reset) <- NULL

  a <- prepare_model_data(Y ~ A + X1 + X2, A ~ X1 + X2, permuted)
  b <- prepare_model_data(Y ~ A + X1 + X2, A ~ X1 + X2, reset)

  expect_equal(unname(a$Y), unname(b$Y))
  expect_equal(unname(a$A), unname(b$A))
  expect_equal(unname(a$Z.lm), unname(b$Z.lm))

  # Row names running past nrow() after a subset must not index out of range.
  subset_rows <- d[d$X3 > 0, ]
  expect_gt(max(as.numeric(rownames(subset_rows))), nrow(subset_rows))
  expect_silent(prepare_model_data(Y ~ A + X1 + X2, A ~ X1 + X2, subset_rows))
})

test_that("a formula without an intercept is rejected rather than corrupted", {
  d <- sim_data()
  expect_error(
    prepare_model_data(Y ~ A + X1 + X2 - 1, A ~ X1 + X2, d),
    "must include an intercept"
  )
  expect_error(
    prepare_model_data(Y ~ A + X1 + X2, A ~ X1 + X2 + 0, d),
    "must include an intercept"
  )
})

test_that("counterfactual designs recompute every treatment-derived column", {
  d <- sim_data()
  md <- prepare_model_data(Y ~ A + X1 + A:X1, A ~ X1, d)

  expect_identical(colnames(md$Z.lm1), colnames(md$Z.lm))
  expect_identical(colnames(md$Z.lm0), colnames(md$Z.lm))

  # Under A = 1 the interaction column equals X1; under A = 0 it is zero. The
  # pre-fix code copied the observed treatment into these columns, which
  # collapsed m1(X) - m0(X) to the coefficient on A for every unit.
  expect_equal(unname(md$Z.lm1[, "A"]), rep(1, nrow(d)))
  expect_equal(unname(md$Z.lm0[, "A"]), rep(0, nrow(d)))
  expect_equal(unname(md$Z.lm1[, "A:X1"]), unname(d$X1))
  expect_equal(unname(md$Z.lm0[, "A:X1"]), rep(0, nrow(d)))
})

test_that("data-dependent bases keep the fitted parameterisation", {
  d <- sim_data()
  md <- prepare_model_data(Y ~ A + poly(X3, 2), A ~ X1, d)

  # Only the treatment differs between the fitted and counterfactual frames, so
  # the orthogonal polynomial basis must be untouched.
  expect_equal(md$Z.lm[, "poly(X3, 2)1"], md$Z.lm1[, "poly(X3, 2)1"])
  expect_equal(md$Z.lm[, "poly(X3, 2)2"], md$Z.lm0[, "poly(X3, 2)2"])
})

test_that("the treatment variable accepts the usual binary encodings", {
  d <- sim_data()

  numeric_md <- prepare_model_data(Y ~ A + X1, A ~ X1, d)
  expect_identical(numeric_md$A, as.integer(d$A))

  d_logical <- d
  d_logical$A <- as.logical(d$A)
  expect_identical(prepare_model_data(Y ~ A + X1, A ~ X1, d_logical)$A,
                   as.integer(d$A))

  d_factor <- d
  d_factor$A <- factor(d$A, levels = c(0, 1))
  expect_message(
    factor_md <- prepare_model_data(Y ~ A + X1, A ~ X1, d_factor),
    "Treating '1' as treated"
  )
  expect_identical(factor_md$A, as.integer(d$A))

  d_bad <- d
  d_bad$A <- d$A + 1
  expect_error(prepare_model_data(Y ~ A + X1, A ~ X1, d_bad), "must be binary")

  d_three <- d
  d_three$A <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  expect_error(prepare_model_data(Y ~ A + X1, A ~ X1, d_three),
               "exactly two levels")
})

test_that("only rows complete in both models are used", {
  d <- sim_data()
  d$X1[1:5]  <- NA          # enters both models
  d$X3[6:10] <- NA          # enters neither of the formulas below

  md <- prepare_model_data(Y ~ A + X1 + X2, A ~ X1 + X2, d)
  expect_equal(md$data_info$n_observations, nrow(d) - 5)
  expect_equal(md$data_info$missing_observations, 5)
})

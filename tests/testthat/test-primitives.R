# The primitives the compiled moment conditions are built from, checked against
# the R functions they exist to reproduce. Each is reachable from R only so that
# its agreement with R is a test rather than a claim.

skip_without_compiled <- function() {
  skip_if_not(exists("drb_mean_one_pass", envir = asNamespace("DRBayes")),
              "no compiled kernel in this build")
}

# The width R sums in, which the kernels are told rather than deduce. Passing
# the same answer the package passes is what makes these comparisons against
# colMeans, mean and rowMeans meaningful; passing the other one would compare
# two different reductions and fail on any build where the two differ.
wide <- DRBayes:::long_double_sums()

adversarial <- list(
  ordinary      = c(1.5, -2.25, 0.125, 7),
  single        = 3.5,
  empty         = numeric(0),
  with_na       = c(0.5, NA_real_, 2),
  with_nan      = c(0.5, NaN, 2),
  with_inf      = c(0.5, Inf, 2),
  with_neg_inf  = c(0.5, -Inf, 2),
  both_inf      = c(-Inf, Inf),
  all_zero      = rep(0, 8),
  large_spread  = c(1e300, -1e300, 1)
)


test_that("the one-pass mean is the mean colMeans computes", {
  skip_without_compiled()
  for (nm in names(adversarial)) {
    x <- adversarial[[nm]]
    if (length(x) == 0L) next
    expect_identical(DRBayes:::drb_mean_one_pass(x, wide),
                     .colMeans(x, length(x), 1L), info = nm)
  }

  skip_on_cran()
  set.seed(1)
  for (i in 1:200) {
    x <- stats::rnorm(3000, 5, 3)
    expect_identical(DRBayes:::drb_mean_one_pass(x, wide), .colMeans(x, 3000L, 1L))
  }
})


test_that("the two-pass mean is the mean mean() computes", {
  skip_without_compiled()
  for (nm in names(adversarial)) {
    x <- adversarial[[nm]]
    if (length(x) == 0L) next
    expect_identical(DRBayes:::drb_mean_two_pass(x, wide), mean(x), info = nm)
  }

  skip_on_cran()
  set.seed(2)
  for (i in 1:200) {
    x <- stats::rnorm(3000, 5, 3)
    expect_identical(DRBayes:::drb_mean_two_pass(x, wide), mean(x))
  }
})


test_that("each compiled reduction lands where the R one it stands for lands", {
  skip_on_cran()
  skip_without_compiled()

  # colMeans sums once and divides; mean() adds a correction pass. The package
  # uses colMeans for the parametric moment conditions and mean() for the
  # nonparametric one, so each compiled reduction has to follow its own.
  set.seed(3)
  for (i in 1:200) {
    x <- stats::rnorm(3000, 5, 3)
    expect_identical(DRBayes:::drb_mean_one_pass(x, wide), .colMeans(x, 3000L, 1L))
    expect_identical(DRBayes:::drb_mean_two_pass(x, wide), mean(x))
  }
})


test_that("the two reductions are not interchangeable where R sums in double", {
  skip_on_cran()
  skip_without_compiled()
  # Whether R's two reductions differ is a property of the build, not of the
  # package: accumulating in double they part company on 185 of these 200
  # vectors, and accumulating in long double the correction rounds to nothing
  # and they agree on all of them. This is the case that makes keeping two
  # reductions necessary; the test above is the one that holds everywhere.
  skip_if(wide, "this build accumulates in long double, where the two agree")

  set.seed(3)
  differ <- 0L
  for (i in 1:200) {
    x <- stats::rnorm(3000, 5, 3)
    if (!identical(mean(x), .colMeans(x, 3000L, 1L))) differ <- differ + 1L
  }
  expect_gt(differ, 100L)
})


test_that("blocked row means equal unblocked row means exactly", {
  skip_without_compiled()
  set.seed(4)
  x <- matrix(stats::rnorm(60 * 37), 60, 37)

  # The accumulator carries across blocks rather than averaging block means, so
  # the partition cannot change the answer. This is what lets the propensity
  # score means be computed a block at a time instead of materialising the whole
  # draws-by-observations matrix.
  reference <- .rowMeans(x, 60L, 37L)
  for (block in c(1L, 2L, 7L, 36L, 37L, 100L)) {
    expect_identical(DRBayes:::drb_row_means(x, block, wide), reference,
                     info = paste("block", block))
  }
})


test_that("the range tracker carries missingness rather than stepping over it", {
  skip_without_compiled()

  # Every comparison against NaN is false, so a running minimum written the
  # obvious way skips it silently. R's min() and max() return NA instead, and
  # the positivity guard depends on that: a propensity score that cannot be
  # evaluated has to refuse, not pass.
  for (nm in names(adversarial)) {
    x <- adversarial[[nm]]
    got <- DRBayes:::drb_range(x)
    if (length(x) == 0L) {
      expect_identical(got, numeric(0), info = nm)
    } else {
      expect_identical(got, c(min(x), max(x)), info = nm)
    }
  }
})


test_that("the inverse links are the ones stats uses", {
  skip_without_compiled()
  eta <- c(-800, -709.8, -100, -38.5, -8.3, -1, 0, 1, 8.3, 36.7, 36.8, 800,
           NA_real_, NaN, Inf, -Inf)

  expect_identical(DRBayes:::drb_inverse_link(eta, 0L), eta)
  expect_identical(DRBayes:::drb_inverse_link(eta, 1L), stats::plogis(eta))
  expect_identical(DRBayes:::drb_inverse_link(eta, 2L), stats::pnorm(eta))
})


test_that("the fused weight is the divided weight, on a logistic link", {
  skip_without_compiled()
  eta <- c(-20, -5, -1, 0, 1, 5, 20)

  # Only where the cancellation in 1 - e is harmless. Past that the two forms
  # genuinely disagree, and the fused one is the correct one: at eta = 20 with
  # A = 0 the exact weight is -(1 + exp(20)) = -485165196.40979, which the fused
  # form returns and the divided form misses by eight significant digits,
  # because 1 - plogis(20) has already lost them.
  moderate <- abs(eta) <= 10
  for (a_value in c(0, 1)) {
    a  <- rep(a_value, length(eta))
    ps <- stats::plogis(eta)
    expect_equal(DRBayes:::drb_ipw_weight(eta, a, 1L)[moderate],
                 ((a - ps) / (ps * (1 - ps)))[moderate],
                 tolerance = 1e-12, info = paste("A =", a_value))
  }

  far <- 20
  expect_equal(DRBayes:::drb_ipw_weight(far, 0, 1L), -(1 + exp(far)),
               tolerance = 1e-15)

  # Under probit there is no such collapse, so the kernel must take the divided
  # form. Asking for the fused one here would be wrong by about 15% at eta = 1,
  # with no error and no NaN to show for it.
  a  <- c(1, 0, 1, 0, 1, 0, 1)
  ps <- stats::pnorm(eta)
  expect_identical(DRBayes:::drb_ipw_weight(eta, a, 2L),
                   (a - ps) / (ps * (1 - ps)))
})

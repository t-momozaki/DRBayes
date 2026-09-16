# The diagnostics of Vehtari et al. (2021). The Figure 2 tests are the ones that
# matter: they are the scenarios the traditional statistic is known to miss.

ar1_chains <- function(n, chains, rho = 0.3, sd = 1, shift = 0, seed = 1) {
  set.seed(seed)
  out <- matrix(0, n, chains)
  for (m in seq_len(chains)) {
    x <- numeric(n)
    x[1] <- rnorm(1, sd = sd / sqrt(1 - rho^2))
    for (t in 2:n) x[t] <- rho * x[t - 1] + rnorm(1, sd = sd)
    out[, m] <- x
  }
  out[, chains] <- out[, chains] + shift
  out
}

test_that("rank normalisation and folding catch what plain R-hat misses", {
  skip_on_cran()

  # Scenario 1 of Figure 2: one of four chains has a third of the variance. The
  # chains agree on location, so a statistic that only compares locations sees
  # nothing; folding turns the scale difference into a location difference.
  x <- ar1_chains(1000, 4, seed = 11)
  x[, 4] <- x[, 4] / 3
  expect_gt(rhat(x), 1.01)
  expect_lt(rhat_basic(split_chains(x)), 1.01)

  # Scenario 2: a Cauchy target with one chain shifted. The variance does not
  # exist, so the moment-based statistic estimates nothing; ranks are always
  # well defined.
  num <- ar1_chains(1000, 4, seed = 21)
  den <- ar1_chains(1000, 4, seed = 22)
  cauchy <- num / den
  cauchy[, 4] <- cauchy[, 4] + 2
  expect_gt(rhat(cauchy), 1.01)
  expect_lt(rhat_basic(split_chains(cauchy)), 1.01)
})

test_that("well mixed chains are not flagged", {
  skip_on_cran()
  expect_lt(rhat(ar1_chains(1000, 4, seed = 31)), 1.01)

  num <- ar1_chains(1000, 4, seed = 41)
  den <- ar1_chains(1000, 4, seed = 42)
  expect_lt(rhat(num / den), 1.01)
})

test_that("the statistics agree with the posterior package", {
  skip_on_cran()
  skip_if_not_installed("posterior")

  x <- ar1_chains(1000, 4, rho = 0.9, seed = 51)
  expect_equal(rhat(x), posterior::rhat(x), tolerance = 1e-8)

  # The effective sample sizes differ by a factor N/(N-1) per lag, because
  # equation (3.10) is implemented literally here while Stan and posterior
  # average the biased autocovariance and set rho_0 to one to compensate. The
  # gap grows with autocorrelation, and rho = 0.9 is close to the worst case.
  expect_equal(ess_bulk(x), posterior::ess_bulk(x), tolerance = 0.05)
  expect_equal(ess_tail(x), posterior::ess_tail(x), tolerance = 0.05)
  expect_equal(mcse_mean(x), posterior::mcse_mean(x), tolerance = 0.05)

  # Through the public entry point the agreement is much tighter, because the
  # array is sliced per parameter exactly as posterior expects.
  arr <- array(x, c(dim(x), 1L), dimnames = list(NULL, NULL, "theta"))
  d <- convergence_diagnostics(arr)
  expect_equal(d$rhat, posterior::rhat(x), tolerance = 1e-8)
  expect_equal(d$ess_bulk, posterior::ess_bulk(x), tolerance = 0.01)
})

test_that("convergence_diagnostics summarises every parameter", {
  skip_on_cran()
  x <- array(rnorm(500 * 4 * 3), c(500, 4, 3),
             dimnames = list(NULL, NULL, c("a", "b", "c")))
  d <- convergence_diagnostics(x)

  expect_identical(d$parameter, c("a", "b", "c"))
  expect_named(d, c("parameter", "mean", "sd", "rhat", "ess_bulk", "ess_tail",
                    "mcse_mean"))
  expect_true(all(d$rhat < 1.01))
  expect_true(all(d$ess_bulk > 400))
})

test_that("drbayes_pc warns when either posterior has not converged", {
  skip_on_cran()
  set.seed(7)
  n <- 300
  d <- data.frame(X1 = rnorm(n))
  d$A <- rbinom(n, 1, plogis(0.3 + 0.7 * d$X1))
  d$Y <- 1 + 2 * d$A + 0.8 * d$X1 + rnorm(n)

  run <- function(...) {
    utils::capture.output(
      res <- suppressMessages(drbayes_pc(Y ~ A + X1, A ~ X1, d,
                                         family = "gaussian",
                                         outcome.model = bayes_lm,
                                         ps.model = bayes_logit, ...)),
      file = nullfile())
    res
  }

  set.seed(99)
  expect_warning(fit <- run(mc = 60, bn = 10, thin = 1),
                 "do not meet the convergence criteria")
  # The message must name the model, the statistic and the parameters, and say
  # that a small R-hat cannot be trusted when the effective sample size is small.
  expect_warning(run(mc = 60, bn = 10, thin = 1), "propensity score model")
  expect_warning(run(mc = 60, bn = 10, thin = 1), "R-hat itself unreliable")

  expect_s3_class(fit$diagnostics, "data.frame")
  expect_true(all(c("model", "parameter", "rhat", "ess_bulk", "ess_tail") %in%
                    names(fit$diagnostics)))

  set.seed(99)
  expect_error(run(mc = 60, bn = 10, thin = 1, diagnostics = "error"),
               "do not meet the convergence criteria")

  set.seed(99)
  quiet <- run(mc = 60, bn = 10, thin = 1, diagnostics = "none")
  expect_null(quiet$diagnostics)
})

test_that("rank plot counts are uniform over the pooled ranks", {
  set.seed(2)
  x <- array(rnorm(400 * 4), c(400, 4, 1), dimnames = list(NULL, NULL, "theta"))
  rp <- rank_plot(x, "theta", bins = 10, plot = FALSE)

  # Every bin holds the same number of pooled ranks by construction, so the
  # rows must sum to S / bins; departure from uniformity shows up per chain.
  expect_identical(dim(rp$counts), c(10L, 4L))
  expect_true(all(rowSums(rp$counts) == 400 * 4 / 10))
  expect_true(all(colSums(rp$counts) == 400))
})

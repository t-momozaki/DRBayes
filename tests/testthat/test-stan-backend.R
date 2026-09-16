# The Stan backend has to fit the same model as the Gibbs sampler it mirrors,
# so every comparison below asks whether the two posterior means differ by more
# than Monte Carlo error can explain. The criterion is the one used elsewhere in
# this package: |mean_a - mean_b| / sqrt(mcse_a^2 + mcse_b^2) < 3.
#
# These tests need CmdStan, which the CRAN machines do not have, so they are
# skipped twice over.

skip_without_cmdstan <- function() {
  skip_on_cran()
  skip_if_not_installed("instantiate")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  skip_if_not(instantiate::stan_cmdstan_exists(),
              "CmdStan is not installed on this machine")
}

drop_burn <- function(draws, burn) {
  draws[seq(burn + 1L, dim(draws)[1]), , , drop = FALSE]
}

# One z score per parameter, comparing two draws arrays of the same model.
mcse_z <- function(a, b) {
  parameters <- dimnames(a)[[3]]
  vapply(parameters, function(p) {
    draws_a <- a[, , p, drop = FALSE]
    draws_b <- b[, , p, drop = FALSE]
    abs(mean(draws_a) - mean(draws_b)) /
      sqrt(posterior::mcse_mean(draws_a)^2 + posterior::mcse_mean(draws_b)^2)
  }, numeric(1))
}

continuous_data <- function(n = 300, seed = 11) {
  set.seed(seed)
  A  <- rbinom(n, 1, plogis(0.4 * rnorm(n)))
  X  <- cbind(A = A, X1 = rnorm(n), X2 = rnorm(n))
  list(X = X,
       Y  = 1 + 2 * A - 0.5 * X[, "X1"] + 0.3 * X[, "X2"] + rnorm(n),
       Yl = rbinom(n, 1, plogis(-0.3 + 0.9 * A + 0.6 * X[, "X1"] -
                                  0.4 * X[, "X2"])),
       Yp = rbinom(n, 1, pnorm(-0.2 + 0.6 * A + 0.4 * X[, "X1"])))
}

sparse_data <- function(n = 250, p = 10, seed = 12) {
  set.seed(seed)
  V <- matrix(rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("V", seq_len(p))))
  A <- rbinom(n, 1, plogis(0.5 * V[, 1]))
  list(X = cbind(A = A, V),
       Y  = 1 + 1.5 * A + 2 * V[, 1] - 1.5 * V[, 2] + rnorm(n),
       Yl = rbinom(n, 1, plogis(-0.2 + 1.0 * A + 1.2 * V[, 1])),
       Yp = rbinom(n, 1, pnorm(-0.2 + 0.7 * A + 0.8 * V[, 1])))
}


test_that("bayes_stan returns the layout the Gibbs samplers return", {
  skip_without_cmdstan()
  d <- continuous_data(n = 200)

  draws <- bayes_stan("lm", d$Y, d$X, mc = 500, warmup = 500, chains = 4)

  expect_true(is.array(draws))
  expect_identical(dim(draws), c(500L, 4L, 4L))
  expect_identical(dimnames(draws)[[3]],
                   c("(Intercept)", "A", "X1", "X2"))
  expect_null(dimnames(draws)[[1]])
  expect_null(dimnames(draws)[[2]])
  expect_true(all(is.finite(draws)))

  # The array is what drbayes_pc expects, so the diagnostics travel alongside
  # it rather than inside it.
  nuts <- attr(draws, "nuts_diagnostics")
  expect_s3_class(nuts, "data.frame")
  expect_identical(nrow(nuts), 4L)
  expect_named(nuts, c("chain", "divergent", "max_treedepth", "ebfmi",
                       "step_size", "accept_stat"))
  expect_true(all(is.finite(nuts$ebfmi)))

  info <- attr(draws, "stan_info")
  expect_identical(info$model, "lm")
  expect_identical(info$warmup, 500L)
  expect_true(nzchar(info$cmdstan_version))
})


test_that("bayes_stan and bayes_lm agree within Monte Carlo error", {
  skip_without_cmdstan()
  d <- continuous_data()

  set.seed(101)
  gibbs <- drop_burn(bayes_lm(d$Y, d$X, mc = 8000, chains = 4), 2000)
  set.seed(101)
  stan  <- bayes_stan("lm", d$Y, d$X, mc = 2000, warmup = 1000, chains = 4)

  expect_true(all(mcse_z(gibbs, stan) < 3))
})


test_that("bayes_stan and bayes_logit agree within Monte Carlo error", {
  skip_without_cmdstan()
  d <- continuous_data()

  set.seed(102)
  gibbs <- drop_burn(bayes_logit(d$Yl, d$X, mc = 12000, chains = 4), 2000)
  set.seed(102)
  stan  <- bayes_stan("logit", d$Yl, d$X, mc = 2000, warmup = 1000, chains = 4)

  expect_true(all(mcse_z(gibbs, stan) < 3))
})


test_that("bayes_stan and bayes_probit agree within Monte Carlo error", {
  skip_without_cmdstan()
  d <- continuous_data()

  set.seed(103)
  gibbs <- drop_burn(bayes_probit(d$Yp, d$X, mc = 12000, chains = 4), 2000)
  set.seed(103)
  stan  <- bayes_stan("probit", d$Yp, d$X, mc = 2000, warmup = 1000,
                      chains = 4)

  expect_true(all(mcse_z(gibbs, stan) < 3))
})


test_that("bayes_stan and bayes_lm_hs agree within Monte Carlo error", {
  skip_without_cmdstan()
  d <- sparse_data()

  set.seed(201)
  gibbs <- drop_burn(
    bayes_lm_hs(d$Y, d$X, mc = 20000, chains = 4, unshrunk = 1L), 5000)
  set.seed(201)
  # A few divergences survive even at adapt_delta 0.99, which is what the
  # warning is for; the point of this test is that the means still line up.
  stan <- suppressWarnings(
    bayes_stan("lm_hs", d$Y, d$X, mc = 2000, warmup = 1000, chains = 4,
               unshrunk = 1L))

  expect_true(all(mcse_z(gibbs, stan) < 3))
})


test_that("bayes_stan and bayes_logit_hs agree within Monte Carlo error", {
  skip_without_cmdstan()
  d <- sparse_data()

  set.seed(202)
  gibbs <- drop_burn(
    bayes_logit_hs(d$Yl, d$X, mc = 20000, chains = 4, unshrunk = 1L), 5000)
  set.seed(202)
  stan <- suppressWarnings(
    bayes_stan("logit_hs", d$Yl, d$X, mc = 2000, warmup = 1000, chains = 4,
               unshrunk = 1L))

  expect_true(all(mcse_z(gibbs, stan) < 3))
})


test_that("bayes_stan and bayes_probit_hs agree within Monte Carlo error", {
  skip_without_cmdstan()
  d <- sparse_data()

  set.seed(203)
  gibbs <- drop_burn(
    bayes_probit_hs(d$Yp, d$X, mc = 20000, chains = 4, unshrunk = 1L), 5000)
  set.seed(203)
  stan <- suppressWarnings(
    bayes_stan("probit_hs", d$Yp, d$X, mc = 2000, warmup = 1000, chains = 4,
               unshrunk = 1L))

  expect_true(all(mcse_z(gibbs, stan) < 3))
})


test_that("unshrunk exempts exactly the coefficients it names", {
  skip_without_cmdstan()
  # A weakly identified treatment effect, so that shrinking it is visible.
  set.seed(21)
  n <- 250
  V <- matrix(rnorm(n * 15), n, 15,
              dimnames = list(NULL, paste0("V", seq_len(15))))
  A <- rbinom(n, 1, plogis(0.5 * V[, 1]))
  X <- cbind(A = A, V)
  Y <- 1 + 0.5 * A + 2 * V[, 1] - 1.5 * V[, 2] + rnorm(n, sd = 3)

  set.seed(301)
  exempt <- suppressWarnings(
    bayes_stan("lm_hs", Y, X, mc = 1500, warmup = 1000, chains = 4,
               unshrunk = 1L))
  set.seed(301)
  shrunk <- suppressWarnings(
    bayes_stan("lm_hs", Y, X, mc = 1500, warmup = 1000, chains = 4,
               unshrunk = integer(0)))

  expect_gt(abs(mean(exempt[, , "A"])), abs(mean(shrunk[, , "A"])))
  # A genuinely null covariate is shrunk whichever set is exempted.
  expect_lt(abs(mean(exempt[, , "V10"])), 0.5)
})


test_that("bayes_stan validates its arguments", {
  # Every check here runs before the model is compiled, so CmdStan is not
  # needed to reach any of these messages.
  d <- continuous_data(n = 50)

  expect_error(bayes_stan("poisson", d$Y, d$X), "must be one of")
  expect_error(bayes_stan(c("lm", "logit"), d$Y, d$X), "single character")
  expect_error(bayes_stan("lm", d$Y, d$X, warmup = 0), "positive integer")
  expect_error(bayes_stan("lm", d$Y, d$X, adapt_delta = 1),
               "strictly between 0 and 1")
  expect_error(bayes_stan("lm", d$Y, d$X, theta_prior = -1), "PRECISION")
  expect_error(bayes_stan("lm", d$Y, d$X, sigma_prior = c(1, 0)),
               "two positive numbers")
  expect_error(bayes_stan("lm_hs", d$Y, d$X, unshrunk = 1:3),
               "nothing for the horseshoe prior to shrink")
  expect_error(bayes_stan("lm_hs", d$Y, d$X, beta0_prior = 0), "PRECISION")

  # The function name is accepted as well as the bare model name, because that
  # is what a reader of the rest of the package has in mind.
  expect_identical(DRBayes:::validate_stan_model_name("bayes_probit_hs"),
                   "probit_hs")
})


test_that("a Stan file is shipped for every model bayes_stan offers", {
  for (model in DRBayes:::stan_model_names()) {
    file <- system.file("stan", paste0("bayes_", model, ".stan"),
                        package = "DRBayes")
    expect_true(nzchar(file))
  }
})


test_that("bayes_stan says how to install CmdStan when it is missing", {
  # assert_cmdstan() asks for instantiate before it asks for CmdStan, so
  # without instantiate the refusal under test is never reached and a different
  # one is returned instead.
  skip_if_not_installed("instantiate")

  # Mocked, so that the message a user without CmdStan sees is checked on a
  # machine that has it.
  local_mocked_bindings(cmdstan_is_available = function() FALSE,
                        .package = "DRBayes")

  expect_error(DRBayes:::assert_cmdstan(), "install_cmdstan")
  expect_error(DRBayes:::assert_cmdstan(), "stan-dev.r-universe.dev")
  # A user who asked for Stan must not be handed a different sampler quietly.
  expect_error(DRBayes:::assert_cmdstan(), "Nothing is substituted for Stan")
})

# The nonparametric outcome path of Section 5.4, where the outcome model is
# supplied as draws of its fitted means rather than as coefficients.

# Draws from both models, plus the fitted means those outcome draws imply. The
# outcome model omits X2, so the moment condition is genuinely violated and the
# tilting has something to do; the effect of X2 is small enough that Algorithm 1
# reaches the constraint from the draws it is given.
np_fixture <- function(n = 300, mc = 1500, seed = 11) {
  set.seed(seed)
  dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  dat$A <- rbinom(n, 1, plogis(0.3 + 0.8 * dat$X1 - 0.5 * dat$X2))
  dat$Y <- 1 + 2 * dat$A + 0.9 * dat$X1 + 0.3 * dat$X2 + rnorm(n)

  md  <- prepare_model_data(Y ~ A + X1, A ~ X1 + X2, dat)
  otc <- bayes_lm(Y = md$Y, X = md$X.lm, mc = mc, chains = 1L)
  ps  <- bayes_logit(Y = md$A, X = md$X.ps, mc = mc, chains = 1L)

  keep      <- seq(mc / 3 + 1, mc, by = 2)
  betas.otc <- matrix(otc[keep, 1L, ], ncol = dim(otc)[3])
  betas.ps  <- matrix(ps[keep, 1L, ], ncol = dim(ps)[3])

  list(data = dat, md = md, betas.otc = betas.otc, betas.ps = betas.ps,
       mu  = tcrossprod(betas.otc, md$Z.lm),
       mu1 = tcrossprod(betas.otc, md$Z.lm1),
       mu0 = tcrossprod(betas.otc, md$Z.lm0))
}

# One observation, an intercept-only propensity score pinned at one half and a
# single fitted mean per draw, so that equation (3.4) reduces to -2 mu and a
# test can hand the tilting whatever moment condition it wants to see handled.
np_moment_fixture <- function(BB) {
  list(d = tilting_data(Z.lm = matrix(1, 1L, 1L), Z.ps = matrix(1, 1L, 1L),
                        A = 1, Y = 0,
                        inverse_link = function(x) x,
                        ps_inverse_link = stats::plogis,
                        ps_formula_text = "A ~ 1",
                        outcome.mu = matrix(-BB / 2, ncol = 1L)),
       betas.ps = matrix(0, nrow = length(BB), ncol = 1L))
}

fit_quietly <- function(...) {
  suppressWarnings(drbayes_pc(..., diagnostics = "none", verbose = FALSE))
}

test_that("fitted means of a linear model reproduce the parametric path", {
  skip_on_cran()
  fx <- np_fixture()

  # The same draws, the same propensity score model and the same seed, reaching
  # the outcome mean by the two different routes. This is the only comparison
  # that can show the two paths compute one estimand rather than two.
  parametric <- fit_quietly(Y ~ A + X1, A ~ X1 + X2, fx$data,
                            family = "gaussian",
                            outcome.samples = fx$betas.otc,
                            ps.samples = fx$betas.ps,
                            bn = 0, thin = 1, method = "is", seed = 1)
  nonparametric <- fit_quietly(Y ~ A + X1, A ~ X1 + X2, fx$data,
                               outcome.mu  = fx$mu,
                               outcome.mu1 = fx$mu1,
                               outcome.mu0 = fx$mu0,
                               ps.samples = fx$betas.ps,
                               bn = 0, thin = 1, seed = 1)

  expect_true(parametric$smc$converged)
  expect_true(nonparametric$smc$converged)

  # The untilted g-formula is an arithmetic identity between the two paths.
  expect_equal(nonparametric$g.comp, parametric$g.comp, tolerance = 1e-10)

  # The tilted one is resampled, so it is held to the Monte Carlo error of the
  # posterior mean it estimates.
  mcse <- stats::sd(parametric$g.comp) / sqrt(length(parametric$g.comp))
  expect_lt(abs(mean(nonparametric$pc) - mean(parametric$pc)), 0.5 * mcse)
  expect_lt(abs(stats::sd(nonparametric$pc) /
                  stats::sd(parametric$pc) - 1), 0.05)

  # And the tilting itself is the same tilting: same parameter, same solution
  # of equation (3.9), same weights behind it.
  expect_equal(nonparametric$smc$lambda, parametric$smc$lambda,
               tolerance = 1e-6)
  expect_equal(nonparametric$smc$B_mean, parametric$smc$B_mean,
               tolerance = 1e-6)
  expect_equal(nonparametric$smc$ess, parametric$smc$ess, tolerance = 1e-6)
})

test_that("the moment condition from fitted means matches the parametric", {
  skip_on_cran()
  fx <- np_fixture()
  d  <- tilting_data(Z.lm = fx$md$Z.lm, Z.ps = fx$md$Z.ps, A = fx$md$A,
                     Y = fx$md$Y, inverse_link = function(x) x,
                     ps_inverse_link = stats::plogis,
                     ps_formula_text = "A ~ X1 + X2",
                     outcome.mu = fx$mu)

  expect_equal(moment_np(fx$betas.ps, d),
               moment_ipw(fx$betas.ps, fx$betas.otc, d), tolerance = 1e-12)

  # Equation (F.2) reads its fitted means from the same place.
  ds <- d
  ds$subclass_weight <- subclass_weights(fx$betas.ps, ds, 5L)
  expect_equal(moment_np(fx$betas.ps, ds),
               moment_subclass(fx$betas.ps, fx$betas.otc, ds),
               tolerance = 1e-12)
})

test_that("a propensity score model can be sampled alongside supplied means", {
  skip_on_cran()
  fx   <- np_fixture()
  draw <- seq(501, 1500, by = 2)

  fit <- fit_quietly(Y ~ A + X1, A ~ X1 + X2, fx$data,
                     outcome.mu  = fx$mu,
                     outcome.mu1 = fx$mu1,
                     outcome.mu0 = fx$mu0,
                     mc = 1500, bn = 500, thin = 2, chains = 1L, seed = 3)

  # The default propensity score sampler is chosen as usual, and the fit holds
  # one draw of the treatment effect per retained propensity score draw.
  expect_length(fit$pc, length(draw))
  expect_length(fit$g.comp, length(draw))
  expect_identical(fit$link, "identity")

  # Only the propensity score model was sampled here, so it is the only one
  # there are convergence diagnostics to report.
  diagnosed <- suppressWarnings(
    drbayes_pc(Y ~ A + X1, A ~ X1 + X2, fx$data,
               outcome.mu = fx$mu, outcome.mu1 = fx$mu1,
               outcome.mu0 = fx$mu0,
               mc = 1500, bn = 500, thin = 2, chains = 1L, seed = 3,
               verbose = FALSE))
  expect_identical(unique(diagnosed$diagnostics$model), "propensity score")
  expect_identical(nrow(diagnosed$diagnostics), ncol(fx$md$Z.ps))
})

test_that("the number of draws must match the propensity score side", {
  skip_on_cran()
  fx <- np_fixture()

  # Half the propensity score draws are thinned away here, so the fitted means
  # no longer line up with them and the pairing would be silently wrong.
  expect_error(
    fit_quietly(Y ~ A + X1, A ~ X1 + X2, fx$data,
                outcome.mu = fx$mu, outcome.mu1 = fx$mu1,
                outcome.mu0 = fx$mu0,
                ps.samples = fx$betas.ps, bn = 0, thin = 2),
    "one row of fitted means per retained propensity score draw")
})

test_that("the three matrices are validated against each other and the data", {
  n   <- 20L
  dat <- data.frame(X1 = seq_len(n) / n, A = rep(0:1, each = n / 2),
                    Y = seq_len(n) / n)
  mu  <- matrix(0, nrow = 8L, ncol = n)
  fit <- function(...) {
    fit_quietly(Y ~ A + X1, A ~ X1, dat, ...)
  }

  expect_error(fit(outcome.mu = mu, outcome.mu1 = mu),
               "outcome.mu0 is missing")
  expect_error(fit(outcome.mu1 = mu), "outcome.mu, outcome.mu0 are missing")
  expect_error(fit(outcome.mu = as.data.frame(mu), outcome.mu1 = mu,
                   outcome.mu0 = mu),
               "must be a numeric matrix")
  expect_error(fit(outcome.mu = mu, outcome.mu1 = mu,
                   outcome.mu0 = mu[, -1L, drop = FALSE]),
               "must have the same dimensions")
  expect_error(fit(outcome.mu = mu[, -1L, drop = FALSE],
                   outcome.mu1 = mu[, -1L, drop = FALSE],
                   outcome.mu0 = mu[, -1L, drop = FALSE]),
               "has 19 columns but the models are fitted on 20 observations")

  na_mu       <- mu
  na_mu[3, 4] <- NA_real_
  expect_error(fit(outcome.mu = mu, outcome.mu1 = na_mu, outcome.mu0 = mu),
               "outcome.mu1 holds missing or infinite fitted means")

  # An outcome model cannot be supplied twice over.
  expect_error(fit(outcome.mu = mu, outcome.mu1 = mu, outcome.mu0 = mu,
                   outcome.prior = "horseshoe"),
               "outcome.prior has nothing left to describe")
  expect_error(fit(outcome.mu = mu, outcome.mu1 = mu, outcome.mu0 = mu,
                   outcome.model = bayes_lm, link = "identity"),
               "outcome.model, link have nothing left to describe")
})

test_that("the sweep is refused for draws it cannot rejuvenate", {
  skip_on_cran()
  fx <- np_fixture()
  supplied <- list(outcome.mu = fx$mu, outcome.mu1 = fx$mu1,
                   outcome.mu0 = fx$mu0, ps.samples = fx$betas.ps,
                   bn = 0, thin = 1)

  expect_error(
    do.call(fit_quietly, c(list(Y ~ A + X1, A ~ X1 + X2, fx$data),
                           supplied, list(method = "smc"))),
    "Gaussian kernel on the outcome model's parameter vector")

  # Left unasked for, the method is Algorithm 1 rather than the default sweep,
  # and the caller is not made to say so.
  fit <- do.call(fit_quietly, c(list(Y ~ A + X1, A ~ X1 + X2, fx$data),
                                supplied))
  expect_identical(fit$smc$max_steps, drbayes_control()$newton_steps)
  expect_true(is.na(fit$smc$n_ancestors_min))
  expect_false(is.na(fit$smc$ess))
})

test_that("collapsed weights are not reported as convergence", {
  # One draw a hair below zero and the rest well above it. As lambda falls the
  # weighted mean of B_n tends to that one draw, so the moment condition is met
  # in the limit by putting the whole posterior on a single point. Nothing
  # rejuvenates the cloud on this path, so the effective sample size is the
  # only thing standing between that and a reported answer.
  BB <- c(-1e-8, seq(0.30, 0.60, length.out = 999))
  fx <- np_moment_fixture(BB)
  expect_equal(moment_np(fx$betas.ps, fx$d), BB)

  set.seed(1)
  fit <- withCallingHandlers(
    run_tilting_np(fx$betas.ps, fx$d, drbayes_control(lambda_max = 1e6)),
    warning = function(cond) invokeRestart("muffleWarning"))

  expect_lt(abs(fit$smc$B_mean), fit$smc$tol)
  expect_lt(fit$smc$ess, 2)
  expect_lt(fit$smc$n_distinct, 0.01 * length(BB))
  expect_false(fit$smc$converged)

  set.seed(1)
  expect_warning(
    run_tilting_np(fx$betas.ps, fx$d, drbayes_control(lambda_max = 1e6)),
    "effective sample size")

  # The floor is what rejects it, and lowering the floor accepts it again.
  set.seed(1)
  lenient <- run_tilting_np(fx$betas.ps, fx$d,
                            drbayes_control(lambda_max = 1e6,
                                            ess_frac = 1e-6))
  expect_true(lenient$smc$converged)
})

test_that("a moment condition of one sign is refused with usable advice", {
  # Algorithm 1 has no root to find here. The parametric path is told to fall
  # back on Algorithm 2, and this path cannot, so it must not say so.
  fx <- np_moment_fixture(seq(0.1, 1, length.out = 40))

  expect_error(run_tilting_np(fx$betas.ps, fx$d, drbayes_control()),
               "same sign at all 40 posterior draws")
  expect_error(run_tilting_np(fx$betas.ps, fx$d, drbayes_control()),
               "pass them as outcome.samples with method = \"smc\"")
})

test_that("an untilted fit keeps every draw of the fitted means", {
  # The constraint already holds at lambda = 0 here, which is Lemma 1: nothing
  # is reweighted, so the index is the identity and the two estimands agree.
  fx  <- np_moment_fixture(seq(-0.5, 0.5, length.out = 500))
  fit <- run_tilting_np(fx$betas.ps, fx$d, drbayes_control())

  expect_identical(fit$idx, seq_len(500L))
  expect_identical(fit$smc$lambda, 0)
  expect_identical(fit$smc$n_steps, 0L)
  expect_true(fit$smc$converged)
})

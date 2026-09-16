# The two tilting algorithms of Section 3.3, the sample pruning of Section
# 5.3.1, and the subclassification moment condition of Appendix F.

# Draws from the two models, plus everything the tilting engine needs. With
# correct = FALSE the outcome model omits X2, so only the propensity score
# model is correctly specified and the moment condition is genuinely non-zero.
tilting_fixture <- function(n = 300, mc = 1500, seed = 5, correct = TRUE) {
  set.seed(seed)
  dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  dat$A <- rbinom(n, 1, plogis(0.3 + 0.8 * dat$X1 - 0.5 * dat$X2))
  dat$Y <- 1 + 2 * dat$A + 0.9 * dat$X1 + 0.7 * dat$X2 + rnorm(n)

  md <- prepare_model_data(if (correct) Y ~ A + X1 + X2 else Y ~ A + X1,
                           A ~ X1 + X2, dat)
  otc <- bayes_lm(Y = md$Y, X = md$Z.lm[, -1, drop = FALSE], mc = mc,
                  chains = 1L)
  ps  <- bayes_logit(Y = md$A, X = md$Z.ps[, -1, drop = FALSE], mc = mc,
                     chains = 1L)

  keep <- seq(501, mc, by = 2)
  list(d = tilting_data(Z.lm = md$Z.lm, Z.ps = md$Z.ps, A = md$A, Y = md$Y,
                        inverse_link = function(x) x,
                        ps_inverse_link = stats::plogis,
                        ps_formula_text = "A ~ X1 + X2"),
       betas.otc = matrix(otc[keep, 1L, ], ncol = dim(otc)[3]),
       betas.ps  = matrix(ps[keep, 1L, ], ncol = dim(ps)[3]),
       md = md)
}

# A fixture built the other way round: one observation, an intercept-only
# propensity score pinned at one half and an intercept-only outcome model, so
# that equation (3.4) reduces to -2 beta at each draw and a test can hand the
# tilting whatever moment condition it wants to see it handle.
moment_fixture <- function(BB) {
  list(d = tilting_data(Z.lm = matrix(1, 1L, 1L), Z.ps = matrix(1, 1L, 1L),
                        A = 1, Y = 0,
                        inverse_link = function(x) x,
                        ps_inverse_link = stats::plogis,
                        ps_formula_text = "A ~ 1"),
       betas.ps  = matrix(0, nrow = length(BB), ncol = 1L),
       betas.otc = matrix(-BB / 2, ncol = 1L))
}

# The g-computation contrast of equation (3.7), one value per draw.
ate_draws <- function(betas, md) {
  rowMeans(tcrossprod(betas, md$Z.lm1)) - rowMeans(tcrossprod(betas, md$Z.lm0))
}

# Several of the paths below warn more than once, and each warning is about a
# different fact, so they are collected rather than matched one at a time.
with_warnings <- function(expr) {
  warns <- character()
  value <- withCallingHandlers(expr, warning = function(cond) {
    warns <<- c(warns, conditionMessage(cond))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warns)
}

test_that("the moment condition is the mean of equation (3.4)", {
  skip_on_cran()
  fx <- tilting_fixture()

  # Equation (3.4) written out, which is not the expression moment_ipw()
  # evaluates: it forms the inverse probability weight from the linear
  # predictor instead, as one exponential.
  ps <- stats::plogis(tcrossprod(fx$d$Z.ps, fx$betas.ps))
  mu <- tcrossprod(fx$d$Z.lm, fx$betas.otc)
  expect_equal(moment_ipw(fx$betas.ps, fx$betas.otc, fx$d),
               colMeans((fx$d$A - ps) * (fx$d$Y - mu) / (ps * (1 - ps))),
               tolerance = 1e-10)
})

test_that("the logistic weight is the weight of equation (3.4)", {
  eta <- seq(-15, 15, by = 0.05)

  for (a in c(0, 1)) {
    d <- tilting_data(Z.lm = matrix(1, 1L, 1L), Z.ps = matrix(1, 1L, 1L),
                      A = a, Y = 0,
                      inverse_link = function(x) x,
                      ps_inverse_link = stats::plogis,
                      ps_formula_text = "A ~ 1")
    ps <- stats::plogis(eta)
    expect_equal(ipw_weight(d)(eta), (a - ps) / (ps * (1 - ps)),
                 tolerance = 1e-10)
  }
})

test_that("a link that is not the logistic one keeps the general weight", {
  eta <- seq(-4, 4, by = 0.05)
  d <- tilting_data(Z.lm = matrix(1, 1L, 1L), Z.ps = matrix(1, 1L, 1L),
                    A = 1, Y = 0,
                    inverse_link = function(x) x,
                    ps_inverse_link = stats::pnorm,
                    ps_formula_text = "A ~ 1")

  # A probit propensity score is not weighted by the logistic identity, so it
  # has to be recognised as a different link and weighted the general way.
  ps <- stats::pnorm(eta)
  expect_identical(ipw_weight(d)(eta), (1 - ps) / (ps * (1 - ps)))

  expect_true(is_logistic(stats::plogis))
  expect_true(is_logistic(function(eta) stats::plogis(eta)))
  expect_true(is_logistic(function(eta) switch("logit",
                                               logit  = stats::plogis(eta),
                                               probit = stats::pnorm(eta))))
  expect_false(is_logistic(stats::pnorm))
  expect_false(is_logistic(function(eta) eta))
  expect_false(is_logistic(function(eta) stats::plogis(eta / 2)))
})

test_that("draw_blocks partitions the draws in order", {
  # Whole matrices of observations by draws are what the moment condition
  # cannot afford, so the block shrinks as the sample grows and never vanishes.
  expect_identical(draw_blocks(10L, 4L, max_cells = 12), list(1:3, 4:6, 7:9,
                                                              10L))
  expect_identical(draw_blocks(3L, 1e6, max_cells = 100), list(1L, 2L, 3L))
  expect_identical(draw_blocks(10L, 1L), list(1:10))
  expect_identical(unlist(draw_blocks(1000L, 7L, max_cells = 30)), 1:1000)
})

test_that("the blocking moves the moment condition by no more than rounding", {
  skip_on_cran()
  fx <- tilting_fixture()
  S  <- nrow(fx$betas.ps)
  half <- seq_len(S / 2)

  # Equal, not identical, and the difference matters. Each draw's value is a
  # mean over observations that no other draw enters, so the reduction is
  # blocking-independent; but the fitted values it averages come from one
  # matrix product per block, and a BLAS may reach a different kernel for a
  # differently shaped product. Measured in plain R on a reference BLAS,
  # tcrossprod(Z, B) and the same product taken in two halves differ in 1680 of
  # 150000 entries, by at most 1.8e-15, with no compiled code involved. Which
  # blocks a fit uses is fixed, so a fit is reproducible; what cannot be
  # claimed is that two different partitions agree bit for bit.
  expect_equal(
    moment_ipw(fx$betas.ps, fx$betas.otc, fx$d),
    c(moment_ipw(fx$betas.ps[half, ], fx$betas.otc[half, ], fx$d),
      moment_ipw(fx$betas.ps[-half, ], fx$betas.otc[-half, ], fx$d)),
    tolerance = 1e-13)
})

test_that("a propensity score of 0 or 1 is refused wherever it appears", {
  d <- tilting_data(Z.lm = matrix(1, 1L, 1L), Z.ps = matrix(1, 1L, 1L),
                    A = 1, Y = 0,
                    inverse_link = function(x) x,
                    ps_inverse_link = stats::plogis,
                    ps_formula_text = "A ~ X1")

  # The guard reads the extremes of the linear predictor rather than every
  # propensity score, so what has to hold is that it refuses exactly the sets
  # of linear predictors whose scores leave (0, 1). plogis() reaches 1 just
  # above eta = 36.7 and 0 just below eta = -709.8, and the grids below
  # straddle both.
  grids <- list(seq(-700, 36, length.out = 101),
                seq(-700, 36.7, length.out = 101),
                seq(-700, 36.8, length.out = 101),
                seq(-709, 30, length.out = 101),
                seq(-710, 30, length.out = 101),
                c(0.5, NA), c(0.5, NaN), c(0.5, Inf), c(0.5, -Inf),
                matrix(seq(-4, 40, length.out = 100), 10L))
  for (eta in grids) {
    ps <- stats::plogis(eta)
    violated <- any(!is.finite(ps)) || any(ps <= 0 | ps >= 1)
    if (violated) {
      expect_error(check_positivity(eta, d), "positivity violation")
    } else {
      expect_null(check_positivity(eta, d))
    }
  }
})

test_that("a separating covariate stops the tilting with the same advice", {
  set.seed(3)
  n  <- 200
  X1 <- stats::rnorm(n)
  A  <- as.numeric(X1 > 0)
  d  <- tilting_data(Z.lm = cbind(1, A, X1), Z.ps = cbind(1, X1), A = A,
                     Y = 1 + 2 * A + X1 + stats::rnorm(n),
                     inverse_link = function(x) x,
                     ps_inverse_link = stats::plogis,
                     ps_formula_text = "A ~ X1")

  # A slope this steep on a covariate that separates the treatment groups puts
  # the fitted score of 60 of the 200 observations at 1, where equation (3.4)
  # has no finite weight. Written as a division it shows up there as 0/0; the
  # weight moment_ipw() evaluates has no division in it and returns an
  # unremarkable 1.4 instead, so nothing but this guard would notice.
  betas.ps  <- matrix(c(0, 60, 0.1, 61), nrow = 2L, byrow = TRUE)
  betas.otc <- matrix(c(1, 2, 1, 1.1, 1.9, 0.9), nrow = 2L, byrow = TRUE)

  expect_error(moment_ipw(betas.ps, betas.otc, d),
               "This is a positivity violation")
  expect_error(moment_ipw(betas.ps, betas.otc, d),
               "range(fitted(glm(A ~ X1, binomial, data)))", fixed = TRUE)
  expect_error(run_tilting(betas.ps, betas.otc, d, drbayes_control()),
               "positivity violation")
  expect_error(subclass_weights(betas.ps, d, 2L), "positivity violation")
})

test_that("importance sampling and the sweep reach the same tilted posterior", {
  skip_on_cran()
  fx   <- tilting_fixture()
  ctrl <- drbayes_control(n_steps = 200)

  set.seed(1)
  is  <- run_tilting(fx$betas.ps, fx$betas.otc, fx$d, ctrl, method = "is")
  set.seed(1)
  smc <- run_tilting(fx$betas.ps, fx$betas.otc, fx$d, ctrl, method = "smc")

  expect_true(is$smc$converged)
  expect_true(smc$smc$converged)
  expect_lt(abs(is$smc$B_mean), is$smc$tol)
  expect_lt(abs(smc$smc$B_mean), smc$smc$tol)

  # Algorithm 1 reweights the draws in place while Algorithm 2 moves them along
  # a grid, so the two reach the constraint at quite different values of the
  # tilting parameter and only the tilted posterior itself is comparable.
  untilted <- ate_draws(fx$betas.otc, fx$md)
  mcse     <- stats::sd(untilted) / sqrt(length(untilted))
  ate_is   <- mean(ate_draws(is$betas.otc, fx$md))
  ate_smc  <- mean(ate_draws(smc$betas.otc, fx$md))

  expect_lt(abs(ate_is - ate_smc), 5 * mcse)
  expect_lt(abs(stats::sd(ate_draws(is$betas.otc, fx$md)) /
                  stats::sd(ate_draws(smc$betas.otc, fx$md)) - 1), 0.25)
})

test_that("the tilting diagnostics carry the effective sample size", {
  skip_on_cran()
  fx   <- tilting_fixture()
  ctrl <- drbayes_control(n_steps = 200)
  S    <- nrow(fx$betas.ps)

  set.seed(1)
  is <- run_tilting(fx$betas.ps, fx$betas.otc, fx$d, ctrl, method = "is")
  set.seed(1)
  smc <- run_tilting(fx$betas.ps, fx$betas.otc, fx$d, ctrl, method = "smc")

  expect_named(is$smc, c("lambda", "B_mean", "B_mean_weighted", "tol",
                         "n_steps", "max_steps", "converged", "ess",
                         "n_ancestors_min", "n_distinct"))
  # A well behaved reweighting keeps nearly every draw, and the sweep's
  # resampling keeps roughly 1 - 1/e of the particles distinct.
  expect_gt(is$smc$ess, 0.5 * S)
  expect_lte(is$smc$ess, S)
  expect_gt(smc$smc$n_ancestors_min, 0.4 * S)
  expect_lte(smc$smc$n_ancestors_min, S)

  # Each method reports the degeneracy measure it has and no other. The sweep
  # resamples to equal weights at every step, so an importance weight
  # effective sample size would say nothing about it; and its kernel moves
  # every particle, so counting distinct ones would always return S.
  expect_true(is.na(smc$smc$ess))
  expect_true(is.na(smc$smc$n_distinct))
  expect_true(is.na(is$smc$n_ancestors_min))
  expect_identical(is$smc$max_steps, ctrl$newton_steps)
  expect_identical(smc$smc$max_steps, ctrl$n_steps)
})

test_that("the moment condition reported is the one at the draws returned", {
  skip_on_cran()
  fx   <- tilting_fixture()
  ctrl <- drbayes_control(n_steps = 200)

  for (method in c("is", "smc")) {
    set.seed(1)
    fit <- run_tilting(fx$betas.ps, fx$betas.otc, fx$d, ctrl, method = method)

    # Whatever the algorithm solved for internally, B_mean has to describe the
    # object the caller is handed, or the number reported is not the number
    # the user's credible intervals were built from.
    expect_equal(fit$smc$B_mean,
                 mean(moment_ipw(fit$betas.ps, fit$betas.otc, fx$d)),
                 tolerance = 1e-12)
  }

  set.seed(1)
  is <- run_tilting(fx$betas.ps, fx$betas.otc, fx$d, ctrl, method = "is")
  # Algorithm 1 solves equation (3.9) in the weights; the resampled draws
  # approximate that, and systematic resampling keeps the difference well
  # inside the tolerance so the two answers agree.
  expect_false(is.na(is$smc$B_mean_weighted))
  expect_lt(abs(is$smc$B_mean - is$smc$B_mean_weighted), 0.5 * is$smc$tol)

  set.seed(1)
  smc <- run_tilting(fx$betas.ps, fx$betas.otc, fx$d, ctrl, method = "smc")
  expect_true(is.na(smc$smc$B_mean_weighted))
})

test_that("collapsed importance weights are not reported as convergence", {
  # One draw a hair below zero and the rest well above it. As lambda falls the
  # weighted mean of B_n tends to that one draw, so the moment condition is
  # met in the limit by putting the whole posterior on a single point.
  BB <- c(-1e-8, seq(0.30, 0.60, length.out = 4999))
  fx <- moment_fixture(BB)
  expect_equal(moment_ipw(fx$betas.ps, fx$betas.otc, fx$d), BB)

  set.seed(1)
  run <- with_warnings(
    run_tilting(fx$betas.ps, fx$betas.otc, fx$d,
                drbayes_control(lambda_max = 1e6), method = "is"))
  smc <- run$value$smc

  # The moment condition on its own says yes, and the answer is still no.
  expect_lt(abs(smc$B_mean), smc$tol)
  expect_lt(smc$ess, 2)
  expect_false(smc$converged)
  expect_match(paste(run$warnings, collapse = " "), "effective sample size")

  # Raising the floor is enough to reject a merely thin sample as well.
  set.seed(1)
  strict <- suppressWarnings(
    run_tilting(fx$betas.ps, fx$betas.otc, fx$d,
                drbayes_control(lambda_max = 1e6, ess_frac = 1e-6),
                method = "is"))
  expect_true(strict$smc$converged)
})

test_that("importance sampling keeps the tilting parameter inside lambda_max", {
  BB <- c(-1e-8, seq(0.30, 0.60, length.out = 4999))
  fx <- moment_fixture(BB)

  for (bound in c(0.5, 10)) {
    set.seed(1)
    run <- with_warnings(
      run_tilting(fx$betas.ps, fx$betas.otc, fx$d,
                  drbayes_control(lambda_max = bound), method = "is"))
    # Section 5.3.1 treats a large tilting parameter as the problem, not as
    # the answer, so the bound holds and the caller is told it was reached.
    expect_equal(run$value$smc$lambda, -bound)
    expect_false(run$value$smc$converged)
    expect_match(paste(run$warnings, collapse = " "), "lambda_max")
  }
})

test_that("the Newton iteration is capped separately from the sweep", {
  BB   <- c(-1e-8, seq(0.30, 0.60, length.out = 4999))
  fx   <- moment_fixture(BB)
  ctrl <- drbayes_control(lambda_max = 1e6, newton_steps = 1L, n_steps = 500L)

  set.seed(1)
  fit <- suppressWarnings(
    run_tilting(fx$betas.ps, fx$betas.otc, fx$d, ctrl, method = "is"))
  expect_identical(fit$smc$n_steps, 1L)
  expect_identical(fit$smc$max_steps, 1L)

  # One iteration from lambda = 0 is the step newton_lambda() takes, because
  # the weights are still uniform there and the two ridged denominators are
  # the same average of B_n squared. A ridge on a sum instead would put them
  # orders of magnitude apart.
  expect_equal(fit$smc$lambda,
               damped_newton(0, newton_lambda(BB, ctrl$ridge), BB,
                             abs(mean(BB)), ctrl$lambda_max))
})

test_that("importance sampling warns when it cannot solve equation (3.9)", {
  skip_on_cran()
  fx <- tilting_fixture(correct = FALSE)

  set.seed(2)
  run <- with_warnings(
    run_tilting(fx$betas.ps, fx$betas.otc, fx$d,
                drbayes_control(n_steps = 200), method = "is"))

  # Only the propensity score model is right here, so the constraint sits far
  # from where the untilted draws are and Algorithm 1 runs to the bound. What
  # matters is that the fit is not passed off as doubly robust.
  expect_false(run$value$smc$converged)
  expect_gt(abs(run$value$smc$B_mean), run$value$smc$tol)
  expect_match(paste(run$warnings, collapse = " "),
               "did not satisfy the moment condition")
})

test_that("importance sampling refuses a moment condition of one sign", {
  betas <- matrix(0, nrow = 40, ncol = 2)
  BB    <- seq(0.1, 1, length.out = 40)

  # Equation (3.9) weights the draws it has; when none of them puts the moment
  # condition below zero, no weighting of them averages to zero either. That
  # is the drawback of Algorithm 1 that Section 3.3 names, so the refusal has
  # to point at Algorithm 2 rather than merely fail.
  expect_error(tilt_is(betas, betas, BB, drbayes_control(), tol = 1e-3),
               "same sign at all 40 posterior draws")
  expect_error(tilt_is(betas, betas, BB, drbayes_control(), tol = 1e-3),
               "method = \"smc\"")
})

test_that("pruning at prune.q = 0 leaves the sweep untouched", {
  skip_on_cran()
  fx <- tilting_fixture()
  as_control <- function(...) drbayes_control(n_steps = 200, ...)

  set.seed(9)
  plain <- run_tilting(fx$betas.ps, fx$betas.otc, fx$d,
                       as_control(pruning = FALSE))
  set.seed(9)
  inert <- run_tilting(fx$betas.ps, fx$betas.otc, fx$d,
                       as_control(pruning = TRUE, prune.rule = "quantile",
                                  prune.q = 0))
  set.seed(9)
  pruned <- run_tilting(fx$betas.ps, fx$betas.otc, fx$d,
                        as_control(pruning = TRUE, prune.rule = "quantile",
                                   prune.q = 0.1))

  # Identical, not merely equal: pruning nothing must also consume the same
  # random numbers, otherwise the option is not provably inert at zero.
  expect_identical(plain, inert)
  expect_false(identical(plain, pruned))
})

test_that("pruning meets the constraint by deleting rather than by tilting", {
  skip_on_cran()
  fx <- tilting_fixture()
  sweep_with <- function(...) {
    set.seed(9)
    run_tilting(fx$betas.ps, fx$betas.otc, fx$d,
                drbayes_control(n_steps = 200, ...))
  }

  plain  <- sweep_with(pruning = FALSE)
  pruned <- sweep_with(pruning = TRUE, prune.rule = "quantile",
                       prune.q = 0.1)

  # Discarding a tenth of the cloud moves the posterior mean of B_n further
  # than an increment of the tilting parameter does, so the sweep ends almost
  # at once and at a much smaller lambda. This is the behaviour ?drbayes_control
  # describes under the pruning threshold, and it is pinned here because it is
  # surprising rather than because it is wrong.
  expect_true(plain$smc$converged)
  expect_true(pruned$smc$converged)
  expect_lt(pruned$smc$n_steps, plain$smc$n_steps)
  expect_lt(abs(pruned$smc$lambda), abs(plain$smc$lambda))

  # The cloud that survives is a truncated one, so its spread understates the
  # tilted posterior's rather than measuring it better.
  expect_lt(stats::sd(ate_draws(pruned$betas.otc, fx$md)),
            stats::sd(ate_draws(plain$betas.otc, fx$md)))
})

test_that("pruning at the shipped defaults is not a no-op", {
  skip_on_cran()
  # The sweep has to actually run, so the outcome model must be misspecified;
  # under a correct one Lemma 1 skips the sweep and nothing can prune.
  fx <- tilting_fixture(correct = FALSE)

  set.seed(9)
  plain <- run_tilting(fx$betas.ps, fx$betas.otc, fx$d, drbayes_control())
  set.seed(9)
  pruned <- run_tilting(fx$betas.ps, fx$betas.otc, fx$d,
                        drbayes_control(pruning = TRUE))

  # A threshold below the smallest weight the sweep ever produces would make
  # pruning bit-identical to leaving it off, so a user asking for section
  # 5.3.1 would silently get the plain sweep. The default has to be inside the
  # range where the rule can act.
  expect_false(identical(plain$betas.otc, pruned$betas.otc))
})

test_that("pruning refuses to leave fewer than two particles", {
  # exp(-800) underflows to zero, which is what pruning meets in practice once
  # the sweep takes a large step in the tilting parameter.
  quantile_at <- function(q) {
    drbayes_control(pruning = TRUE, prune.rule = "quantile", prune.q = q)
  }
  collapsed <- normalised_weights(c(0, -800, -800, -800))
  expect_error(prune_weights(collapsed, collapsed, quantile_at(0.5)),
               "Lower prune.q")

  spread <- normalised_weights(seq(0, 3, length.out = 20))
  kept   <- prune_weights(spread, spread, quantile_at(0.25))
  expect_identical(sum(kept == 0), 5L)
  expect_equal(sum(kept), 1)
  expect_identical(prune_weights(spread, spread, quantile_at(0)), spread)

  # The weight rule discards on an absolute threshold rather than a fixed
  # share, so it removes nothing at all from weights that are already even.
  even <- rep(1 / 20, 20)
  expect_identical(
    prune_weights(even, even, drbayes_control(pruning = TRUE, prune.w = 0.3)),
    even)
})

test_that("systematic resampling reproduces the weighted mean it targets", {
  set.seed(4)
  BB <- stats::rnorm(2000, mean = 0.2, sd = 0.5)
  ww <- normalised_weights(-0.4 * BB)

  draws <- vapply(1:200, function(seed) {
    set.seed(seed)
    mean(BB[systematic_resample(ww, BB)])
  }, numeric(1))
  multinomial <- vapply(1:200, function(seed) {
    set.seed(seed)
    mean(BB[sample.int(length(ww), length(ww), replace = TRUE, prob = ww)])
  }, numeric(1))

  # Both target the weighted mean, but a multinomial draw scatters around it
  # by about the tolerance the tilting is judged against, which would make
  # convergence a coin toss; the stratified sweep costs far less than that.
  expect_lt(abs(mean(draws) - sum(ww * BB)), 1e-3)
  expect_lt(stats::sd(draws), 0.2 * stats::sd(multinomial))
  expect_length(systematic_resample(ww, BB), length(ww))
})

test_that("the subclassification moment tracks the inverse weighted one", {
  skip_on_cran()
  fx <- tilting_fixture()
  ds <- fx$d
  ds$subclass_weight <- subclass_weights(fx$betas.ps, ds, 5L)

  ipw <- moment_ipw(fx$betas.ps, fx$betas.otc, fx$d)
  sub <- moment_subclass(fx$betas.ps, fx$betas.otc, ds)

  # Equation (F.2) is equation (3.4) with the stratum treated fraction in place
  # of the individual propensity score, so the two move together across draws
  # and have a similar spread.
  expect_gt(stats::cor(ipw, sub), 0.9)
  expect_lt(abs(stats::sd(sub) / stats::sd(ipw) - 1), 0.3)
  expect_lt(abs(mean(sub)), 0.5 * stats::sd(sub))

  # Where the constraint is genuinely violated the two also agree on how far
  # from zero it is, which is the claim of Appendix F.
  wrong <- tilting_fixture(correct = FALSE)
  dw    <- wrong$d
  dw$subclass_weight <- subclass_weights(wrong$betas.ps, dw, 5L)
  ipw_w <- mean(moment_ipw(wrong$betas.ps, wrong$betas.otc, wrong$d))
  sub_w <- mean(moment_subclass(wrong$betas.ps, wrong$betas.otc, dw))

  expect_identical(sign(sub_w), sign(ipw_w))
  expect_gt(sub_w / ipw_w, 0.5)
  expect_lt(sub_w / ipw_w, 2)
})

test_that("the sweep holds the propensity score strata fixed", {
  skip_on_cran()
  fx   <- tilting_fixture()
  ctrl <- drbayes_control(n_steps = 200, moment = "subclass", n_subclass = 5L)

  set.seed(1)
  fit <- run_tilting(fx$betas.ps, fx$betas.otc, fx$d, ctrl)
  expect_true(fit$smc$converged)

  # The strata built from the untilted draws must still be the ones the sweep
  # was working with at its last step; restratifying as the particles moved
  # would leave the reported mean of B_n unreproducible.
  ds <- fx$d
  ds$subclass_weight <- subclass_weights(fx$betas.ps, ds, 5L)
  expect_equal(mean(moment_subclass(fit$betas.ps, fit$betas.otc, ds)),
               fit$smc$B_mean)
})

test_that("a stratum without both treatment arms is refused", {
  set.seed(3)
  n  <- 200
  X1 <- rnorm(n)
  A  <- as.numeric(X1 > 0)
  d  <- tilting_data(Z.lm = cbind(1, A, X1), Z.ps = cbind(1, X1), A = A,
                     Y = 1 + 2 * A + X1 + rnorm(n),
                     inverse_link = function(x) x,
                     ps_inverse_link = stats::plogis,
                     ps_formula_text = "A ~ X1")
  betas.ps  <- matrix(c(0, 0.5, 0.05, 0.55), nrow = 2L, byrow = TRUE)
  betas.otc <- matrix(c(1, 2, 1, 1.1, 1.9, 0.9), nrow = 2L, byrow = TRUE)

  # Treatment is decided by X1 here, so the propensity score strata at the
  # bottom of the range hold controls only and equation (F.2) has no treated
  # unit to weight.
  expect_error(subclass_weights(betas.ps, d, 5L),
               "stratum 1 of 5.*no treated units")
  expect_error(
    run_tilting(betas.ps, betas.otc, d,
                drbayes_control(moment = "subclass", n_subclass = 5L)),
    "positivity violation within the stratum"
  )
})

test_that("the new tilting settings are validated", {
  expect_identical(drbayes_control()$newton_steps, 100L)
  expect_identical(drbayes_control()$ess_frac, 0.1)
  expect_identical(drbayes_control(newton_steps = 5)$newton_steps, 5L)

  expect_error(drbayes_control(newton_steps = 0), "Newton iterations")
  expect_error(drbayes_control(newton_steps = 2.5), "newton_steps")
  expect_error(drbayes_control(ess_frac = 0), "ess_frac")
  expect_error(drbayes_control(ess_frac = 1.5), "ess_frac")
  expect_error(drbayes_control(ess_frac = "half"), "ess_frac")
  expect_identical(drbayes_control(ess_frac = 1)$ess_frac, 1)
})

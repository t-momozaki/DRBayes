# Behaviour of drbayes_pc that the pre-fix code got wrong: the propensity score
# link, the stopping rule, the diagnostics it reports, and argument validation.

sim_data <- function(n = 200, seed = 7, ps_link = "logit") {
  set.seed(seed)
  d <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  eta <- 0.3 + 0.7 * d$X1 - 0.4 * d$X2
  d$A <- rbinom(n, 1, if (ps_link == "probit") pnorm(eta) else plogis(eta))
  d$Y <- 1 + 2 * d$A + 0.8 * d$X1 + 0.5 * d$X2 + rnorm(n)
  d
}

fit_quietly <- function(...) {
  utils::capture.output(res <- suppressWarnings(drbayes_pc(..., verbose = FALSE)), file = nullfile())
  res
}

test_that("the propensity score link is honoured and must match its sampler", {
  skip_on_cran()
  d <- sim_data(ps_link = "probit")

  set.seed(1)
  fit <- fit_quietly(Y ~ A + X1 + X2, A ~ X1 + X2, d, family = "gaussian",
                     ps.link = "probit", outcome.model = bayes_lm,
                     ps.model = bayes_probit, mc = 800, bn = 200, thin = 2)
  expect_identical(fit$ps.link, "probit")

  # Probit draws pushed through the logistic inverse link would solve a
  # different constraint, so the mismatch is refused rather than tolerated.
  expect_error(
    fit_quietly(Y ~ A + X1 + X2, A ~ X1 + X2, d, family = "gaussian",
                ps.link = "logit", outcome.model = bayes_lm, ps.model = bayes_probit,
                mc = 400, bn = 100, thin = 2),
    "probit scale"
  )
  expect_error(
    fit_quietly(Y ~ A + X1 + X2, A ~ X1 + X2, d, family = "gaussian",
                link = "identity", outcome.model = bayes_logit, ps.model = bayes_logit,
                mc = 400, bn = 100, thin = 2),
    "logit scale"
  )
})

test_that("the sweep reports whether it met the moment condition", {
  skip_on_cran()
  d <- sim_data()

  set.seed(1)
  fit <- fit_quietly(Y ~ A + X1 + X2, A ~ X1 + X2, d, family = "gaussian",
                     outcome.model = bayes_lm, ps.model = bayes_logit,
                     mc = 800, bn = 200, thin = 2,
                     control = drbayes_control(keep_particles = TRUE))

  expect_named(fit$smc, c("lambda", "B_mean", "B_mean_weighted", "tol",
                          "n_steps", "max_steps", "converged", "ess",
                          "n_ancestors_min", "n_distinct"))
  expect_true(fit$smc$converged)
  # B_mean is recomputed from the draws that are returned, so the number
  # reported is the one describing the output rather than the weights that
  # produced it. The two differ under importance sampling.
  expect_equal(fit$smc$B_mean,
               mean(moment_ipw(fit$particles$ps, fit$particles$outcome,
                               fit$particles$model_data)))
  # Equation (3.6) requires the posterior mean of the moment condition to be
  # zero; converged means it is within tolerance, not merely sign-flipped.
  expect_lt(abs(fit$smc$B_mean), fit$smc$tol)
  expect_lte(fit$smc$n_steps, fit$smc$max_steps)
})

test_that("an unmeetable constraint warns instead of returning silently", {
  skip_on_cran()
  d <- sim_data()

  # A tolerance far below the Monte Carlo error of the mean is a legal request,
  # but no sweep can reach it, so the grid is exhausted.
  expect_warning(
    fit <- drbayes_pc(Y ~ A + X1 + X2, A ~ X1 + X2, d, family = "gaussian",
                      outcome.model = bayes_lm, ps.model = bayes_logit,
                      mc = 500, bn = 100, thin = 2, verbose = FALSE,
                      control = drbayes_control(tol = 1e-12)),
    "did not satisfy the moment condition"
  )
  expect_false(fit$smc$converged)
})

test_that("g.comp and pc average over the same empirical distribution", {
  skip_on_cran()
  d <- sim_data()

  set.seed(3)
  fit <- fit_quietly(Y ~ A + X1 + X2, A ~ X1 + X2, d, family = "gaussian",
                     outcome.model = bayes_lm, ps.model = bayes_logit,
                     mc = 800, bn = 200, thin = 2)

  expect_length(fit$g.comp, length(fit$pc))
  # Both are the correctly specified estimand here, so they should agree closely;
  # a Dirichlet-weighted g.comp against a 1/n pc would inflate the spread of one
  # of them and make the two incomparable.
  expect_lt(abs(mean(fit$g.comp) - mean(fit$pc)), 0.5)
  expect_lt(abs(sd(fit$g.comp) / sd(fit$pc) - 1), 0.5)
})

test_that("the g-computation averages the fitted values, not the design", {
  set.seed(21)
  n <- 130
  Z1 <- cbind(1, 1, matrix(rnorm(n * 3), n, 3))
  Z0 <- Z1
  Z0[, 2] <- 0
  # A draw count that is not a whole number of blocks, so that the last block
  # is a short one.
  betas <- matrix(rnorm(37 * 5), 37L, 5L)

  links <- list(identity = function(eta) eta,
                logit    = stats::plogis,
                probit   = stats::pnorm)

  for (nm in names(links)) {
    f <- links[[nm]]
    # Blocking the draws regroups nothing inside a row mean, so the answer is
    # the one the whole draws by observations matrix gives, bit for bit.
    expect_identical(fitted_row_means(betas, Z1, f),
                     rowMeans(f(tcrossprod(betas, Z1))), info = nm)
  }

  # Only the identity link lets the average over observations move inside the
  # link, which is why the other two cannot be reduced to a product against
  # the column means of the design matrix.
  expect_equal(fitted_row_means(betas, Z1, links$identity),
               drop(betas %*% colMeans(Z1)))
  for (nm in c("logit", "probit")) {
    expect_false(isTRUE(all.equal(fitted_row_means(betas, Z1, links[[nm]]),
                                  drop(links[[nm]](betas %*% colMeans(Z1))))),
                 info = nm)
  }

  # And the contrast of the two treatment values is what drbayes_pc() reports.
  ate <- fitted_row_means(betas, Z1, links$logit) -
    fitted_row_means(betas, Z0, links$logit)
  expect_identical(ate,
                   rowMeans(stats::plogis(tcrossprod(betas, Z1))) -
                     rowMeans(stats::plogis(tcrossprod(betas, Z0))))
})

test_that("pooling the chains keeps every draw where it belongs", {
  skip_on_cran()
  d <- sim_data(n = 60)

  set.seed(8)
  fit <- fit_quietly(Y ~ A + X1 + X2, A ~ X1 + X2, d, family = "gaussian",
                     outcome.model = bayes_lm, ps.model = bayes_logit,
                     mc = 400, bn = 100, thin = 2, chains = 3L,
                     control = drbayes_control(keep_particles = TRUE))

  set.seed(8)
  draws <- bayes_lm(Y = d$Y, X = cbind(A = d$A, X1 = d$X1, X2 = d$X2),
                    mc = 400, chains = 3L)
  kept  <- draws[seq(101L, 400L, by = 2L), , , drop = FALSE]

  # The pooled matrix is the array read column by column: chain 1 first, then
  # chain 2 under it, for each parameter in turn.
  expect_identical(fit$particles$outcome_untilted,
                   matrix(kept, ncol = dim(kept)[3],
                          dimnames = list(NULL, dimnames(kept)[[3]])))
  expect_identical(fit$particles$outcome_untilted[, 1L],
                   as.vector(kept[, , 1L]))
})

test_that("MCMC control arguments are validated up front", {
  d <- sim_data(n = 80)
  args <- list(Y ~ A + X1 + X2, A ~ X1 + X2, d, family = "gaussian",
               outcome.model = bayes_lm, ps.model = bayes_logit)

  expect_error(do.call(fit_quietly, c(args, list(mc = 100, bn = 100))),
               "mc must exceed bn")
  expect_error(do.call(fit_quietly, c(args, list(mc = 0))),
               "mc must be a single positive integer")
  expect_error(do.call(fit_quietly, c(args, list(mc = 500, bn = -1))),
               "bn must be a single non-negative integer")
  expect_error(do.call(fit_quietly, c(args, list(mc = 500, thin = 0))),
               "thin must be a single positive integer")
})

test_that("internal and external sampling retain the same draws", {
  skip_on_cran()
  d <- sim_data()

  set.seed(11)
  internal <- fit_quietly(Y ~ A + X1 + X2, A ~ X1 + X2, d, family = "gaussian",
                          outcome.model = bayes_lm, ps.model = bayes_logit,
                          mc = 800, bn = 200, thin = 2)

  set.seed(11)
  outcome_draws <- bayes_lm(Y = d$Y, X = cbind(A = d$A, X1 = d$X1, X2 = d$X2), mc = 800)
  ps_draws      <- bayes_logit(Y = d$A, X = cbind(X1 = d$X1, X2 = d$X2), mc = 800)
  external <- fit_quietly(Y ~ A + X1 + X2, A ~ X1 + X2, d, family = "gaussian",
                          outcome.samples = outcome_draws, ps.samples = ps_draws,
                          bn = 200, thin = 2)

  # The two branches must keep the same iteration indices, otherwise feeding a
  # sampler's own draws back in cannot reproduce the internal result.
  expect_length(external$g.comp, length(internal$g.comp))
})

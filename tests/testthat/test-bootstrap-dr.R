# Regression tests for the Bayesian bootstrap estimator. The interaction test
# fails against the pre-fix code by roughly -3 in a true ATE of 5.

heterogeneous_data <- function(n = 1500, seed = 99) {
  set.seed(seed)
  d <- data.frame(X1 = rnorm(n, mean = 1), X2 = rnorm(n))
  d$A <- rbinom(n, 1, plogis(0.5 * d$X1 + 0.3 * d$X2 - 0.5))
  # Effect modification by X1, which has a non-zero mean so that omitting the
  # interaction from the counterfactual shifts the ATE rather than cancelling.
  d$Y <- 1 + 2 * d$A + d$X1 + d$X2 + 3 * d$A * d$X1 + rnorm(n)
  d
}

true_ate <- function(d) 2 + 3 * mean(d$X1)

# The per-unit fitted contrast m1(X) - m0(X) at equal weights, which is what
# the counterfactual design matrices exist to produce.
unit_contrast <- function(outcome.formula, ps.formula, data) {
  md <- prepare_model_data(outcome.formula, ps.formula, data)
  beta <- fit_weighted(md$Z.lm, md$Y, rep(1, length(md$Y)), NULL, "outcome")
  drop((md$Z.lm1 - md$Z.lm0) %*% beta)
}

test_that("treatment interactions are recomputed in the counterfactuals", {
  skip_on_cran()
  d <- heterogeneous_data()

  set.seed(99)
  draws <- suppressMessages(
    drbayes_bb(Y ~ A + X1 + X2 + A:X1, A ~ X1 + X2, d,
               num_iterations = 300, verbose = FALSE))

  ci <- unname(quantile(draws, c(0.025, 0.975)))
  # Copying the observed treatment into the interaction columns collapses
  # m1(X) - m0(X) to the coefficient on A, which lands near 2 instead of 5.
  expect_gt(ci[1], 4)
  expect_lt(ci[2], 6)
  expect_lt(abs(mean(draws) - true_ate(d)), 0.3)
})

test_that("the fitted contrast is heterogeneous and follows the truth", {
  d <- heterogeneous_data()

  contrast <- unit_contrast(Y ~ A + X1 + X2 + A:X1, A ~ X1 + X2, d)

  # The true unit level effect is 2 + 3 * X1. Building the counterfactual
  # design matrix by overwriting the treatment column alone leaves the observed
  # treatment in the A:X1 column, which makes this contrast the same number for
  # every unit. Comparing only the average would not notice, because a single
  # constant can still land near the average effect.
  expect_gt(sd(contrast), 1)
  recovered <- unname(coef(lm(contrast ~ d$X1)))
  expect_equal(recovered, c(2, 3), tolerance = 0.05)
})

test_that("the estimate does not depend on the order of the formula terms", {
  d <- heterogeneous_data(n = 400)

  set.seed(5)
  a <- suppressMessages(drbayes_bb(Y ~ A + X1 + X2, A ~ X1 + X2, d,
                                   num_iterations = 100, verbose = FALSE))
  set.seed(5)
  b <- suppressMessages(drbayes_bb(Y ~ X1 + A + X2, A ~ X1 + X2, d,
                                   num_iterations = 100, verbose = FALSE))
  expect_equal(a, b)
})

test_that("a rank deficient design is refused rather than returning NA", {
  d <- heterogeneous_data(n = 300)
  d$X1dup <- d$X1
  d$X2dup <- d$X2
  d$Ybin <- rbinom(nrow(d), 1, plogis(0.5 * d$A + 0.3 * d$X1))

  # Each solver reports the model it was fitting. Naming the wrong one sends
  # the reader to a formula that is not the problem, and both solvers have
  # historically returned a vector of NAs here instead of stopping.
  expect_error(
    suppressMessages(drbayes_bb(Y ~ A + X1 + X1dup, A ~ X1, d,
                                num_iterations = 5, verbose = FALSE)),
    "The outcome design matrix is rank deficient: 4 columns but rank 3")
  expect_error(
    suppressMessages(drbayes_bb(Y ~ A + X1, A ~ X1 + X1dup, d,
                                num_iterations = 5, verbose = FALSE)),
    paste("The propensity score design matrix is rank deficient:",
          "3 columns but rank 2"))
  expect_error(
    suppressMessages(drbayes_bb(Ybin ~ A + X2 + X2dup, A ~ X1, d,
                                num_iterations = 5, verbose = FALSE)),
    "The outcome design matrix is rank deficient: 4 columns but rank 3")
})

test_that("a model that will not converge names itself in the error", {
  set.seed(3)
  n <- 200
  d <- data.frame(X = rnorm(n), Z = rnorm(n))
  d$A <- as.integer(d$X > 0)
  d$Y <- 1 + 2 * d$A + d$X + rnorm(n)
  d$Ysep <- as.integer(d$Z > 0)

  # X separates the treatment groups completely, so the score model has no
  # finite maximiser. Reporting an arbitrary stopping point as a fit is how a
  # runaway coefficient reaches the returned draws.
  expect_error(
    suppressWarnings(suppressMessages(
      drbayes_bb(Y ~ A + X, A ~ X, d, num_iterations = 2, verbose = FALSE))),
    "The propensity score model did not converge")

  # Z separates the binary outcome instead, so the failure is in the other fit.
  expect_error(
    suppressWarnings(suppressMessages(
      drbayes_bb(Ysep ~ A + Z, A ~ Z, d, num_iterations = 2,
                 verbose = FALSE))),
    "The outcome model did not converge")
})

# Treatment groups that are separated except for a handful of crossovers. The
# fitted scores then run to the edge of the unit interval without the fit
# failing to converge, which is what the trim argument exists for.
poor_overlap_data <- function(n = 200, seed = 3, crossovers = 6) {
  set.seed(seed)
  d <- data.frame(X = rnorm(n))
  d$A <- as.integer(d$X > 0)
  flipped <- sample(n, crossovers)
  d$A[flipped] <- 1L - d$A[flipped]
  d$Y <- 1 + 2 * d$A + d$X + rnorm(n)
  d
}

test_that("a propensity score at the boundary is reported, not propagated", {
  d <- poor_overlap_data()

  # Under a probit link the fitted score underflows to exactly zero here, so
  # the augmentation term would be Inf or NaN. The message has to say which
  # assumption failed and what to do, since the arithmetic itself is silent.
  set.seed(61)
  expect_error(
    suppressWarnings(suppressMessages(
      drbayes_bb(Y ~ A + X, A ~ X, d, num_iterations = 20,
                 ps.link = "probit", verbose = FALSE))),
    "positivity violation")

  set.seed(61)
  rescued <- suppressWarnings(suppressMessages(
    drbayes_bb(Y ~ A + X, A ~ X, d, num_iterations = 20, ps.link = "probit",
               trim = 0.01, verbose = FALSE)))
  expect_length(rescued, 20)
  expect_true(all(is.finite(rescued)))
})

test_that("trim acts exactly when a score is outside the interval", {
  set.seed(21)
  n <- 400
  d <- data.frame(X = rnorm(n))
  d$A <- rbinom(n, 1, plogis(3 * d$X))
  d$Y <- 1 + 2 * d$A + d$X + rnorm(n)

  # A threshold below every fitted score has to leave the draws untouched. A
  # threshold set where nothing can reach it is how an option ends up being
  # bit-identical to leaving it off while still looking as if it did something.
  set.seed(55)
  untrimmed <- suppressMessages(
    drbayes_bb(Y ~ A + X, A ~ X, d, num_iterations = 4, verbose = FALSE))
  set.seed(55)
  inert <- suppressMessages(
    drbayes_bb(Y ~ A + X, A ~ X, d, num_iterations = 4, trim = 1e-12,
               verbose = FALSE))
  expect_identical(inert, untrimmed)

  set.seed(55)
  trimmed <- suppressMessages(
    drbayes_bb(Y ~ A + X, A ~ X, d, num_iterations = 4, trim = 0.05,
               verbose = FALSE))
  expect_gt(max(abs(trimmed - untrimmed)), 0.01)

  # The reported count is the reason to trust the bias it trades for variance,
  # so it is checked against the same count taken from an independent fit.
  set.seed(55)
  expected_count <- sum(replicate(4, {
    weighted <- d
    xi <- rexp(n)
    weighted$bb_weight <- n * (xi / sum(xi))
    e <- predict(glm(A ~ X, data = weighted, weights = bb_weight,
                     family = quasibinomial()), type = "response")
    sum(e < 0.05 | e > 0.95)
  }))

  set.seed(55)
  reported <- capture.output(
    invisible(drbayes_bb(Y ~ A + X, A ~ X, d, num_iterations = 4, trim = 0.05,
                         verbose = TRUE)),
    type = "message")
  expect_match(reported,
               paste0(expected_count, " of ", 4 * n,
                      " fitted propensity scores were truncated"),
               fixed = TRUE, all = FALSE)

  expect_error(
    suppressMessages(drbayes_bb(Y ~ A + X, A ~ X, d, num_iterations = 5,
                                trim = 0.6, verbose = FALSE)),
    "trim must be")
})

test_that("missing values are dropped instead of poisoning every draw", {
  d <- heterogeneous_data(n = 300)
  d$X1[1:5] <- NA

  draws <- suppressMessages(
    drbayes_bb(Y ~ A + X1 + X2, A ~ X1 + X2, d,
               num_iterations = 50, verbose = FALSE))

  expect_length(draws, 50)
  expect_true(all(is.finite(draws)))

  expect_error(
    suppressMessages(drbayes_bb(Y ~ A + X1 + X2, A ~ X1 + X2, d,
                                num_iterations = 5, na.action = "na.fail",
                                verbose = FALSE)),
    "missing values")

  # Both fits have to see the same rows, so a variable that is missing in only
  # one of the two models still removes the observation from both.
  d$X2[10:20] <- NA
  complete <- sum(!is.na(d$X1) & !is.na(d$X2))
  expect_warning(
    used <- capture.output(
      invisible(drbayes_bb(Y ~ A + X1, A ~ X2, d, num_iterations = 2,
                           verbose = TRUE)),
      type = "message"),
    "different missing value patterns")
  expect_match(used, paste0("on ", complete, " observations"),
               fixed = TRUE, all = FALSE)
})

test_that("rows are matched by name, not by position", {
  d <- heterogeneous_data(n = 300)
  set.seed(17)
  # A subset carries the row names of the parent frame, so the names are no
  # longer 1:nrow. Treating them as positions would select different rows, or
  # none at all, without any error.
  scrambled <- d[sample(nrow(d), 200), ]
  scrambled$X1[c(3, 50, 120)] <- NA
  renumbered <- scrambled
  rownames(renumbered) <- NULL

  set.seed(8)
  a <- suppressMessages(drbayes_bb(Y ~ A + X1 + X2, A ~ X1, scrambled,
                                   num_iterations = 20, verbose = FALSE))
  set.seed(8)
  b <- suppressMessages(drbayes_bb(Y ~ A + X1 + X2, A ~ X1, renumbered,
                                   num_iterations = 20, verbose = FALSE))
  expect_identical(a, b)
})

test_that("binary outcomes run through the binomial path", {
  d <- heterogeneous_data(n = 400)
  d$Yb <- rbinom(nrow(d), 1, plogis(-0.3 + 1.2 * d$A + 0.6 * d$X1))

  draws <- suppressMessages(
    drbayes_bb(Yb ~ A + X1 + X2, A ~ X1 + X2, d,
               num_iterations = 100, verbose = FALSE))

  expect_length(draws, 100)
  # A risk difference has to lie in [-1, 1].
  expect_true(all(draws >= -1 & draws <= 1))

  # The contrast is a difference of probabilities, so it must be computed after
  # the inverse link, not on the linear predictor. The two differ by more than
  # rounding whenever the fitted probabilities are away from one half.
  set.seed(71)
  probability_scale <- suppressMessages(
    drbayes_bb(Yb ~ A + X1, A ~ X1, d, num_iterations = 200, verbose = FALSE))
  fit <- glm(Yb ~ A + X1, data = d, family = binomial())
  treated <- d
  control <- d
  treated$A <- 1
  control$A <- 0
  reference <- mean(predict(fit, treated, type = "response") -
                      predict(fit, control, type = "response"))
  expect_lt(abs(mean(probability_scale) - reference), 0.02)
  expect_gt(abs(reference - unname(coef(fit)["A"])), 0.1)
})

test_that("num_iterations is validated", {
  d <- heterogeneous_data(n = 100)
  expect_error(
    suppressMessages(drbayes_bb(Y ~ A + X1, A ~ X1, d, num_iterations = 0)),
    "num_iterations must be"
  )
  expect_error(
    suppressMessages(drbayes_bb(Y ~ A + X1, A ~ X1, d, num_iterations = 2.5)),
    "num_iterations must be"
  )
})

test_that("the closed form gaussian solver agrees with least squares", {
  set.seed(80)
  n <- 200
  Z <- cbind(1, matrix(rnorm(n * 3), n, 3))
  y <- drop(Z %*% c(1, 2, -1, 0.5)) + rnorm(n)
  w <- n * rexp(n)

  # The gaussian path takes a shortcut through a pivoted QR on the square root
  # weighted problem. If the square root were dropped, or applied to only one
  # side, the coefficients would still look reasonable but would answer a
  # different weighted least squares problem.
  closed_form <- fit_weighted(Z, y, w, NULL, "outcome")
  reweighted <- unname(glm.fit(Z, y, weights = w,
                               family = gaussian())$coefficients)
  expect_equal(closed_form, reweighted, tolerance = 1e-10)
  expect_equal(closed_form, unname(coef(lm(y ~ Z - 1, weights = w))),
               tolerance = 1e-10)

  # The engine multiplies the Dirichlet weights by the sample size before
  # fitting. That is presented as cosmetic, so neither solver may move when the
  # weights are rescaled.
  binary <- rbinom(n, 1, plogis(drop(Z %*% c(0, 1, -0.5, 0.2))))
  logistic <- quasibinomial(link = "logit")
  expect_equal(fit_weighted(Z, y, w / n, NULL, "outcome"), closed_form,
               tolerance = 1e-10)
  expect_equal(fit_weighted(Z, binary, w / n, logistic, "outcome"),
               fit_weighted(Z, binary, w, logistic, "outcome"),
               tolerance = 1e-8)
})

test_that("the bootstrap weights are Dirichlet(1, ..., 1)", {
  n <- 40
  draws <- 20000

  # The concentration of these weights is what sets the width of the whole
  # posterior, so what is pinned here is the distribution and not merely the
  # shape. A generator that divided by n rather than by the realised sum, or
  # that drew from the wrong shape parameter, still returns positive numbers of
  # the right length and would move the posterior spread silently.
  set.seed(70)
  w <- t(replicate(draws, dirichlet_weights(n)))

  expect_equal(dim(w), c(draws, n))
  expect_true(all(w > 0))
  expect_equal(rowSums(w), rep(1, draws), tolerance = 1e-12)

  # Marginally each coordinate is Beta(1, n - 1), and every coordinate has the
  # same law, so the first and the last are both checked against it.
  expect_gt(stats::ks.test(w[, 1], "pbeta", 1, n - 1)$p.value, 0.001)
  expect_gt(stats::ks.test(w[, n], "pbeta", 1, n - 1)$p.value, 0.001)

  # By the aggregation property of the Dirichlet, any k coordinates sum to a
  # Beta(k, n - k). That is a statement about the joint law rather than the
  # margins, so it is what rules out coordinates that are individually right
  # but wrongly dependent.
  expect_gt(stats::ks.test(rowSums(w[, 1:10]), "pbeta", 10, n - 10)$p.value,
            0.001)

  # Moments of Dirichlet(1, ..., 1): mean 1 / n, variance
  # (n - 1) / (n^2 (n + 1)), and covariance -1 / (n^2 (n + 1)) between any two
  # distinct coordinates. The covariance is negative because the coordinates are
  # tied together by summing to one; independent weights would put it at zero.
  expect_equal(mean(w[, 1]), 1 / n, tolerance = 0.03)
  expect_equal(var(w[, 1]), (n - 1) / (n^2 * (n + 1)), tolerance = 0.08)
  expect_equal(cov(w[, 1], w[, 2]), -1 / (n^2 * (n + 1)), tolerance = 0.5)
})

test_that("only the links the models were fitted with are inverted", {
  expect_equal(link_function("identity")(c(-1, 0, 2)), c(-1, 0, 2))
  expect_equal(link_function("logit")(0.4), plogis(0.4))
  expect_equal(link_function("probit")(0.4), pnorm(0.4))
  # An unsupported link must stop here rather than fall through to the identity
  # and report linear predictors as if they were probabilities.
  expect_error(link_function("cloglog"), "Unsupported link: cloglog")
})

test_that("the family and link are resolved from the outcome or refused", {
  continuous <- rnorm(20)
  binary <- rep(c(0, 1), 10)

  expect_message(resolve_family_link(continuous, NULL, NULL, "logit", TRUE),
                 "Auto-detected continuous outcome")
  expect_message(resolve_family_link(binary, NULL, NULL, "logit", TRUE),
                 "Auto-detected binary outcome")
  expect_silent(resolve_family_link(continuous, NULL, NULL, "logit", FALSE))

  expect_equal(resolve_family_link(binary, NULL, NULL, "probit", FALSE),
               list(family = "binomial", link = "logit", ps.link = "probit"))
  expect_equal(resolve_family_link(continuous, NULL, NULL, "logit", FALSE),
               list(family = "gaussian", link = "identity", ps.link = "logit"))

  # The link belongs to a family, and the propensity score link is separate
  # from the outcome link. Accepting a mismatch is how probit coefficients end
  # up being pushed through a logistic inverse.
  expect_error(
    resolve_family_link(continuous, "gaussian", "logit", "logit", FALSE),
    "Link \"logit\" is not supported for family \"gaussian\"")
  expect_error(
    resolve_family_link(binary, "binomial", "cloglog", "logit", FALSE),
    "Valid links: logit, probit")
  expect_error(
    resolve_family_link(continuous, "binomial", NULL, "logit", FALSE),
    "the outcome must be binary")
  expect_error(resolve_family_link(binary, NULL, NULL, "cauchit", FALSE),
               "should be one of")
})

test_that("the old dotted names still work and warn", {
  d <- heterogeneous_data(n = 300)

  set.seed(5)
  current <- suppressMessages(drbayes_bb(Y ~ A + X1, A ~ X1, d,
                                         num_iterations = 50, verbose = FALSE))

  # Every renamed export forwards to its replacement for one release. The
  # nleqslv variant forwards here too: it computed the same estimator to about
  # 1e-08, so there is nothing of its own left to keep.
  for (old in list(DRBayes.BB, DRBayes.BB.nleqslv)) {
    set.seed(5)
    expect_warning(
      forwarded <- suppressMessages(old(Y ~ A + X1, A ~ X1, d,
                                        num_iterations = 50, verbose = FALSE)),
      "deprecated"
    )
    expect_identical(forwarded, current)
  }

  expect_warning(
    bayes_draws <- B.LM(rnorm(50), matrix(rnorm(50), 50, 1), mc = 20),
    "deprecated")
  expect_identical(dim(bayes_draws), c(20L, 4L, 2L))
})

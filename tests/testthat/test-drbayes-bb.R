# Tests for the exported Bayesian bootstrap estimator of Saarela et al. (2016).
# The estimand tests come first because they are the ones that matter: this
# function does not target the superpopulation average treatment effect, and a
# later "fix" that made it do so would be a regression, not an improvement.

# Both working models are misspecified here: the true propensity score and the
# true outcome regression are quadratic in X while the formulas used below are
# linear. That is the only situation in which the doubly robust bias term
# survives, and therefore the only situation in which the estimand defined by
# the estimated propensity score and the superpopulation ATE come apart. The
# superpopulation ATE is 1 by construction.
mixed_estimand_data <- function(n = 3000, seed = 7) {
  set.seed(seed)
  x <- runif(n, -1.5, 1.5)
  d <- data.frame(X = x)
  d$A <- rbinom(n, 1, plogis(1.2 * x^2 + 0.4 * x - 1.2))
  d$Y <- 1 * d$A + 3 * x^2 + 0.5 * x + rnorm(n)
  d
}

# The estimator of section 4.3 evaluated at equal weights, reached through the
# formula interface and stats::predict() rather than through the package's
# design matrices. The Bayesian bootstrap posterior concentrates on this value,
# so it is an independent route to whatever drbayes_bb() is targeting.
dr_reference <- function(data, outcome.formula, ps.formula) {
  ps_fit  <- glm(ps.formula, data = data, family = binomial())
  out_fit <- lm(outcome.formula, data = data)
  treated <- data
  control <- data
  treated$A <- 1
  control$A <- 0
  e  <- as.numeric(predict(ps_fit, type = "response"))
  m1 <- as.numeric(predict(out_fit, newdata = treated))
  m0 <- as.numeric(predict(out_fit, newdata = control))
  mu <- as.numeric(predict(out_fit))
  y  <- data[[all.vars(outcome.formula)[1]]]
  psi <- m1 - m0 + (data$A - e) / (e * (1 - e)) * (y - mu)
  list(estimate = mean(psi), se = sd(psi) / sqrt(nrow(data)))
}

test_that("the target is the estimand induced by the estimated score", {
  skip_on_cran()
  d <- mixed_estimand_data()

  set.seed(101)
  draws <- suppressMessages(
    drbayes_bb(Y ~ A + X, A ~ X, d, num_iterations = 1200, verbose = FALSE))

  mixed <- dr_reference(d, Y ~ A + X, A ~ X)$estimate
  expect_lt(abs(mean(draws) - mixed), 0.02)

  # Section 4.3: "in general, this estimand is not consistent with the ATE".
  # Anything that quietly turned this function into an estimator of the
  # superpopulation ATE would pull the posterior down towards 1.
  interval <- unname(quantile(draws, c(0.005, 0.995)))
  expect_gt(interval[1], 2)
  expect_gt(mixed, 2)
})

test_that("a correct propensity score model reconciles the two estimands", {
  skip_on_cran()
  d <- mixed_estimand_data()

  # Same misspecified outcome model, but the propensity score model now matches
  # the data generating process. The documented caveat is that the two
  # estimands coincide then, so the posterior has to cover the true ATE of 1.
  set.seed(102)
  draws <- suppressMessages(
    drbayes_bb(Y ~ A + X, A ~ X + I(X^2), d, num_iterations = 1200,
               verbose = FALSE))

  interval <- unname(quantile(draws, c(0.005, 0.995)))
  expect_lt(interval[1], 1)
  expect_gt(interval[2], 1)
})

test_that("the posterior spread matches the estimator's sampling error", {
  skip_on_cran()
  d <- mixed_estimand_data()

  # With both models correct the influence function of the equal-weight
  # estimator is known, so its standard error is an independent yardstick for
  # the Bayesian bootstrap posterior. This is what pins the weights to
  # Dirichlet(1, ..., 1): weights with any other concentration would rescale
  # the posterior spread by a constant while leaving the location alone.
  set.seed(103)
  draws <- suppressMessages(
    drbayes_bb(Y ~ A + X + I(X^2), A ~ X + I(X^2), d, num_iterations = 1200,
               verbose = FALSE))

  ref <- dr_reference(d, Y ~ A + X + I(X^2), A ~ X + I(X^2))
  expect_lt(abs(mean(draws) - ref$estimate), 0.02)
  expect_gt(sd(draws) / ref$se, 0.85)
  expect_lt(sd(draws) / ref$se, 1.18)
})

# One Dirichlet draw of the estimator of section 4.3, written out with
# stats::glm() and stats::predict() instead of the package's prebuilt design
# matrices. Weights are drawn exactly as bayes_bootstrap_dr() draws them, so
# that both computations see the same weight vector.
bb_reference <- function(data, outcome.formula, ps.formula, family, link,
                         ps.link, num_iterations, seed) {
  n       <- nrow(data)
  treated <- data
  control <- data
  treated$A <- 1
  control$A <- 0
  y <- data[[all.vars(outcome.formula)[1]]]

  set.seed(seed)
  out <- numeric(num_iterations)
  for (iter in seq_len(num_iterations)) {
    # Written out rather than taken from the package, so that the reference
    # stays an independent route, but it must consume the stream identically or
    # the two computations would be comparing different weight vectors.
    xi <- rexp(n)
    xi <- xi / sum(xi)
    weighted <- data
    weighted$bb_weight <- n * xi

    ps_fit <- glm(ps.formula, data = weighted, weights = bb_weight,
                  family = quasibinomial(link = ps.link))
    out_fit <- if (family == "gaussian") {
      lm(outcome.formula, data = weighted, weights = bb_weight)
    } else {
      glm(outcome.formula, data = weighted, weights = bb_weight,
          family = quasibinomial(link = link))
    }

    e  <- as.numeric(predict(ps_fit, type = "response"))
    m1 <- as.numeric(predict(out_fit, newdata = treated, type = "response"))
    m0 <- as.numeric(predict(out_fit, newdata = control, type = "response"))
    mu <- as.numeric(predict(out_fit, type = "response"))

    out[iter] <- sum(xi * (m1 - m0)) +
      sum(xi * (y - mu) * (data$A - e) / (e * (1 - e)))
  }
  out
}

link_grid_data <- function(n = 300, seed = 7) {
  set.seed(seed)
  d <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  d$A <- rbinom(n, 1, plogis(0.6 * d$X1 - 0.4 * d$X2))
  d$Y <- 1 + 2 * d$A + d$X1 + 0.8 * d$A * d$X1 + rnorm(n)
  d$Yb <- rbinom(n, 1, plogis(-0.2 + 0.9 * d$A + 0.5 * d$X1))
  d
}

test_that("every family and link reproduces the estimator of section 4.3", {
  d <- link_grid_data()

  # A link that is accepted and then ignored is the failure this guards: the
  # propensity score link was once fixed at logit while probit coefficients
  # were handed to it, which returns a plausible number rather than an error.
  grid <- list(
    list(Y ~ A + X1 + X2 + A:X1, A ~ X1 + X2, "gaussian", "identity", "logit"),
    list(Y ~ A + X1 + X2, A ~ X1 + X2, "gaussian", "identity", "probit"),
    list(Yb ~ A + X1 + X2, A ~ X1 + X2, "binomial", "logit", "logit"),
    list(Yb ~ A + X1 + A:X1, A ~ X1 + X2, "binomial", "probit", "probit")
  )

  got <- list()
  want <- list()
  for (k in seq_along(grid)) {
    case <- grid[[k]]
    label <- paste(case[[3]], case[[4]], "outcome with", case[[5]], "score")
    set.seed(200 + k)
    got[[label]] <- suppressMessages(
      drbayes_bb(case[[1]], case[[2]], d, num_iterations = 3,
                 family = case[[3]], link = case[[4]], ps.link = case[[5]],
                 verbose = FALSE))
    want[[label]] <- bb_reference(d, case[[1]], case[[2]], case[[3]],
                                  case[[4]], case[[5]], num_iterations = 3,
                                  seed = 200 + k)
  }
  expect_equal(got, want, tolerance = 1e-8)
})

test_that("changing a link changes the answer", {
  d <- link_grid_data(n = 200)

  # A cheaper guard than the grid above, and one that survives on CRAN: an
  # argument that is accepted, recorded and then never consulted leaves the
  # draws bit-identical.
  set.seed(31)
  logit_ps <- suppressMessages(
    drbayes_bb(Y ~ A + X1, A ~ X1 + X2, d, num_iterations = 20,
               ps.link = "logit", verbose = FALSE))
  set.seed(31)
  probit_ps <- suppressMessages(
    drbayes_bb(Y ~ A + X1, A ~ X1 + X2, d, num_iterations = 20,
               ps.link = "probit", verbose = FALSE))
  expect_gt(max(abs(logit_ps - probit_ps)), 1e-4)

  set.seed(32)
  logit_y <- suppressMessages(
    drbayes_bb(Yb ~ A + X1, A ~ X1, d, num_iterations = 20,
               link = "logit", verbose = FALSE))
  set.seed(32)
  probit_y <- suppressMessages(
    drbayes_bb(Yb ~ A + X1, A ~ X1, d, num_iterations = 20,
               link = "probit", verbose = FALSE))
  expect_gt(max(abs(logit_y - probit_y)), 1e-4)
})

test_that("a saturated noiseless design returns the effect exactly", {
  # The outcome model fits without residual, so the augmentation term is
  # identically zero and every draw reduces to sum(xi) times the treatment
  # effect. Because the Dirichlet weights sum to one, that is 3 to machine
  # precision for every draw, whatever the weights happen to be. Unnormalised
  # weights, or a counterfactual design matrix that does not switch the
  # treatment, both break this.
  d <- data.frame(X = c(0, 0, 1, 1, 0, 0, 1, 1),
                  A = c(0, 1, 0, 1, 0, 1, 0, 1))
  d$Y <- 10 + 3 * d$A + 5 * d$X

  set.seed(41)
  draws <- suppressMessages(
    drbayes_bb(Y ~ A * X, A ~ X, d, num_iterations = 25, verbose = FALSE))

  expect_equal(draws, rep(3, 25), tolerance = 1e-12)
})

test_that("a factor treatment is read the same way as a 0/1 treatment", {
  set.seed(11)
  n <- 400
  d <- data.frame(X = rnorm(n))
  d$A <- rbinom(n, 1, plogis(0.5 * d$X))
  d$Y <- 1 + 2 * d$A + d$X + rnorm(n)
  d$A_factor <- factor(ifelse(d$A == 1, "treated", "control"),
                       levels = c("control", "treated"))
  d$A_flipped <- factor(d$A_factor, levels = c("treated", "control"))

  set.seed(2)
  numeric_draws <- suppressMessages(
    drbayes_bb(Y ~ A + X, A ~ X, d, num_iterations = 40, verbose = FALSE))
  set.seed(2)
  factor_draws <- suppressMessages(
    drbayes_bb(Y ~ A_factor + X, A_factor ~ X, d, num_iterations = 40,
               verbose = FALSE))
  expect_equal(factor_draws, numeric_draws)

  # Reversing the level order swaps which arm counts as treated, so the
  # contrast has to change sign. A mapping that ignored the level order would
  # return the same draws and no error.
  set.seed(2)
  flipped_draws <- suppressMessages(
    drbayes_bb(Y ~ A_flipped + X, A_flipped ~ X, d, num_iterations = 40,
               verbose = FALSE))
  expect_equal(flipped_draws, -numeric_draws)

  expect_message(
    invisible(drbayes_bb(Y ~ A_factor + X, A_factor ~ X, d,
                         num_iterations = 2, verbose = FALSE)),
    "Treating 'treated' as treated and 'control' as control")
})

test_that("a factor covariate keeps its contrasts in the counterfactuals", {
  set.seed(12)
  n <- 1200
  d <- data.frame(G = factor(sample(c("a", "b", "c"), n, TRUE)), X = rnorm(n))
  shift <- c(a = 0, b = 2, c = -1)[as.character(d$G)]
  d$A <- rbinom(n, 1, plogis(0.4 * d$X + 0.6 * shift))
  # The effect is 3 in group b and 2 elsewhere, so the counterfactual design
  # matrix has to recompute both dummy interaction columns, not just one.
  d$Y <- 1 + 2 * d$A + 3 * shift + d$X + 1 * d$A * (d$G == "b") + rnorm(n)
  true_ate <- 2 + mean(d$G == "b")

  md <- prepare_model_data(Y ~ A * G + X, A ~ G + X, d)
  beta <- fit_weighted(md$Z.lm, md$Y, rep(1, n), NULL, "outcome")
  contrast <- drop((md$Z.lm1 - md$Z.lm0) %*% beta)
  expect_equal(as.numeric(tapply(contrast, d$G, mean)), c(2, 3, 2),
               tolerance = 0.1)
  expect_lt(max(tapply(contrast, d$G, sd)), 1e-10)

  set.seed(4)
  draws <- suppressMessages(
    drbayes_bb(Y ~ A * G + X, A ~ G + X, d, num_iterations = 300,
               verbose = FALSE))
  interval <- unname(quantile(draws, c(0.025, 0.975)))
  expect_lt(interval[1], true_ate)
  expect_gt(interval[2], true_ate)
})

test_that("the progress report describes the data actually used", {
  d <- link_grid_data(n = 120)
  d$X1[1:5] <- NA

  messages <- capture.output(
    invisible(drbayes_bb(Y ~ A + X1, A ~ X1, d, num_iterations = 2,
                         verbose = TRUE)),
    type = "message")
  report <- grep("Bayesian bootstrap on", messages, value = TRUE)
  expect_length(report, 1)

  # The counts have to describe the rows the estimator kept, not the rows it
  # was handed, or the report understates how much of the data went missing.
  used <- d[!is.na(d$X1), ]
  expect_match(report, paste0("on ", nrow(used), " observations"),
               fixed = TRUE)
  expect_match(report, paste0("(", sum(used$A), " treated, ",
                              sum(1 - used$A), " control)"), fixed = TRUE)
  expect_match(report, "2 Dirichlet draws", fixed = TRUE)
  expect_match(messages, "Auto-detected continuous outcome", all = FALSE)
})

test_that("an unsupported family and link combination is refused", {
  d <- link_grid_data(n = 60)

  expect_error(
    suppressMessages(drbayes_bb(Y ~ A + X1, A ~ X1, d, num_iterations = 2,
                                family = "binomial", verbose = FALSE)),
    "outcome must be binary")
  expect_error(
    suppressMessages(drbayes_bb(Y ~ A + X1, A ~ X1, d, num_iterations = 2,
                                link = "logit", verbose = FALSE)),
    "not supported for family \"gaussian\"")
  expect_error(
    suppressMessages(drbayes_bb(Yb ~ A + X1, A ~ X1, d, num_iterations = 2,
                                link = "cloglog", verbose = FALSE)),
    "Valid links: logit, probit")
  expect_error(
    suppressMessages(drbayes_bb(Y ~ A + X1, A ~ X1, d, num_iterations = 2,
                                ps.link = "cauchit", verbose = FALSE)),
    "should be one of")
})

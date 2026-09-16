# Algorithm 3 and Algorithm A.1: the weight formula against a brute force
# product of densities, the order the paper puts the reweighting and the
# coupling in, the scale the shift is on and the size it comes out, what is
# reported when the reweighting spends the posterior, and the arithmetic that
# keeps a product over hundreds of observations representable.

sim_data <- function(n = 400, seed = 7) {
  set.seed(seed)
  d <- data.frame(X1 = stats::rnorm(n), X2 = stats::rnorm(n))
  d$A <- stats::rbinom(n, 1, stats::plogis(0.3 + 0.7 * d$X1 - 0.4 * d$X2))
  d$Y <- 1 + 2 * d$A + 0.8 * d$X1 + 0.5 * d$X2 + stats::rnorm(n)
  d
}

fit_quietly <- function(...) {
  utils::capture.output(
    res <- suppressWarnings(drbayes_pc(..., verbose = FALSE)),
    file = nullfile())
  res
}

# Fitting is by far the slowest part of these tests, so each fit is built once
# and every test that needs it gets the same object. A coarse lambda grid keeps
# the coupling of the reweighted draws, which every call now runs, quick.
cached_fit <- local({
  fits <- list()
  function(name, build) {
    if (is.null(fits[[name]])) fits[[name]] <<- build()
    fits[[name]]
  }
})

gaussian_fit <- function() {
  cached_fit("gaussian", function() {
    d <- sim_data()
    set.seed(1)
    list(data = d,
         fit  = fit_quietly(Y ~ A + X1 + X2, A ~ X1 + X2, d,
                            family = "gaussian", outcome.model = bayes_lm,
                            ps.model = bayes_logit, mc = 2000, bn = 500,
                            thin = 2,
                            control = drbayes_control(keep_particles = TRUE,
                                                      n_steps = 50)))
  })
}

# A tolerance this loose is already met at lambda = 0, so neither the fit nor
# the coupling of the reweighted draws tilts anything. The reweighted posterior
# of a coefficient is then exactly the original one shifted by -xi, which is
# what the two tests below measure.
untilted_fit <- function() {
  cached_fit("untilted", function() {
    d <- sim_data()
    set.seed(1)
    list(data = d,
         fit  = fit_quietly(Y ~ A + X1 + X2, A ~ X1 + X2, d,
                            family = "gaussian", outcome.model = bayes_lm,
                            ps.model = bayes_logit, mc = 2000, bn = 500,
                            thin = 2,
                            control = drbayes_control(keep_particles = TRUE,
                                                      n_steps = 50, tol = 1)))
  })
}

# Two thousand observations, of which about half are treated: the product of
# that many likelihood ratios is what a direct implementation cannot represent.
large_fit <- function() {
  cached_fit("large", function() {
    d <- sim_data(n = 2000, seed = 11)
    set.seed(2)
    list(data = d,
         fit  = fit_quietly(Y ~ A + X1 + X2, A ~ X1 + X2, d,
                            family = "gaussian", outcome.model = bayes_lm,
                            ps.model = bayes_logit, mc = 1000, bn = 250,
                            thin = 2, chains = 2, diagnostics = "none",
                            control = drbayes_control(keep_particles = TRUE,
                                                      n_steps = 40)))
  })
}

# The paper's own setting: an uncommon binary outcome, as the dementia
# incidence of Section 6.3 is, so the average treatment effect is a risk
# difference of 0.04 and xi is a shift of the same risk. About one treated
# observation in six has the outcome, which leaves room below 1 for a positive
# shift of a few hundredths but not for the whole of the paper's Tri(0, 0.5).
binary_fit <- function() {
  cached_fit("binary", function() {
    d <- sim_data(n = 500, seed = 5)
    d$Y <- stats::rbinom(500, 1,
                         stats::plogis(-2.2 + 0.5 * d$A + 0.5 * d$X1))
    set.seed(13)
    list(data = d,
         fit  = fit_quietly(Y ~ A + X1 + X2, A ~ X1 + X2, d,
                            family = "binomial", outcome.model = bayes_logit,
                            ps.model = bayes_logit, mc = 800, bn = 200,
                            thin = 2, chains = 2, diagnostics = "none",
                            control = drbayes_control(keep_particles = TRUE,
                                                      n_steps = 50)))
  })
}

# The moment condition (3.4) at the draws a sensitivity analysis reweights.
untilted_moment <- function(fit) {
  moment_ipw(fit$particles$ps_untilted, fit$particles$outcome_untilted,
             fit$particles$model_data)
}


test_that("a sensitivity parameter of zero leaves the posterior alone", {
  skip_on_cran()
  example <- untilted_fit()

  for (xi in list(0, function(n) rep(0, n))) {
    set.seed(21)
    sens <- drbayes_sensitivity(example$fit, xi = xi, M = 8)

    # A zero sensitivity parameter is strong ignorability, so every likelihood
    # ratio is one and every weight is 1/S exactly. Anything else means the
    # weight formula has picked up a term that does not belong to it.
    expect_equal(sens$log_weights, rep(0, sens$n_draws))
    expect_equal(sens$weights, rep(1 / sens$n_draws, sens$n_draws))
    expect_equal(sens$ess, sens$n_draws)

    # Uniform weights are what systematic resampling leaves alone, so the
    # particles that go into the coupling are the ones that came out of the
    # fit, and this fit does not tilt them. Only their order changes.
    expect_equal(sort(sens$ate), sort(example$fit$pc))
    expect_equal(sens$smc$lambda, 0)

    reported <- summary(sens)$estimands
    expect_equal(reported["sensitivity", ], reported["original", ],
                 ignore_attr = TRUE)

    # With every draw of xi at the same value the integral over g() is one
    # atom repeated, which is M effective atoms and not one.
    expect_equal(unname(sens$xi_ess), c(8, 8))
  }

  # The same must hold observation by observation, where the ratios are
  # averaged before rather than after the product over observations.
  set.seed(21)
  sens <- drbayes_sensitivity(example$fit, xi = 0, M = 8,
                              method = "per-observation")
  expect_equal(sens$weights, rep(1 / sens$n_draws, sens$n_draws))
  expect_equal(sort(sens$ate), sort(example$fit$pc))
})


test_that("the reweighted draws are coupled again, so they stay robust", {
  skip_on_cran()
  example <- gaussian_fit()
  BB      <- untilted_moment(example$fit)

  set.seed(22)
  sens <- drbayes_sensitivity(example$fit, xi = xi_triangular(0, 0.05),
                              M = 100)

  # Reweighting alone moves the particle cloud off the moment condition (3.4)
  # by many times the tolerance the fit itself was held to, and the effective
  # sample size does not notice: these weights are perfectly healthy.
  expect_gt(sens$ess, 0.5 * sens$n_draws)
  expect_gt(abs(sum(sens$weights * BB)), 5 * sens$smc$tol)

  # Coupling the reweighted draws is what puts it back, and it is reported.
  expect_true(sens$smc$converged)
  expect_true(sens$smc$converged_coupling)
  expect_lt(abs(sens$smc$B_mean), sens$smc$tol)
  expect_equal(sens$smc$tol, example$fit$smc$tol)
  expect_equal(length(sens$ate), sens$n_draws)
  expect_output(print(sens), "Moment condition")
  expect_output(print(sens), "Still doubly robust:   yes")
})


test_that("the verdict counts the sensitivity weights, not just the tilting", {
  # The coupling cannot see how few draws are behind the cloud it is given:
  # the resampling and the sweep leave every particle distinct again. So a
  # coupling that met the moment condition is only half the verdict.
  met <- list(lambda = -9, B_mean = 1e-4, tol = 2e-3, converged = TRUE)

  ample <- sensitivity_verdict(met, ess = 900, n_draws = 3000, ess_frac = 0.1)
  expect_true(ample$converged)
  expect_true(ample$converged_coupling)
  expect_equal(ample$ess_floor, 300)

  spent <- sensitivity_verdict(met, ess = 78, n_draws = 3000, ess_frac = 0.1)
  expect_false(spent$converged)
  expect_true(spent$converged_coupling)
  expect_equal(spent$ess_sensitivity, 78)

  # And a coupling that failed stays failed however many draws are behind it.
  missed <- sensitivity_verdict(list(converged = FALSE), 3000, 3000, 0.1)
  expect_false(missed$converged)
  expect_false(missed$converged_coupling)
})


test_that("a spent posterior is not reported as doubly robust", {
  skip_on_cran()
  example <- gaussian_fit()

  # E[xi] of 0.2 is nearly two posterior standard deviations of an effect of
  # 2.05, which leaves a few dozen effective draws out of three thousand. The
  # coupling of those draws may or may not meet the moment condition -- at this
  # effective sample size the mean it tests carries several times the
  # tolerance as Monte Carlo error -- and either way the analysis is not
  # doubly robust and says so.
  set.seed(77)
  sens <- suppressWarnings(
    drbayes_sensitivity(example$fit, xi = xi_triangular(0, 0.3), M = 200))

  expect_lt(sens$ess, sens$smc$ess_floor)
  expect_false(sens$smc$converged)
  expect_output(print(sens), "Still doubly robust:   NO")

  # The credible interval of a posterior that thin is not tabulated either.
  reported <- summary(sens)$estimands
  expect_true(all(is.na(reported["sensitivity", c("2.5%", "97.5%")])))
  expect_false(anyNA(reported["original", ]))
  expect_output(print(sens), "credible interval of the sensitivity row is not")
})


test_that("a positive sensitivity parameter moves the effect down by E[xi]", {
  skip_on_cran()
  example  <- untilted_fit()
  original <- mean(example$fit$pc)
  expect_equal(example$fit$smc$lambda, 0)

  # The bias enters as m_A(X; beta) + A xi, so the part of the observed
  # contrast that is confounding has to come off the treatment effect: the
  # posterior moves down by about E[xi], and up again when xi is negative.
  for (direction in c(1, -1)) {
    bound <- direction * 0.05
    set.seed(4)
    sens  <- drbayes_sensitivity(example$fit,
                                 xi = xi_triangular(min(0, bound),
                                                    max(0, bound)),
                                 M = 400)
    shift <- mean(sens$ate) - original

    expect_equal(sign(shift), -direction)
    expect_lt(abs(shift + mean(sens$xi)), 0.35 * abs(mean(sens$xi)))

    # Nothing is coupled here, so the two rows of the table are the same draws.
    expect_equal(sens$ate, sens$ate_gcomp)
  }
})


test_that("the coupling gives most of the shift back", {
  skip_on_cran()
  example  <- gaussian_fit()
  original <- mean(example$fit$pc)

  # The moment condition (3.4) is written with the unshifted mean function and
  # the propensity score of the same fit, both of which describe the data with
  # the confounding in it, so imposing it again pulls the treated fitted values
  # back to what the propensity score weighted data say. What the reweighting
  # moved by about E[xi], the coupling returns all but a tenth of.
  set.seed(71)
  sens <- drbayes_sensitivity(example$fit, xi = xi_triangular(0, 0.05),
                              M = 200)

  before <- mean(sens$ate_gcomp) - original
  after  <- mean(sens$ate) - original

  expect_lt(before, 0)
  expect_lt(abs(before + mean(sens$xi)), 0.35 * abs(mean(sens$xi)))
  expect_lt(abs(after), 0.5 * abs(before))
  expect_equal(rownames(summary(sens)$estimands),
               c("sensitivity", "reweighted", "original"))
})


test_that("the two algorithms differ in how they carry the spread of g()", {
  skip_on_cran()
  example <- untilted_fit()
  original_sd <- stats::sd(example$fit$pc)
  wide <- function(n) stats::runif(n, -0.1, 0.1)

  set.seed(4)
  common <- drbayes_sensitivity(example$fit, xi = wide, M = 400)
  set.seed(4)
  each   <- drbayes_sensitivity(example$fit, xi = wide, M = 100,
                                method = "per-observation")

  # Algorithm A.1 gives every observation the same bias, so the uncertainty
  # about that one bias is added to the posterior: the variance grows by the
  # variance of g().
  expect_equal(stats::sd(common$ate),
               sqrt(original_sd^2 + stats::var(as.numeric(common$xi))),
               tolerance = 0.05)

  # Algorithm 3 draws a bias for each observation independently, and hundreds
  # of independent biases average away, so a distribution centred on zero
  # leaves the posterior where it was.
  expect_equal(stats::sd(each$ate), original_sd, tolerance = 0.05)
  expect_gt(stats::sd(common$ate), stats::sd(each$ate))
})


test_that("the weights are the product of densities the paper writes down", {
  # A sample small enough that the products of densities themselves are
  # representable, so the two algorithms can be evaluated exactly as they are
  # written: raw densities over the whole sample, controls included, no logs
  # and no accumulator. Controls drop out of the code's product because the
  # shift is A xi, and this is what says that dropping them is right.
  set.seed(101)
  n <- 60
  Z <- cbind(1, stats::rbinom(n, 1, 0.5), stats::rnorm(n), stats::rnorm(n))
  A <- Z[, 2]
  Y <- drop(Z %*% c(1, 2, 0.8, 0.5)) + stats::rnorm(n)
  S <- 7
  betas <- matrix(stats::rnorm(S * 4, c(1, 2, 0.8, 0.5), 0.2), nrow = S,
                  byrow = TRUE)
  eta   <- Z[A == 1, , drop = FALSE] %*% t(betas)
  M     <- 6
  shared <- stats::runif(M, -0.4, 0.4)
  each   <- matrix(stats::runif(n * M, -0.4, 0.4), nrow = n)
  sigma  <- 1.3

  # Algorithm A.1: the sum over observations is inside the integral.
  by_draw <- function(f) vapply(seq_len(S), f, numeric(1))
  gaussian_common <- by_draw(function(s) {
    mu <- drop(Z %*% betas[s, ])
    log(mean(vapply(shared,
                    function(x) prod(stats::dnorm(Y, mu + A * x, sigma)),
                    numeric(1))) / prod(stats::dnorm(Y, mu, sigma)))
  })
  expect_equal(sensitivity_log_weights(eta, Y[A == 1], shared, "identity",
                                       sigma, "common")$log_weights,
               gaussian_common)

  # Algorithm 3: each observation is averaged over its own draws of xi first.
  gaussian_each <- by_draw(function(s) {
    mu <- drop(Z %*% betas[s, ])
    sum(log(vapply(seq_len(n), function(i) {
      mean(stats::dnorm(Y[i], mu[i] + A[i] * each[i, ], sigma)) /
        stats::dnorm(Y[i], mu[i], sigma)
    }, numeric(1))))
  })
  expect_equal(sensitivity_log_weights(eta, Y[A == 1],
                                       each[A == 1, , drop = FALSE],
                                       "identity", sigma,
                                       "per-observation")$log_weights,
               gaussian_each)

  # One residual standard deviation per draw, which is the default, has to
  # reach the right column of the matrix of linear predictors.
  sigmas <- stats::runif(S, 0.9, 1.6)
  expect_equal(
    sensitivity_log_weights(eta, Y[A == 1], shared, "identity", sigmas,
                            "common")$log_weights,
    by_draw(function(s) {
      mu <- drop(Z %*% betas[s, ])
      log(mean(vapply(shared,
                      function(x) prod(stats::dnorm(Y, mu + A * x, sigmas[s])),
                      numeric(1))) / prod(stats::dnorm(Y, mu, sigmas[s])))
    }))

  # And the same against binomial densities. The shifted outcome model is
  # m_A(X; beta) + A xi with m_A the conditional mean, so for a binary outcome
  # the shift is added to the fitted probability and not to the log odds. The
  # coefficients here keep every shifted probability inside (0, 1), where that
  # model exists.
  betas_b <- matrix(stats::rnorm(S * 4, c(-0.5, 0.5, 0.1, -0.1), 0.05),
                    nrow = S, byrow = TRUE)
  eta_b   <- Z[A == 1, , drop = FALSE] %*% t(betas_b)
  Yb      <- stats::rbinom(n, 1, stats::plogis(drop(Z %*% betas_b[1, ])))
  probability <- function(s) stats::plogis(drop(Z %*% betas_b[s, ]))

  # Every shifted probability has to stay inside (0, 1), where the shifted
  # model exists and the densities below are defined.
  shared_b  <- shared / 2
  each_b    <- each / 2
  treated_p <- vapply(seq_len(S), function(s) range(probability(s)[A == 1]),
                      numeric(2))
  expect_gt(min(treated_p) + min(shared_b, each_b), 0)
  expect_lt(max(treated_p) + max(shared_b, each_b), 1)

  expect_equal(
    sensitivity_log_weights(eta_b, Yb[A == 1], shared_b, "logit", NA,
                            "common")$log_weights,
    by_draw(function(s) {
      mu <- probability(s)
      log(mean(vapply(shared_b, function(x) {
        prod(stats::dbinom(Yb, 1, mu + A * x))
      }, numeric(1))) / prod(stats::dbinom(Yb, 1, mu)))
    }))
  expect_equal(
    sensitivity_log_weights(eta_b, Yb[A == 1], each_b[A == 1, , drop = FALSE],
                            "logit", NA, "per-observation")$log_weights,
    by_draw(function(s) {
      mu <- probability(s)
      sum(log(vapply(seq_len(n), function(i) {
        mean(stats::dbinom(Yb[i], 1, mu[i] + A[i] * each_b[i, ])) /
          stats::dbinom(Yb[i], 1, mu[i])
      }, numeric(1))))
    }))
})


test_that("the two algorithms agree exactly at a point mass", {
  # With one value of xi for everybody, Algorithm 3 and Algorithm A.1 are the
  # same integral, so the two accumulators must return the same number however
  # far the log ratios run from zero. The per-observation one sums over
  # observations last, and at a small residual scale its individual terms fall
  # below the smallest representable exponential.
  set.seed(101)
  n <- 60
  Z <- cbind(1, stats::rbinom(n, 1, 0.5), stats::rnorm(n), stats::rnorm(n))
  A <- Z[, 2]
  Y <- drop(Z %*% c(1, 2, 0.8, 0.5)) + stats::rnorm(n)
  betas <- matrix(stats::rnorm(7 * 4, c(1, 2, 0.8, 0.5), 0.2), nrow = 7,
                  byrow = TRUE)
  eta <- Z[A == 1, , drop = FALSE] %*% t(betas)
  y   <- Y[A == 1]

  for (sigma in c(1, 0.1, 0.05, 0.02)) {
    shared <- sensitivity_log_weights(eta, y, rep(0.5, 3), "identity", sigma,
                                      "common")$log_weights
    apiece <- sensitivity_log_weights(eta, y, matrix(0.5, nrow(eta), 3),
                                      "identity", sigma,
                                      "per-observation")$log_weights
    expect_equal(apiece, shared)
    expect_true(all(is.finite(apiece)))
  }
  # At the smallest of those scales 80 of the 217 log ratios are below -745,
  # which is where an exponential taken against a reference of zero underflows.
  expect_lt(min(outcome_log_ratio(outcome_scale(eta, "identity", 0.02), y,
                                  0.5)), -745)

  # A binary outcome reaches the same place without any small scale, through
  # fitted probabilities that double precision cannot tell from zero and a
  # shift that takes the mean out of the unit interval, both of which are held
  # at the boundary and have to be held there identically in the two.
  eta_far <- matrix(c(-100, -60, 0, 2), nrow = 4, ncol = 3)
  y_far   <- c(1, 1, 0, 1)
  apiece  <- sensitivity_log_weights(eta_far, y_far, matrix(-0.5, 4, 3),
                                     "probit", NA,
                                     "per-observation")$log_weights
  expect_equal(apiece,
               sensitivity_log_weights(eta_far, y_far, rep(-0.5, 3), "probit",
                                       NA, "common")$log_weights)
  expect_true(all(is.finite(apiece)))
})


test_that("the weights stay finite where the product of ratios underflows", {
  skip_on_cran()
  example <- large_fit()

  set.seed(5)
  sens <- suppressWarnings(
    drbayes_sensitivity(example$fit, xi = function(n) stats::runif(n, 1.8, 2.2),
                        M = 5, method = "per-observation"))

  expect_true(all(is.finite(sens$log_weights)))
  expect_true(all(is.finite(sens$weights)))
  expect_true(all(sens$weights >= 0))
  expect_equal(sum(sens$weights), 1)

  # What the log scale is for: at this sample size the weight itself, the
  # product over observations of the averaged ratios, is zero for every draw,
  # so weights taken from it would be 0/0.
  expect_true(all(exp(sens$log_weights) == 0))
})


test_that("a binary outcome is shifted on the scale of the risk", {
  skip_on_cran()
  example <- binary_fit()

  set.seed(31)
  expect_equal(drbayes_sensitivity(example$fit, xi = 0, M = 4)$weights,
               rep(1 / length(example$fit$pc), length(example$fit$pc)))

  # m_A(X; beta, xi) = m_A(X; beta) + A xi is a shift of the conditional mean,
  # which for a binary outcome is the risk. So xi comes off the risk difference
  # one for one, as it does for a gaussian outcome, and not through the link:
  # added to the log odds instead this shift would move the effect by about an
  # eighth of E[xi] on this fit.
  set.seed(13)
  sens  <- drbayes_sensitivity(example$fit, xi = xi_triangular(0, 0.02),
                               M = 200)
  shift <- mean(sens$ate_gcomp) - mean(example$fit$g.comp)

  expect_lt(shift, 0)
  expect_lt(abs(shift + mean(sens$xi)), 0.35 * mean(sens$xi))
  expect_equal(sens$n_saturated[["saturated"]], 0)
  expect_true(is.na(sens$sigma))
  set.seed(13)
  expect_warning(drbayes_sensitivity(example$fit, xi = 0, M = 4, sigma = 2),
                 "logit link, so it is ignored")
})


test_that("a shift the outcome has no room for is held and reported", {
  skip_on_cran()
  example <- binary_fit()

  # A risk of 0.15 cannot be lowered by 0.5, so the shifted mean is held at
  # zero for most of the treated observations and the analysis says so rather
  # than reporting a bias it did not apply.
  set.seed(35)
  expect_warning(
    sens <- drbayes_sensitivity(example$fit, xi = xi_triangular(-0.5, 0),
                                M = 100),
    "leaves \\[0, 1\\]")
  expect_gt(sens$n_saturated[["saturated"]], 0.5 * sens$n_saturated[["of"]])
  expect_equal(sens$n_saturated[["of"]],
               sum(example$fit$particles$model_data$A == 1L) * sens$n_draws)
  expect_output(print(sens), "Shift out of range")

  # A positive shift of the same size has room in the other direction, and an
  # identity link is unbounded and never held anywhere.
  set.seed(35)
  up <- suppressWarnings(
    drbayes_sensitivity(example$fit, xi = xi_triangular(0, 0.1), M = 100))
  expect_lt(up$n_saturated[["saturated"]], 0.01 * up$n_saturated[["of"]])
  set.seed(35)
  gaussian <- suppressWarnings(
    drbayes_sensitivity(gaussian_fit()$fit, xi = xi_triangular(0, 2), M = 20))
  expect_equal(gaussian$n_saturated[["saturated"]], 0L)
})


test_that("the distributions of Section 6.3 are reported honestly", {
  skip_on_cran()

  # Paper line 338: a triangular distribution on (0, 0.5) and one on
  # (-0.5, 0). Both ask for a posterior that sits E[xi] = 1/3 away from the
  # original one, which is three posterior standard deviations of the gaussian
  # effect here and ten of the binary one, whose risk difference is 0.04, and
  # at these numbers of draws neither leaves a usable reweighted posterior.
  # The requirement is not that they work; it is that nothing is called doubly
  # robust, that whichever way the analysis fell apart is on the screen, and
  # that the direction of the shift is still the right one.
  for (example in list(gaussian_fit(), binary_fit())) {
    for (bounds in list(c(0, 0.5), c(-0.5, 0))) {
      set.seed(32)
      sens <- suppressWarnings(
        drbayes_sensitivity(example$fit,
                            xi = xi_triangular(bounds[1], bounds[2]), M = 200))

      expect_false(sens$smc$converged)
      expect_equal(sign(mean(sens$ate_gcomp) - mean(example$fit$g.comp)),
                   -sign(sum(bounds)))
      expect_true(all(is.na(
        summary(sens)$estimands["sensitivity", c("2.5%", "97.5%")])))

      printed <- utils::capture.output(print(sens))
      expect_true(any(grepl("Still doubly robust:   NO", printed,
                            fixed = TRUE)))
      expect_true(any(grepl("credible interval of the sensitivity row",
                            printed)))
    }
  }

  # On the gaussian fit both distributions spend the posterior outright, which
  # is the failure the weights themselves report.
  set.seed(32)
  expect_warning(drbayes_sensitivity(gaussian_fit()$fit,
                                     xi = xi_triangular(0, 0.5), M = 200),
                 "effective sample size")
})


test_that("a credible interval is not reported off a handful of particles", {
  reweighted <- function(ess) {
    structure(list(ate          = seq(1, 2, length.out = 3000),
                   ate_gcomp    = seq(1, 2, length.out = 3000),
                   ate_original = seq(1, 2, length.out = 3000),
                   ess = ess, n_draws = 3000,
                   xi = c(0.4, 0.5), xi_ess = c(median = 2, min = 2),
                   n_saturated = c(saturated = 0L, of = 3000L),
                   method = "common", M = 2, sigma = NA_real_,
                   smc = list(lambda = 0, B_mean = 0, tol = 1,
                              converged = TRUE, converged_coupling = TRUE,
                              ess_sensitivity = ess, ess_frac = 0.1,
                              ess_floor = 300)),
              class = c("DRBayes_sensitivity", "list"))
  }

  # Forty effective draws out of three thousand: one of them lies beyond the
  # 2.5% point, so what is reported there would be that draw and not a
  # quantile, while the median, which twenty lie either side of, is one.
  scarce <- summary(reweighted(40))$estimands
  expect_true(all(is.na(scarce["sensitivity", c("2.5%", "97.5%")])))
  expect_false(is.na(scarce["sensitivity", "50%"]))
  expect_false(anyNA(scarce["original", ]))
  expect_output(print(reweighted(40)),
                "credible interval of the sensitivity row is not")

  # Four hundred put ten beyond each end, which is the least a 95% interval
  # can be built from.
  expect_false(anyNA(summary(reweighted(400))$estimands))
  expect_true(all(is.na(
    summary(reweighted(399))$estimands["sensitivity", c("2.5%", "97.5%")])))

  # Nineteen leave the median an order statistic too.
  expect_true(is.na(summary(reweighted(19))$estimands["sensitivity", "50%"]))

  plenty <- utils::capture.output(print(reweighted(2000)))
  expect_false(any(grepl("credible interval of the sensitivity row", plenty)))
})


test_that("the integral over g() reports how many draws it really rests on", {
  skip_on_cran()
  example <- gaussian_fit()

  # The integrand of Algorithm A.1 is a likelihood over hundreds of
  # observations, far narrower than g() is, so most of the M draws contribute
  # nothing at all and the weight is a much rougher Monte Carlo average than M
  # suggests.
  set.seed(34)
  coarse <- suppressWarnings(
    drbayes_sensitivity(example$fit, xi = xi_triangular(0, 0.5), M = 100))
  expect_lt(coarse$xi_ess[["median"]], 0.3 * coarse$M)
  expect_lte(coarse$xi_ess[["min"]], coarse$xi_ess[["median"]])

  set.seed(34)
  expect_warning(
    drbayes_sensitivity(example$fit, xi = xi_triangular(0, 0.5), M = 12),
    "Monte Carlo integral over the sensitivity distribution")

  set.seed(34)
  fine <- suppressWarnings(
    drbayes_sensitivity(example$fit, xi = xi_triangular(0, 0.5), M = 1000))
  expect_gt(fine$xi_ess[["median"]], coarse$xi_ess[["median"]])

  # Per observation the integrals are of one term each, and there are as many
  # of them as there are treated observations and draws.
  set.seed(34)
  apiece <- drbayes_sensitivity(example$fit, xi = xi_triangular(0, 0.02),
                                M = 20, method = "per-observation")
  expect_true(all(apiece$xi_ess <= 20))
  expect_gt(apiece$xi_ess[["median"]], 1)
})


test_that("the residual scale is drawn rather than plugged in", {
  skip_on_cran()
  example <- gaussian_fit()

  set.seed(41)
  drawn <- drbayes_sensitivity(example$fit, xi = xi_triangular(0, 0.02),
                               M = 50)
  expect_length(drawn$sigma, drawn$n_draws)

  # sigma^2 is drawn from its full conditional given each set of coefficients,
  # so its spread is the posterior one, sigma / sqrt(2n), and not the much
  # smaller spread of the residual sum of squares across the draws alone.
  expect_equal(mean(drawn$sigma), 1, tolerance = 0.05)
  expect_equal(stats::sd(drawn$sigma), 1 / sqrt(2 * 400), tolerance = 0.25)

  # Which matters because the answer moves with it: one posterior standard
  # deviation of sigma is worth a fifth of the shift being reported.
  shift <- function(sigma) {
    set.seed(42)
    fixed <- suppressWarnings(
      drbayes_sensitivity(example$fit, xi = xi_triangular(0, 0.2), M = 200,
                          sigma = sigma))
    mean(fixed$ate) - mean(example$fit$pc)
  }
  low  <- shift(1 - 1 / sqrt(2 * 400))
  high <- shift(1 + 1 / sqrt(2 * 400))
  expect_gt(abs(low - high), 0.1 * abs(low))

  # A vector of draws from a sampler that returns its own sigma2 chain is
  # taken as it stands.
  set.seed(43)
  supplied <- drbayes_sensitivity(
    example$fit, xi = 0.02, M = 1,
    sigma = stats::runif(length(example$fit$pc), 0.9, 1.1))
  expect_length(supplied$sigma, supplied$n_draws)

  expect_error(drbayes_sensitivity(example$fit, xi = 0, sigma = c(1, 2)),
               "one value for each of the")
})


test_that("a larger residual scale weakens the reweighting", {
  skip_on_cran()
  example <- gaussian_fit()

  # The gaussian log likelihood ratio is divided by the error variance, so the
  # noisier the outcome model is taken to be, the less the same bias tells us
  # about the coefficients.
  set.seed(9)
  tight <- suppressWarnings(
    drbayes_sensitivity(example$fit, xi = 0.1, M = 1, sigma = 1))
  set.seed(9)
  loose <- drbayes_sensitivity(example$fit, xi = 0.1, M = 1, sigma = 10)

  expect_lt(tight$ess, loose$ess)
  expect_gt(loose$ess, 0.99 * loose$n_draws)
})


test_that("the weights are computed from the data the fit was made on", {
  skip_on_cran()
  example <- gaussian_fit()

  # Nothing here is looked up by name, so rebinding the name the fit was made
  # under cannot change the answer by a single bit.
  set.seed(51)
  before <- drbayes_sensitivity(example$fit, xi = xi_triangular(0, 0.05),
                                M = 100)
  rebound <- sim_data(n = 400, seed = 99)
  assign("d", rebound, envir = environment())
  set.seed(51)
  after <- drbayes_sensitivity(example$fit, xi = xi_triangular(0, 0.05),
                               M = 100)
  expect_equal(after$log_weights, before$log_weights)
  expect_equal(after$ate, before$ate)

  expect_error(drbayes_sensitivity(example$fit, xi = 0, data = rebound),
               "takes no data argument")
})


test_that("a fit without its particles says how to get them", {
  skip_on_cran()
  example <- gaussian_fit()

  # What drbayes_pc() returns without keep_particles.
  without <- example$fit
  without$particles <- NULL
  expect_error(drbayes_sensitivity(without, xi = 0), "keep_particles = TRUE")

  # And a fit that kept the coupled particles but not the draws they came from
  # or the data they were fitted to.
  partial <- example$fit
  partial$particles <- list(outcome = example$fit$particles$outcome)
  expect_error(drbayes_sensitivity(partial, xi = 0), "keep_particles = TRUE")
})


test_that("the sensitivity distribution and its settings are validated", {
  skip_on_cran()
  example <- gaussian_fit()
  args <- list(example$fit)

  expect_error(do.call(drbayes_sensitivity, c(args, list(xi = "0.2"))),
               "must be either a function")
  expect_error(do.call(drbayes_sensitivity, c(args, list(xi = numeric(0)))),
               "must be either a function")
  expect_error(do.call(drbayes_sensitivity,
                       c(args, list(xi = function(n) 0.2))),
               "returned 1 value")
  expect_error(do.call(drbayes_sensitivity, c(args, list(xi = c(0.1, NA)))),
               "must all be finite")
  expect_error(do.call(drbayes_sensitivity, c(args, list(xi = 0, M = 0))),
               "M must be a single positive integer")
  expect_error(do.call(drbayes_sensitivity, c(args, list(xi = 0, sigma = -1))),
               "sigma must be positive and finite")
  expect_error(do.call(drbayes_sensitivity, c(args, list(xi = 0, Method = 1))),
               "does not use the argument")
  expect_error(drbayes_sensitivity(example$fit$pc, xi = 0),
               "must be a fit of class")
})


test_that("the shifted log likelihood ratio agrees with the densities", {
  set.seed(12)
  eta   <- matrix(stats::rnorm(12, sd = 0.3), nrow = 4)
  y     <- c(1, 0, 1, 0)
  shift <- c(0.15, -0.1, 0.08, 0.05)

  # m_A(X; beta) + A xi is the conditional mean, so for a binary outcome the
  # shift moves the fitted probability. This is the arithmetic of that, against
  # the density R itself provides for the model.
  ratio_of <- function(link, inverse_link) {
    got  <- outcome_log_ratio(outcome_scale(eta, link, NA), y, shift)
    want <- stats::dbinom(y, 1, inverse_link(eta) + shift, log = TRUE) -
      stats::dbinom(y, 1, inverse_link(eta), log = TRUE)
    expect_false(anyNA(want))
    expect_equal(got, matrix(want, nrow = nrow(eta)))
  }
  ratio_of("logit", stats::plogis)
  ratio_of("probit", stats::pnorm)

  # For an identity link the mean is the linear predictor and nothing changes.
  expect_equal(
    outcome_log_ratio(outcome_scale(eta, "identity", 1.7),
                      c(2.5, -1, 0.4, 3), shift),
    matrix(stats::dnorm(c(2.5, -1, 0.4, 3), eta + shift, 1.7, log = TRUE) -
             stats::dnorm(c(2.5, -1, 0.4, 3), eta, 1.7, log = TRUE),
           nrow = nrow(eta)))

  # A shift that takes the mean out of [0, 1] leaves no distribution to
  # evaluate, so the mean is held at the boundary: the ratio stays finite where
  # the density route returns NaN, and an observation contributes at most
  # log(1 / eps) = 36 to a log weight.
  over <- outcome_log_ratio(outcome_scale(matrix(c(-1, 2), nrow = 2), "logit",
                                          NA), c(0, 0), 0.9)
  expect_true(all(is.finite(over)))
  expect_gte(min(over), -36.1)
  expect_true(anyNA(suppressWarnings(
    stats::dbinom(0, 1, stats::plogis(2) + 0.9))))

  # Far out in the tail the density route loses the answer entirely, because
  # the probability underflows before its logarithm is taken.
  far <- matrix(c(-40, 40), nrow = 2)
  expect_true(all(is.finite(
    outcome_log_ratio(outcome_scale(far, "probit", NA), c(1, 0), 0.5))))
  expect_false(all(is.finite(
    stats::dbinom(c(1, 0), 1, stats::pnorm(far), log = TRUE))))

  # And a fitted probability that double precision cannot tell from zero or one
  # is held at the boundary in the same way, so that a sensitivity parameter of
  # zero still leaves the likelihood ratio exactly zero.
  expect_equal(outcome_log_ratio(outcome_scale(far, "probit", NA), c(1, 0), 0),
               matrix(0, nrow = 2, ncol = 1))

  expect_error(outcome_scale(eta, "cauchit", NA),
               "no likelihood the sensitivity analysis knows how to shift")
})


test_that("xi_triangular reproduces the distributions of Section 6.3", {
  positive <- xi_triangular(0, 0.5)
  negative <- xi_triangular(-0.5, 0)

  set.seed(10)
  up   <- positive(20000)
  down <- negative(20000)

  expect_true(all(up > 0 & up < 0.5))
  expect_true(all(down > -0.5 & down < 0))
  # A triangular density on (a, b) peaking at b has mean (a + 2b)/3, so the
  # peak defaults to the endpoint further from zero, as in the paper.
  expect_equal(mean(up), 1 / 3, tolerance = 0.02)
  expect_equal(mean(down), -1 / 3, tolerance = 0.02)
  expect_equal(mean(xi_triangular(-0.5, 0.5, mode = 0)(20000)), 0,
               tolerance = 0.01)

  expect_error(xi_triangular(0.5, 0), "lower must be smaller than upper")
  expect_error(xi_triangular(0, 0.5, mode = 1), "mode must be a single number")
  expect_error(xi_triangular(0, NA), "single finite number")
})


test_that("the printed output leads with the diagnostics", {
  skip_on_cran()
  example <- gaussian_fit()

  set.seed(61)
  sens <- drbayes_sensitivity(example$fit, xi = xi_triangular(0, 0.02),
                              M = 50)

  expect_s3_class(sens, "DRBayes_sensitivity")
  expect_output(print(sens), "Effective sample size")
  expect_output(print(sens), "Algorithm A.1")
  expect_output(print(sens), "integral resolved by a median")
  expect_output(print(sens), "Residual scale sigma")
  expect_output(print(sens), "Moment condition")
  expect_output(print(sens), "Still doubly robust")
  expect_output(print(summary(sens)), "Still doubly robust")

  summarised <- summary(sens)
  expect_s3_class(summarised, "summary.DRBayes_sensitivity")
  expect_equal(rownames(summarised$estimands),
               c("sensitivity", "reweighted", "original"))
  expect_equal(summarised$estimands["original", "mean"],
               mean(example$fit$pc))
  expect_output(print(summarised), "Sensitivity parameter xi")

  # The number of decimal places comes from the posterior standard deviation,
  # so a wide interval around a large effect is not rounded away.
  expect_output(print(sens), "2\\.0")
})

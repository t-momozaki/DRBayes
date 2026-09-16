# Draws now arrive as iterations x chains x parameters. Most of these tests care
# about the marginal posterior, so pool the chains after discarding a burn-in.
pooled_draws <- function(draws, burn = 0L) {
  draws <- suppressMessages(draws)
  keep  <- seq(burn + 1L, dim(draws)[1])
  matrix(draws[keep, , , drop = FALSE], ncol = dim(draws)[3],
         dimnames = list(NULL, dimnames(draws)[[3]]))
}

# Regression tests for the Gibbs samplers. Three properties are easy to lose
# and hard to notice: which coefficients the horseshoe exempts from shrinkage,
# that the amount of shrinkage does not depend on the units of the data, and
# that the default prior does not inflate the error variance and with it every
# credible interval.

sparse_data <- function(n = 300, p = 40, effect = 0.5, sigma = 3, seed = 7) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("V", seq_len(p))))
  d <- as.data.frame(X)
  d$A <- rbinom(n, 1, plogis(0.6 * d$V1 - 0.6 * d$V2 + 0.4 * d$V3))
  d$Y <- 1 + effect * d$A + 2 * d$V1 - 1.5 * d$V2 + d$V3 + rnorm(n, sd = sigma)
  d
}

# Coefficients and a residual variance far from one, as in section 5.1 of the
# paper. This is where a coefficient prior fixed on an absolute scale bites.
wide_scale_data <- function(n = 200, sigma = 40, seed = 17) {
  set.seed(seed)
  X <- cbind(A = rbinom(n, 1, 0.5), V1 = rnorm(n), V2 = rnorm(n))
  list(X = X,
       Y = 100 + 110 * X[, "A"] + 27 * X[, "V1"] + rnorm(n, sd = sigma))
}

# The outcome model of section 5.1 has no interaction with the treatment, so
# the G-formula estimate of the ATE is the treatment coefficient itself and
# Table 1 is a statement about the posterior this sampler returns.
paper_data <- function(seed) {
  d <- generate_dataset(nn = 500, pp = 0, seed = seed)
  list(d = d,
       correct = as.matrix(d[, c("AA", "WW1", "WW2", "WW3", "WW4")]),
       wrong = as.matrix(d[, c("AA", "WW1")]))
}

# The posterior of theta ~ N(0, V), sigma^2 ~ InvGamma(a, b) has no closed
# form, but it can be integrated in one dimension without reusing anything the
# sampler does. Marginalising theta out of the likelihood leaves
# Y | sigma^2 ~ N(0, sigma^2 I + Z V Z'), whose density is written below
# through the identities det(A) = sigma^2n det(V^-1 + Z'Z/sigma^2) det(V) and
# A^-1 = (I - Z (sigma^2 V^-1 + Z'Z)^-1 Z') / sigma^2. Averaging the normal
# full conditional of theta over that density gives its exact moments.
lm_reference <- function(Y, Z, precision, a = 1, b = 1, points = 1200) {
  nn  <- length(Y)
  dd  <- ncol(Z)
  ZtZ <- crossprod(Z)
  ZtY <- drop(crossprod(Z, Y))
  YtY <- drop(crossprod(Y))

  # A log-spaced grid ten standard deviations either side of the residual mean
  # square, which is where all but a vanishing part of the mass sits.
  rss   <- sum(lm.fit(Z, Y)$residuals^2)
  scale <- rss / max(nn - dd, 1)
  width <- 10 * sqrt(2 / nn)
  s2    <- exp(seq(log(scale) - width, log(scale) + width, length.out = points))

  log_w <- numeric(points)
  mean_k <- matrix(0, points, dd)
  var_k  <- matrix(0, points, dd)
  for (k in seq_len(points)) {
    Qn <- ZtZ + s2[k] * precision
    Un <- chol(Qn)
    mk <- drop(backsolve(Un, backsolve(Un, ZtY, transpose = TRUE)))
    mean_k[k, ] <- mk
    var_k[k, ]  <- s2[k] * diag(chol2inv(Un))
    # Prior, marginal likelihood and the log Jacobian of the log grid.
    log_w[k] <- (-a - 1) * log(s2[k]) - b / s2[k] -
      0.5 * (nn * log(s2[k]) + 2 * sum(log(diag(Un))) - dd * log(s2[k])) -
      0.5 * (YtY - sum(mk * ZtY)) / s2[k] + log(s2[k])
  }
  w <- exp(log_w - max(log_w))
  w <- w / sum(w)

  mean_theta <- drop(crossprod(w, mean_k))
  second     <- drop(crossprod(w, var_k + mean_k^2))
  list(mean = mean_theta, sd = sqrt(second - mean_theta^2),
       sigma2 = sum(w * s2))
}


test_that("bayes_lm targets the posterior of the model it documents", {
  skip_on_cran()
  d <- wide_scale_data()
  Z <- cbind(1, d$X)

  # The prior is informative on purpose: with an effectively flat one every
  # sampler agrees with least squares and the test would prove nothing.
  exact <- lm_reference(d$Y, Z, precision = diag(0.01, ncol(Z)))

  set.seed(21)
  draws <- pooled_draws(bayes_lm(d$Y, d$X, mc = 20000, chains = 4L,
                                 theta_prior = 0.01),
                        burn = 2000)

  # Distances are measured in posterior standard deviations, which is the only
  # scale on which "close" means anything here.
  expect_lt(max(abs(colMeans(draws) - exact$mean) / exact$sd), 0.1)
  expect_lt(max(abs(apply(draws, 2, sd) / exact$sd - 1)), 0.05)
})


test_that("bayes_lm reproduces the interval length of the paper's Table 1", {
  skip_on_cran()
  # Table 1 reports an average credible interval length of 0.397 for the ATE
  # at n = 500 with the outcome model correctly specified. A prior that
  # inflates the error variance leaves the point estimate alone and widens
  # every interval, so this is the number that catches it.
  reps    <- 5L
  lengths <- numeric(reps)
  ratios  <- numeric(reps)
  for (j in seq_len(reps)) {
    p <- paper_data(seed = j)
    set.seed(100 + j)
    draws <- pooled_draws(bayes_lm(p$d$YY, p$correct, mc = 4000, chains = 2L),
                          burn = 1000)
    lengths[j] <- diff(quantile(draws[, "AA"], c(0.025, 0.975)))
    fit <- lm(YY ~ AA + WW1 + WW2 + WW3 + WW4, data = p$d)
    ratios[j] <- max(apply(draws, 2, sd) /
                       summary(fit)$coefficients[, "Std. Error"])
  }

  expect_lt(abs(mean(lengths) - 0.397), 0.02)
  # The same inflation seen a second way. A coefficient prior conditioned on
  # sigma^2 put these ratios at 1.21, because the prior quadratic form was
  # then part of the error variance.
  expect_lt(max(ratios), 1.05)
})


test_that("the default prior barely moves the paper's misspecified fit", {
  skip_on_cran()
  # The other half of Table 1: with only WW1 in the outcome model the residual
  # variance is about 550 and the treatment coefficient about 110. A fixed
  # prior precision of 1/100 pulls the posterior mean about 2.5 below the
  # least squares estimate there, which is most of the absolute bias of 2.535
  # the paper reports for the G-formula, so the default has to leave it alone.
  for (j in 1:3) {
    p <- paper_data(seed = j)
    set.seed(200 + j)
    draws <- pooled_draws(bayes_lm(p$d$YY, p$wrong, mc = 6000, chains = 2L),
                          burn = 1000)
    ols <- unname(coef(lm(YY ~ AA + WW1, data = p$d)))
    expect_lt(max(abs(colMeans(draws) - ols)), 0.1)
  }
})


test_that("the default prior shrinks the same however the data are scaled", {
  skip_on_cran()
  d <- wide_scale_data()

  set.seed(31)
  unit <- pooled_draws(bayes_lm(d$Y, d$X, mc = 6000, chains = 2L), burn = 1000)
  set.seed(31)
  # The same outcome read in units a thousand times smaller, so every number
  # is a thousand times larger. The default is read off sd(Y) and off the
  # spread of each column, so it follows.
  large <- pooled_draws(bayes_lm(1000 * d$Y, d$X, mc = 6000, chains = 2L),
                        burn = 1000)

  rescaled <- colMeans(large) / 1000
  expect_lt(max(abs(rescaled - colMeans(unit))) / max(apply(unit, 2, sd)), 0.05)

  # And the same when a column of X changes units rather than Y.
  X_small <- d$X
  X_small[, "V1"] <- X_small[, "V1"] / 1000
  set.seed(31)
  narrow <- pooled_draws(bayes_lm(d$Y, X_small, mc = 6000, chains = 2L),
                         burn = 1000)
  expect_lt(abs(colMeans(narrow)[["V1"]] / 1000 - colMeans(unit)[["V1"]]) /
              apply(unit, 2, sd)[["V1"]], 0.05)

  # And the posterior mean is the least squares estimate, up to the ridge
  # penalty the prior contributes and Monte Carlo error.
  ols <- unname(coef(lm(d$Y ~ d$X)))
  expect_lt(max(abs(colMeans(unit) - ols) / apply(unit, 2, sd)), 0.1)
})


test_that("the horseshoe shrinks the same however Y is scaled", {
  skip_on_cran()
  d <- sparse_data(n = 200, p = 12, effect = 2, sigma = 1)
  X <- as.matrix(d[, c("A", paste0("V", 1:12))])

  set.seed(41)
  unit <- pooled_draws(bayes_lm_hs(d$Y, X, mc = 8000, chains = 2L,
                                   unshrunk = 1L), burn = 2000)
  set.seed(41)
  large <- pooled_draws(bayes_lm_hs(1000 * d$Y, X, mc = 8000, chains = 2L,
                                    unshrunk = 1L), burn = 2000)

  # The horseshoe scales carry a factor of sigma and the two fixed priors are
  # read off the data, so the whole posterior follows the units of Y and the
  # shrinkage stays where it was.
  rescaled <- colMeans(large) / 1000
  expect_lt(max(abs(rescaled - colMeans(unit))) / max(apply(unit, 2, sd)), 0.1)
})


test_that("the horseshoe does not inflate the error variance either", {
  skip_on_cran()
  # Same test as for bayes_lm, on the paper's correctly specified model. The
  # intercept of 100 and the treatment effect of 110 used to be given prior
  # variances in units of sigma^2, which added 221 to a residual sum of
  # squares of 499 and widened every interval by a fifth.
  p <- paper_data(seed = 1)
  set.seed(61)
  draws <- pooled_draws(bayes_lm_hs(p$d$YY, p$correct, mc = 8000, chains = 2L,
                                    unshrunk = 1L), burn = 2000)
  fit <- lm(YY ~ AA + WW1 + WW2 + WW3 + WW4, data = p$d)

  ratio <- apply(draws, 2, sd) / summary(fit)$coefficients[, "Std. Error"]
  expect_lt(max(ratio), 1.05)
  expect_lt(abs(mean(draws[, "AA"]) - coef(fit)[["AA"]]), 0.05)
})


test_that("unshrunk exempts exactly the coefficients it names", {
  skip_on_cran()
  d <- sparse_data(n = 300, p = 30, effect = 0.5, sigma = 3)
  X <- as.matrix(d[, c("A", paste0("V", 1:30))])

  set.seed(1)
  exempt <- pooled_draws(bayes_lm_hs(d$Y, X, mc = 4000, unshrunk = 1L),
                         burn = 1000)
  set.seed(1)
  shrunk <- pooled_draws(bayes_lm_hs(d$Y, X, mc = 4000, unshrunk = integer(0)),
                         burn = 1000)

  # A weakly identified treatment effect survives when it is exempted and is
  # pulled towards zero when it is not.
  expect_gt(abs(mean(exempt[, "A"])), abs(mean(shrunk[, "A"])))

  # A genuinely null covariate is shrunk either way.
  expect_lt(abs(mean(exempt[, "V20"])), 0.5)
})


test_that("the shrinkage set does not depend on column position", {
  skip_on_cran()
  d <- sparse_data(n = 250, p = 20, effect = 0.5, sigma = 3)
  first <- as.matrix(d[, c("A", paste0("V", 1:20))])
  mid   <- as.matrix(d[, c(paste0("V", 1:10), "A", paste0("V", 11:20))])

  set.seed(3)
  a <- pooled_draws(bayes_lm_hs(d$Y, first, mc = 3000, unshrunk = 1L),
                    burn = 800)
  set.seed(3)
  b <- pooled_draws(bayes_lm_hs(d$Y, mid, mc = 3000, unshrunk = 11L),
                    burn = 800)

  expect_lt(abs(mean(a[, "A"]) - mean(b[, "A"])), 0.25)
})


test_that("the binary horseshoe samplers can exempt the treatment too", {
  skip_on_cran()
  d <- sparse_data(n = 300, p = 20)
  d$Yb <- rbinom(nrow(d), 1, plogis(-0.5 + 0.5 * d$A + 0.8 * d$V1))
  X <- as.matrix(d[, c("A", paste0("V", 1:20))])

  expect_true("unshrunk" %in% names(formals(bayes_logit_hs)))
  expect_true("unshrunk" %in% names(formals(bayes_probit_hs)))

  set.seed(2)
  exempt <- pooled_draws(bayes_logit_hs(d$Yb, X, mc = 3000, unshrunk = 1L),
                         burn = 800)
  set.seed(2)
  shrunk <- pooled_draws(bayes_logit_hs(d$Yb, X, mc = 3000,
                                        unshrunk = integer(0)), burn = 800)

  expect_gt(abs(mean(exempt[, "A"])), abs(mean(shrunk[, "A"])))
})


test_that("the horseshoe prior reaches the diagonal it is meant to", {
  skip_on_cran()
  # The coefficient precision differs from the cross product of the design only
  # along its diagonal, and the horseshoe samplers write the prior straight into
  # those positions. One covariate is the smallest design that can be got wrong:
  # the intercept and the slope then sit at opposite corners of a two by two
  # matrix, so a near-degenerate intercept prior must hold the intercept at zero
  # and leave the slope entirely to the data. A prior that landed one position
  # out would reach an off-diagonal entry instead, which chol() does not
  # complain about because it reads only one triangle.
  set.seed(11)
  n <- 400
  X <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "V1"))
  Y <- 3 + 2 * drop(X) + rnorm(n, sd = 0.5)
  Yb <- rbinom(n, 1, plogis(1 + 2 * drop(X)))

  fits <- list(bayes_lm_hs(Y, X, mc = 3000, chains = 2L, beta0_prior = 1e6),
               bayes_logit_hs(Yb, X, mc = 3000, chains = 2L,
                              beta0_prior = 1e6),
               bayes_probit_hs(Yb, X, mc = 3000, chains = 2L,
                               beta0_prior = 1e6))

  for (draws in fits) {
    pooled <- pooled_draws(draws, burn = 1000)
    # A prior standard deviation of a thousandth leaves the intercept nowhere
    # else to be.
    expect_lt(abs(mean(pooled[, "(Intercept)"])), 0.05)
    # The slope carries a signal several times its own posterior spread, and
    # nothing in this fit is entitled to shrink it away.
    expect_gt(mean(pooled[, "V1"]), 0.5)
  }
})


test_that("init decides where a chain starts", {
  # The scan draws the coefficients from a full conditional that does not
  # depend on their previous value, so a starting value can only reach the
  # sampler through the variance components. If it does not, every chain
  # begins in the same state and R-hat has nothing to compare.
  d <- sparse_data(n = 120, p = 3)
  X <- as.matrix(d[, c("A", "V1", "V2")])
  same <- list(rep(0, 4), rep(0, 4))
  wild <- list(rep(0, 4), c(1e3, -1e3, 1e3, -1e3))

  for (sampler in list(bayes_lm, bayes_lm_hs)) {
    set.seed(9)
    a <- suppressMessages(sampler(d$Y, X, mc = 40, chains = 2L, init = same))
    set.seed(9)
    b <- suppressMessages(sampler(d$Y, X, mc = 40, chains = 2L, init = wild))

    # Chain one is given the same starting value in both fits and reproduces
    # itself exactly; chain two is not, and must not.
    expect_equal(a[, 1, ], b[, 1, ])
    expect_false(isTRUE(all.equal(a[, 2, ], b[, 2, ])))

    # A chain started far from the data starts with a large residual sum of
    # squares, so it begins at a large error variance and its first draws are
    # visibly more dispersed than a chain started at the origin.
    expect_gt(sd(b[1:5, 2, "A"]), 5 * sd(a[1:5, 2, "A"]))
  }
})


test_that("drbayes_pc exempts the treatment and its interactions", {
  skip_on_cran()
  d <- sparse_data(n = 250, p = 10)

  md <- prepare_model_data(Y ~ A + V1 + A:V1, A ~ V1 + V2, d)
  expect_identical(colnames(md$X.lm)[md$treatment_columns], c("A", "A:V1"))

  # The treatment need not be written first.
  md2 <- prepare_model_data(Y ~ V1 + A + V2, A ~ V1, d)
  expect_identical(colnames(md2$X.lm)[md2$treatment_columns], "A")
})


test_that("column names of X survive into the draws", {
  d <- sparse_data(n = 120, p = 3)
  X <- as.matrix(d[, c("A", "V1", "V2")])
  draws <- suppressMessages(bayes_lm(d$Y, X, mc = 50))
  expect_length(dim(draws), 3L)
  expect_identical(dimnames(draws)[[3]], c("(Intercept)", "A", "V1", "V2"))
  expect_identical(typeof(draws), "double")
})


test_that("an intercept-only design is allowed", {
  set.seed(4)
  n <- 150
  A <- rbinom(n, 1, 0.5)
  draws <- suppressMessages(bayes_logit(A, matrix(numeric(0), n, 0), mc = 100,
                                    chains = 2L))
  expect_identical(dim(draws), c(100L, 2L, 1L))
  expect_identical(dimnames(draws)[[3]], "(Intercept)")

  Y <- rnorm(n, mean = 5)
  lm_draws <- suppressMessages(bayes_lm(Y, matrix(numeric(0), n, 0), mc = 500,
                                        chains = 2L))
  expect_identical(dim(lm_draws), c(500L, 2L, 1L))
  expect_lt(abs(mean(lm_draws) - mean(Y)), 0.1)
})


test_that("the data-scaled default survives data that carry no scale", {
  set.seed(12)
  n <- 60
  X <- cbind(V1 = rnorm(n), constant = rep(3, n), zero = numeric(n))
  Y <- 1 + 2 * X[, "V1"] + rnorm(n)

  # A constant column is collinear with the intercept and a zero column
  # carries nothing at all. The default has to stay proper on both, or the
  # coefficient precision matrix cannot be factorised.
  expect_true(all(is.finite(bayes_lm(Y, X, mc = 300, chains = 2L))))
  expect_true(all(is.finite(bayes_lm_hs(Y, X, mc = 300, chains = 2L))))

  # A constant outcome leaves sd(Y) zero, and fewer rows than columns leaves
  # the least squares fit that places the starting values unusable.
  expect_true(all(is.finite(bayes_lm(rep(4, n), X, mc = 200, chains = 2L))))
  expect_true(all(is.finite(bayes_lm(rnorm(5), matrix(rnorm(100), 5, 20),
                                     mc = 200, chains = 2L))))

  # An outcome enormous in size but barely varying: the prior belongs on the
  # spread of Y, not on its location, and the intercept must not be shrunk
  # towards zero. The posterior standard deviation of the slope is about
  # 1/sqrt(n) here whatever the location of Y.
  set.seed(13)
  far <- pooled_draws(bayes_lm(1e8 + rnorm(n), X[, "V1", drop = FALSE],
                               mc = 3000, chains = 2L), burn = 500)
  expect_gt(abs(mean(far[, "(Intercept)"])), 9e7)
  expect_lt(sd(far[, "V1"]), 1)
  expect_gt(sd(far[, "V1"]), 0.01)
})


test_that("a default precision that cannot be represented stops the fit", {
  set.seed(71)
  n <- 60
  X <- cbind(a = rnorm(n))

  # The default intercept prior is 100 (|mean(Y)| + sd(Y)) wide, and its
  # square leaves the range of a double here, so the precision underflowed to
  # exactly zero: an improper prior whose precision matrix has no Cholesky
  # factor. What came out was "the leading minor of order 1 is not positive",
  # which names neither the coefficient nor the reason.
  huge <- tryCatch(bayes_lm(1e153 + rnorm(n), X, mc = 10, chains = 1L),
                   error = conditionMessage)
  expect_match(huge, "default prior for the intercept")
  expect_match(huge, "theta_prior")

  # The other end, and the worse of the two. A column whose spread is
  # astronomically wider than the spread of Y sends its precision to infinity,
  # and an infinite precision raises nothing at all: the sampler ran to
  # completion and returned that coefficient as exactly zero in every draw.
  set.seed(72)
  wide <- tryCatch(bayes_lm(1e-100 * rnorm(n), cbind(a = 1e100 * rnorm(n)),
                            mc = 10, chains = 1L),
                   error = conditionMessage)
  expect_match(wide, "default prior for column a of X")
  expect_false(grepl("collinear", wide))

  # An unnamed design has to say which column too, and say the right one: the
  # first column here is on a perfectly ordinary scale.
  set.seed(73)
  unnamed <- tryCatch(bayes_lm(1e-100 * rnorm(n),
                               cbind(rnorm(n), 1e100 * rnorm(n)),
                               mc = 10, chains = 1L),
                      error = conditionMessage)
  expect_match(unnamed, "default prior for column 2 of X")

  # The message offers theta_prior as the way past the default, so naming it
  # has to get the same call through. The fit means nothing on data this
  # extreme; the point is that the choice of prior is then the user's.
  set.seed(71)
  expect_true(all(is.finite(
    bayes_lm(1e153 + rnorm(n), X, mc = 20, chains = 1L,
             theta_prior = 1e-300))))

  # bayes_lm_hs reads the same scales for every column, exempt or not, so it
  # stops on the same data. A default that silently went improper there would
  # be harder to see, since the horseshoe columns carry no fixed precision of
  # their own to fail on.
  set.seed(74)
  shared <- tryCatch(bayes_lm_hs(1e153 + rnorm(n), cbind(a = rnorm(n),
                                                         b = rnorm(n)),
                                 mc = 10, chains = 1L),
                     error = conditionMessage)
  expect_match(shared, "default prior for the intercept")
})


test_that("a collinear design far from the scale of Y still fits", {
  # The reference fit that keeps the residual sum of squares accurate was
  # placed by solve(), which refuses any matrix whose reciprocal condition
  # number is below the machine epsilon. The data-scaled default is a far
  # weaker ridge than the fixed one it replaced, so a duplicated column
  # measured in units 1e8 times finer than Y took it below that threshold and
  # the fit died with "system is computationally singular" on data least
  # squares handles without complaint.
  set.seed(74)
  n <- 60
  V <- 1e-8 * rnorm(n)
  X <- cbind(V1 = V, V1_duplicate = V)
  Y <- 1e8 + rnorm(n)

  set.seed(75)
  draws <- pooled_draws(bayes_lm(Y, X, mc = 4000, chains = 2L), burn = 1000)
  ols <- lm(Y ~ V)

  # Only the intercept and the sum of the two coefficients are identified, and
  # both have to agree with least squares on the single copy of the column.
  total <- draws[, "V1"] + draws[, "V1_duplicate"]
  expect_lt(abs(mean(draws[, "(Intercept)"]) - coef(ols)[[1]]) /
              sd(draws[, "(Intercept)"]), 0.15)
  expect_lt(abs(mean(total) - coef(ols)[[2]]) / sd(total), 0.15)

  # The spread is what fails if the residual sum of squares has lost its
  # digits to cancellation at 1e8: the intercept would then carry the standard
  # error of a wildly wrong error variance rather than of the data.
  expect_lt(abs(sd(draws[, "(Intercept)"]) /
                  summary(ols)$coefficients[1, 2] - 1), 0.1)
})


test_that("a duplicated column is fitted when the outcome varies hugely", {
  # Same design, with the size of Y in its spread rather than its location.
  # The ridge is then weak enough that the Cholesky factor does not exist
  # either, and the reference falls back on least squares through a pivoted
  # QR, which answers a rank deficient design by aliasing rather than failing.
  set.seed(76)
  n <- 60
  V <- rnorm(n)
  X <- cbind(V1 = V, V1_duplicate = V)
  Y <- 1e7 * rnorm(n)

  set.seed(77)
  draws <- pooled_draws(bayes_lm(Y, X, mc = 4000, chains = 2L), burn = 1000)
  ols <- lm(Y ~ V)

  total <- draws[, "V1"] + draws[, "V1_duplicate"]
  expect_lt(abs(mean(draws[, "(Intercept)"]) - coef(ols)[[1]]) /
              sd(draws[, "(Intercept)"]), 0.15)
  expect_lt(abs(mean(total) - coef(ols)[[2]]) / sd(total), 0.15)
  expect_lt(abs(sd(draws[, "(Intercept)"]) /
                  summary(ols)$coefficients[1, 2] - 1), 0.1)
})


test_that("a column too large to square is named, not called collinear", {
  set.seed(78)
  n <- 60
  X <- cbind(ordinary = rnorm(n), enormous = 1e155 * rnorm(n))
  Y <- rnorm(n)

  # crossprod() returns Inf for the second column, everything after it sees
  # NaN, and the Cholesky factorisation that finally complains reported the
  # design as collinear or a shrinkage parameter as underflowed. Neither is
  # the cause and neither points at the fix.
  for (sampler in list(bayes_lm, bayes_lm_hs)) {
    raised <- tryCatch(sampler(Y, X, mc = 10, chains = 1L),
                       error = conditionMessage)
    expect_match(raised, "enormous")
    expect_match(raised, "not finite")
    expect_false(grepl("collinear", raised))
  }

  # And the fix the message gives has to work.
  X[, "enormous"] <- X[, "enormous"] / 1e155
  expect_true(all(is.finite(bayes_lm(Y, X, mc = 50, chains = 1L))))
})


test_that("an improper error variance prior is refused", {
  set.seed(79)
  n <- 80
  X <- cbind(V1 = rnorm(n), V2 = rnorm(n))
  Y <- rnorm(n)

  # extraDistr::rinvgamma() answers a non-positive shape or scale with NaN
  # rather than an error, so an unchecked sigma_prior would fill every chain
  # with NaN long after the call itself had looked fine.
  for (sampler in list(bayes_lm, bayes_lm_hs)) {
    expect_error(sampler(Y, X, mc = 10, sigma_prior = c(0, 1)),
                 "sigma_prior must be two positive numbers")
    expect_error(sampler(Y, X, mc = 10, sigma_prior = c(1, -1)),
                 "sigma_prior must be two positive numbers")
    expect_error(sampler(Y, X, mc = 10, sigma_prior = 1),
                 "sigma_prior must be two positive numbers")
  }
})


test_that("non-finite input is refused rather than hanging the sampler", {
  set.seed(5)
  n <- 100
  A <- rbinom(n, 1, 0.5)
  X <- cbind(V1 = rnorm(n), V2 = rnorm(n))

  # pgdraw's rejection sampler never returns on a NaN, so the logistic samplers
  # used to hang forever instead of failing.
  X_na <- X; X_na[3, 1] <- NA
  expect_error(bayes_logit(A, X_na, mc = 10), "X must be finite")
  X_inf <- X; X_inf[5, 2] <- Inf
  expect_error(bayes_logit(A, X_inf, mc = 10), "X must be finite")
  expect_error(bayes_logit_hs(A, X_inf, mc = 10), "X must be finite")

  Y <- rnorm(n); Y[2] <- NaN
  expect_error(bayes_lm(Y, X, mc = 10), "Y must be finite")
})


test_that("an improper prior precision is refused", {
  set.seed(6)
  n <- 120
  A <- rbinom(n, 1, 0.5)
  X <- cbind(V1 = rnorm(n), V2 = rnorm(n))

  # A negative precision made bayes_probit diverge to 1e23 while reporting
  # success.
  expect_error(bayes_probit(A, X, mc = 10, theta_prior = -10),
               "must be positive")
  expect_error(bayes_lm(rnorm(n), X, mc = 10, theta_prior = -1),
               "must be positive")

  # A 1x1 matrix used to fall into the scalar branch and die in diag().
  expect_error(bayes_lm(rnorm(n), X, mc = 10, theta_prior = matrix(0.01, 1, 1)),
               "must be a 3 by 3 matrix")

  asymmetric <- diag(0.01, 3); asymmetric[1, 2] <- 3
  expect_error(bayes_lm(rnorm(n), X, mc = 10, theta_prior = asymmetric),
               "must be symmetric")

  indefinite <- diag(c(1, 1, -50))
  expect_error(bayes_lm(rnorm(n), X, mc = 10, theta_prior = indefinite),
               "positive definite")

  expect_error(bayes_lm_hs(rnorm(n), X, mc = 10, beta0_prior = -1),
               "must be a single positive number")
  expect_error(bayes_lm_hs(rnorm(n), X, mc = 10, unshrunk = 1L,
                           unshrunk_prior = 0),
               "must be a single positive number")

  # tau_prior <= 0 drove tau2 to exactly zero and killed the chain thousands of
  # iterations later, after the default mc had already looked fine.
  expect_error(bayes_lm_hs(rnorm(n), X, mc = 10, tau_prior = 0),
               "tau_prior must be a single positive number")
  expect_error(bayes_lm_hs(rnorm(n), X, mc = 10, tau_prior = -1),
               "tau_prior must be a single positive number")
})


test_that("unshrunk is validated", {
  set.seed(8)
  n <- 100
  X <- cbind(V1 = rnorm(n), V2 = rnorm(n))
  Y <- rnorm(n)
  expect_error(bayes_lm_hs(Y, X, mc = 10, unshrunk = 5L), "outside X")
  expect_error(bayes_lm_hs(Y, X, mc = 10, unshrunk = 1.5), "whole numbers")
  expect_error(bayes_lm_hs(Y, X, mc = 10, unshrunk = "V1"), "column indices")
  expect_error(bayes_lm_hs(Y, X, mc = 10, unshrunk = 1:2), "nothing for the")
})


test_that("the installed Stan models are the ones in inst/stan", {
  # bayes_stan() compiles the copy inside the installed package, which
  # system.file() finds, and not the one edited here. An edit that has not
  # been installed leaves the two backends fitting different models in
  # silence: it once left bayes_stan("lm") reporting a treatment coefficient
  # 14 away from the one bayes_lm reported, which no test noticed.
  source_dir <- test_path("..", "..", "inst", "stan")
  skip_if_not(dir.exists(source_dir),
              "installed package only, so there is no source copy to compare")

  stale <- character(0)
  for (file in c("bayes_lm.stan", "bayes_lm_hs.stan")) {
    installed <- system.file("stan", file, package = "DRBayes")
    if (!nzchar(installed) ||
        !identical(readLines(installed, warn = FALSE),
                   readLines(file.path(source_dir, file), warn = FALSE))) {
      stale <- c(stale, file)
    }
  }
  # A failure names the files inst/stan has moved on from. Reinstall DRBayes.
  expect_identical(stale, character(0))
})


skip_without_cmdstan <- function() {
  skip_on_cran()
  skip_if_not_installed("instantiate")
  skip_if_not_installed("cmdstanr")
  skip_if_not(instantiate::stan_cmdstan_exists(),
              "CmdStan is not installed on this machine")
}

# Comparing two backends on data where sigma is one, and coefficients are of
# order one, cannot tell a prior variance measured in units of sigma^2 from one
# measured in units of Y: that is why the cross-backend suite passed against a
# stale model. Both tests below use the scale of the paper instead, where sigma
# is 40, and pass the precisions explicitly so that what is compared is the
# model each backend fits rather than the default each one chooses.

test_that("bayes_stan fits the model bayes_lm fits, at the paper's scale", {
  skip_without_cmdstan()
  d <- wide_scale_data()

  set.seed(51)
  gibbs <- pooled_draws(bayes_lm(d$Y, d$X, mc = 20000, chains = 4L,
                                 theta_prior = 0.01), burn = 5000)
  set.seed(51)
  stan <- bayes_stan("lm", d$Y, d$X, mc = 2000, warmup = 1000, chains = 4,
                     theta_prior = 0.01)

  spread <- apply(gibbs, 2, sd)
  difference <- abs(colMeans(gibbs) - apply(stan, 3, mean)) / spread
  expect_lt(max(difference), 0.05)
})


test_that("bayes_stan fits the model bayes_lm_hs fits, at the paper's scale", {
  skip_without_cmdstan()
  d <- wide_scale_data()

  set.seed(52)
  gibbs <- pooled_draws(bayes_lm_hs(d$Y, d$X, mc = 20000, chains = 4L,
                                    unshrunk = 1L, beta0_prior = 0.01,
                                    unshrunk_prior = 0.01), burn = 5000)
  set.seed(52)
  # A few divergences survive even at adapt_delta 0.99, which is what the
  # warning is for; the point of this test is that the means still line up.
  stan <- suppressWarnings(
    bayes_stan("lm_hs", d$Y, d$X, mc = 2000, warmup = 1000, chains = 4,
               unshrunk = 1L, beta0_prior = 0.01, unshrunk_prior = 0.01))

  spread <- apply(gibbs, 2, sd)
  difference <- abs(colMeans(gibbs) - apply(stan, 3, mean)) / spread
  expect_lt(max(difference), 0.1)
})

# The moment conditions on inputs chosen so that every intermediate value is
# exactly representable: entries drawn from {0, +-1/4, +-1/2, +-1, +-2}, a
# sample size that is a power of two so the division is exact, and a propensity
# score linear predictor that is identically zero so the fused weight is exactly
# plus or minus two.
#
# That matters because these are the only assertions in the suite that can be
# made bit-for-bit and still hold on another machine. Anywhere else, the answer
# depends on the order the BLAS accumulated its products in, on whether the
# compiler contracted a multiply and an add into one instruction, and on whether
# this build of R has long double. Here none of those can change the result, so
# expect_identical is a statement about the code rather than about the machine.

dyadic_case <- function(n = 8L, p_ps = 2L, p_otc = 3L, S = 5L, seed = 1L,
                        outcome_link = "identity", ps_link = "logit",
                        ps_at_zero = TRUE) {
  values <- c(0, 0.25, -0.25, 0.5, -0.5, 1, -1, 2)
  # The stride matters as much as the offset: the cycle has the same length as
  # the value set, so two offsets that differ by a multiple of it give the same
  # matrix, and a design matrix equal to its coefficients would make the fitted
  # means symmetric and hide a transposed read.
  pick <- function(k, offset, stride = 1L) {
    values[(((seq_len(k) * stride) + offset) %% length(values)) + 1L]
  }

  Z.ps  <- matrix(pick(n * p_ps, seed), n, p_ps)
  Z.lm  <- matrix(pick(n * p_otc, seed + 3L), n, p_otc)
  A     <- rep(c(1, 0), length.out = n)
  Y     <- pick(n, seed + 7L, stride = 5L)

  # Coefficients of zero on the propensity score side put every linear
  # predictor at zero, where exp() is exactly one and every value below stays a
  # binary fraction. That exactness is what lets the tests here assert identity
  # rather than nearness, and it costs something: at a linear predictor of zero
  # every inverse link returns one half, so the fused logistic weight, the
  # divided logistic weight and the probit weight all come to plus or minus two
  # and no fixture built this way can tell them apart. ps_at_zero = FALSE moves
  # off that point for the tests that need to.
  betas.ps  <- if (ps_at_zero) {
    matrix(0, S, p_ps)
  } else {
    matrix(pick(S * p_ps, seed + 5L, stride = 7L), S, p_ps)
  }
  betas.otc <- matrix(pick(S * p_otc, seed + 11L, stride = 3L), S, p_otc)

  inverse_link <- switch(outcome_link,
                         identity = function(eta) eta,
                         logit    = stats::plogis,
                         probit   = stats::pnorm)
  ps_inverse_link <- switch(ps_link,
                            logit  = stats::plogis,
                            probit = stats::pnorm,
                            # Numerically the logistic link to within one unit
                            # in the last place, but not identical() to it, so
                            # the general branch runs on the fused branch's
                            # numbers.
                            logit_general = function(eta) exp(eta) / (1 + exp(eta)))

  list(d = tilting_data(Z.lm = Z.lm, Z.ps = Z.ps, A = A, Y = Y,
                        inverse_link = inverse_link,
                        ps_inverse_link = ps_inverse_link,
                        ps_formula_text = "A ~ ."),
       betas.ps = betas.ps, betas.otc = betas.otc, n = n, S = S)
}


test_that("the moment condition is the mean over observations, not over draws", {
  # A wrong divisor rescales every entry of BB by the same factor. The tolerance
  # the sweep stops on is sd(BB)/sqrt(S), which rescales with it, so the
  # convergence decision is unchanged and only the treatment effect is wrong.
  # The fixture the suite had before this one has a single observation, where
  # the mean over observations is the identity, so it could not see this at all.
  fx <- dyadic_case(n = 8L, S = 5L)
  got <- moment_ipw(fx$betas.ps, fx$betas.otc, fx$d)

  # Computed from the definition rather than from either implementation: at a
  # zero linear predictor the propensity score is one half, so the weight is
  # 2 for a treated unit and -2 for a control one.
  weight <- ifelse(fx$d$A == 1, 2, -2)
  expected <- vapply(seq_len(fx$S), function(s) {
    mu <- drop(fx$d$Z.lm %*% fx$betas.otc[s, ])
    sum(weight * (fx$d$Y - mu)) / fx$n
  }, numeric(1))

  expect_identical(got, expected)
  expect_length(got, fx$S)
  expect_null(names(got))
})


test_that("the residual is broadcast down columns, not across rows", {
  # With as many draws as observations the two directions have the same shape,
  # so a transposed broadcast returns a plausible number instead of an error.
  fx <- dyadic_case(n = 8L, S = 8L)
  weight <- ifelse(fx$d$A == 1, 2, -2)
  expected <- vapply(seq_len(8L), function(s) {
    mu <- drop(fx$d$Z.lm %*% fx$betas.otc[s, ])
    sum(weight * (fx$d$Y - mu)) / 8
  }, numeric(1))

  expect_identical(moment_ipw(fx$betas.ps, fx$betas.otc, fx$d), expected)
})


test_that("block boundaries do not move the answer", {
  # draw_blocks fills 2^18 cells, so at n = 8 a block holds 32768 draws. The
  # two sizes here straddle that: one leaves a final block of a single draw,
  # the other divides exactly and must not emit an empty one.
  skip_on_cran()
  for (S in c(32767L, 32768L, 32769L, 65536L)) {
    fx <- dyadic_case(n = 8L, S = S)
    got <- moment_ipw(fx$betas.ps, fx$betas.otc, fx$d)

    weight <- ifelse(fx$d$A == 1, 2, -2)
    mu <- fx$d$Z.lm %*% t(fx$betas.otc)
    expected <- .colMeans(weight * (fx$d$Y - mu), fx$n, S)

    expect_identical(got, expected, info = paste("S =", S))
  }
})


test_that("the subclassification moment never touches the propensity score", {
  fx <- dyadic_case(n = 8L, S = 5L)
  fx$d$subclass_weight <- rep(c(2, -2, 4, -4), length.out = fx$n)

  expected <- vapply(seq_len(fx$S), function(s) {
    mu <- drop(fx$d$Z.lm %*% fx$betas.otc[s, ])
    sum(fx$d$subclass_weight * (fx$d$Y - mu)) / fx$n
  }, numeric(1))
  expect_identical(moment_subclass(fx$betas.ps, fx$betas.otc, fx$d), expected)

  # Equation (F.2) reads the stratum treated fraction, not the model, so a
  # propensity score that could not be evaluated at all must make no difference.
  broken <- fx$d
  broken$ps_inverse_link <- function(eta) stop("must not be called")
  expect_identical(moment_subclass(fx$betas.ps, fx$betas.otc, broken), expected)
})


test_that("the nonparametric moment reads its draws by row", {
  # outcome.mu is draws by observations, the transpose of every other matrix
  # here. With as many draws as observations a transposed read is not a
  # dimension error, just a different number.
  fx <- dyadic_case(n = 8L, S = 8L)
  mu <- fx$betas.otc %*% t(fx$d$Z.lm)
  np <- fx$d
  np$outcome.mu <- mu

  # The two paths reduce differently, with mean() against colMeans, so they are
  # compared to tolerance rather than bit for bit. On this fixture the values
  # are exactly representable and the two reductions happen to agree, which is
  # why the difference between them is pinned on random data elsewhere and not
  # here.
  expect_equal(moment_np(fx$betas.ps, np),
               moment_ipw(fx$betas.ps, fx$betas.otc, fx$d),
               tolerance = 1e-12)

  # A transposed read is a different number, not an error, at this shape.
  transposed <- np
  transposed$outcome.mu <- t(mu)
  expect_false(isTRUE(all.equal(moment_np(fx$betas.ps, transposed),
                                moment_np(fx$betas.ps, np))))
})


test_that("degenerate shapes return rather than error", {
  for (S in c(0L, 1L)) {
    fx <- dyadic_case(n = 8L, S = max(S, 1L))
    betas.ps  <- fx$betas.ps[seq_len(S), , drop = FALSE]
    betas.otc <- fx$betas.otc[seq_len(S), , drop = FALSE]
    got <- moment_ipw(betas.ps, betas.otc, fx$d)
    expect_length(got, S)
  }

  fx <- dyadic_case(n = 1L, p_ps = 1L, p_otc = 1L, S = 3L)
  expect_length(moment_ipw(fx$betas.ps, fx$betas.otc, fx$d), 3L)
})


test_that("the fused and general weights agree on the same numbers", {
  # The general branch is selected by a link that is not identical() to plogis
  # but equals it to within one unit in the last place, so the two branches are
  # compared on the same values rather than on two different fixtures. Applying
  # the fused identity to a probit link is wrong by about 15 per cent at a
  # linear predictor of one, with no error and no NaN to show for it.
  fused   <- dyadic_case(ps_link = "logit")
  general <- dyadic_case(ps_link = "logit_general")

  expect_false(is_logistic(general$d$ps_inverse_link))
  expect_true(is_logistic(fused$d$ps_inverse_link))
  expect_equal(moment_ipw(fused$betas.ps, fused$betas.otc, fused$d),
               moment_ipw(general$betas.ps, general$betas.otc, general$d),
               tolerance = 1e-13)
})


test_that("every link combination reaches the same definition", {
  skip_on_cran()
  for (otc in c("identity", "logit", "probit")) {
    for (ps in c("logit", "probit", "logit_general")) {
      fx <- dyadic_case(outcome_link = otc, ps_link = ps)
      got <- moment_ipw(fx$betas.ps, fx$betas.otc, fx$d)

      ps_value <- fx$d$ps_inverse_link(0)
      weight   <- (fx$d$A - ps_value) / (ps_value * (1 - ps_value))
      expected <- vapply(seq_len(fx$S), function(s) {
        mu <- fx$d$inverse_link(drop(fx$d$Z.lm %*% fx$betas.otc[s, ]))
        sum(weight * (fx$d$Y - mu)) / fx$n
      }, numeric(1))

      expect_equal(got, expected, tolerance = 1e-13,
                   info = paste(otc, ps, sep = " / "))
    }
  }
})


test_that("the propensity score link is honoured away from the origin", {
  # The grid above evaluates every link at a linear predictor of zero, where
  # all three agree, so it would accept a kernel that applied the logistic
  # weight whatever link it was given. Here the coefficients are not zero, the
  # three links disagree by much more than rounding, and the weight varies from
  # draw to draw rather than being fixed across them.
  for (ps in c("logit", "probit", "logit_general")) {
    for (otc in c("identity", "logit")) {
      fx  <- dyadic_case(ps_link = ps, outcome_link = otc, ps_at_zero = FALSE)
      got <- moment_ipw(fx$betas.ps, fx$betas.otc, fx$d)

      expected <- vapply(seq_len(fx$S), function(s) {
        score  <- fx$d$ps_inverse_link(drop(fx$d$Z.ps %*% fx$betas.ps[s, ]))
        weight <- (fx$d$A - score) / (score * (1 - score))
        mu     <- fx$d$inverse_link(drop(fx$d$Z.lm %*% fx$betas.otc[s, ]))
        sum(weight * (fx$d$Y - mu)) / fx$n
      }, numeric(1))

      # Not identity: under a logistic link the kernel forms the weight as
      # s (1 + exp(-s eta)) and the line above forms it as a quotient, which
      # differ in the last bits by construction.
      expect_equal(got, expected, tolerance = 1e-12,
                   info = paste(otc, ps, sep = " / "))
    }
  }
})


test_that("the positivity guard sees every draw of a block, not just the first", {
  # The guard is applied once per block, to the extremes of the whole block.
  # One that read only the first draw would let a violation through wherever
  # it actually arises, which is late: the sweep pushes the propensity score
  # coefficients outwards step by step, so the draw that separates the
  # treatment groups is never the first one.
  fx <- dyadic_case(S = 3L)
  betas <- fx$betas.ps
  betas[3L, 1L] <- 1000

  expect_error(moment_ipw(betas, fx$betas.otc, fx$d),
               "positivity violation")
})


random_case <- function(n, S, ps_inverse_link) {
  Z.ps <- cbind(1, matrix(stats::rnorm(n * 3), n, 3))
  Z.lm <- cbind(1, matrix(stats::rnorm(n * 2), n, 2))
  list(d = tilting_data(Z.lm = Z.lm, Z.ps = Z.ps,
                        A = stats::rbinom(n, 1, 0.5), Y = stats::rnorm(n),
                        inverse_link = function(eta) eta,
                        ps_inverse_link = ps_inverse_link,
                        ps_formula_text = "A ~ ."),
       bps  = matrix(stats::rnorm(S * 4, sd = 0.3), S),
       botc = matrix(stats::rnorm(S * 3, sd = 0.3), S))
}


test_that("each compiled path reduces the way the R it replaces reduces", {
  skip_on_cran()
  # Not a dyadic fixture: the point is the last bits, which exact values would
  # hide. It holds across platforms because both sides accumulate in whatever
  # width this build of R uses, and reach the same BLAS for the same products.
  #
  # The sample sizes are chosen to include several that are not a multiple of
  # the vector width, and both propensity score links are run. An earlier
  # version of this test used one sample size, four hundred, under the logistic
  # link, and that is the one combination in which a compiled kernel that let
  # the compiler fuse its multiply and add still agreed with R: the vectorised
  # body rounds the product, and only the scalar remainder, or the unvectorised
  # probit branch, contracts. It passed against a kernel that disagreed with R
  # at every sample size under a probit link.
  set.seed(21)
  for (ps in list(logit = stats::plogis, probit = stats::pnorm)) {
    for (n in c(17L, 23L, 40L, 400L)) {
      S  <- 120L
      fx <- random_case(n, S, ps)
      d  <- fx$d; bps <- fx$bps; botc <- fx$botc
      label <- paste("n =", n)

      expect_identical(moment_ipw(bps, botc, d), moment_ipw_r(bps, botc, d),
                       info = label)

      sub <- d
      sub$subclass_weight <- stats::rnorm(n)
      expect_identical(moment_subclass(bps, botc, sub),
                       moment_subclass_r(bps, botc, sub), info = label)

      # The nonparametric path reaches its mean through mean() rather than
      # colMeans, and the two differ in the last bits. Reducing it like the
      # others would move every value here by about 1e-14, which the cross-path
      # comparison comes nowhere near seeing.
      np <- d
      np$outcome.mu <- matrix(stats::rnorm(S * n, sd = 0.5), S, n)
      expect_identical(moment_np(bps, np), moment_np_r(bps, np), info = label)

      np_sub <- np
      np_sub$subclass_weight <- stats::rnorm(n)
      expect_identical(moment_np(bps, np_sub), moment_np_r(bps, np_sub),
                       info = label)
    }
  }
})


test_that("the nonparametric moment counts draws from the fitted means", {
  # The two matrices are made to agree by drbayes_pc, so a kernel that took its
  # draw count from the propensity score side instead would be inert everywhere
  # the package itself calls it, and wrong for anyone assembling the pieces.
  set.seed(22)
  n <- 60L; S <- 25L
  Z.ps <- cbind(1, matrix(stats::rnorm(n * 2), n, 2))
  d <- tilting_data(Z.lm = Z.ps, Z.ps = Z.ps,
                    A = stats::rbinom(n, 1, 0.5), Y = stats::rnorm(n),
                    inverse_link = function(eta) eta,
                    ps_inverse_link = stats::plogis, ps_formula_text = "A ~ .")
  d$outcome.mu <- matrix(stats::rnorm(S * n), S, n)
  spare <- matrix(stats::rnorm((S + 7L) * 3, sd = 0.3), S + 7L)

  expect_length(moment_np(spare, d), S)
  expect_identical(moment_np(spare, d), moment_np_r(spare, d))
})


test_that("each dimension check says what it means", {
  # Rcpp::stop takes a format string, so a second string is an argument to it
  # rather than a continuation of the message, and tinyformat reports its own
  # error instead of the sentence that was written. Five of these read
  # "tinyformat: Not enough conversion specifiers in format string" until the
  # commas came out. Asking only whether they threw could not tell the
  # difference, and that is what the sanitiser exercise had been asking.
  wide <- long_double_sums()
  Z <- matrix(1, 8L, 2L)
  b <- matrix(0, 3L, 2L)
  y <- stats::rnorm(8)

  expect_error(drb_moment_ipw(Z, Z, stats::rnorm(3), y, b, b, 1L, 0L,
                              moment_max_cells, wide),
               "must agree on the number of observations")
  expect_error(drb_moment_subclass(Z, y, stats::rnorm(3), b, 0L,
                                   moment_max_cells, wide),
               "stratum weights and Y must agree on the number of observations")
  expect_error(drb_ps_rowmeans(Z, matrix(0, 3L, 5L), 1L, moment_max_cells, wide),
               "one column per design matrix column")
  expect_error(drb_moment_np_ipw(Z, y, y, matrix(0, 3L, 5L), b, 1L, wide),
               "fitted means must agree on the number of observations")
  expect_error(drb_moment_np_ipw(Z, y, y, matrix(0, 2L, 8L),
                                 b[1L, , drop = FALSE], 1L, wide),
               "at least one propensity score draw per row of fitted means")
  expect_error(drb_moment_np_subclass(y, matrix(0, 3L, 5L), y, wide),
               "stratum weights must agree with Y")

  # This one had no check at all: the loop took its length from the treatment
  # and wrote into a buffer made from the linear predictors, so a longer
  # treatment wrote past the end.
  expect_error(drb_ipw_weight(stats::rnorm(3), stats::rnorm(8), 1L),
               "must be the same length")
})

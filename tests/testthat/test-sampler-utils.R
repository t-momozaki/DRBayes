# The helpers every Gibbs sampler is built from. Most of what is checked here
# is a refusal: these functions are where a prior that is not a prior, a
# starting value of the wrong length or a covariate matrix holding an Inf is
# supposed to be turned away, before it becomes a wrong number thousands of
# iterations later. Each refusal is reached and its message read, since a
# message naming the wrong cause sends a user to fix the wrong thing.
#
# The rest pin arithmetic against a second, independent route to the same
# answer: an explicit matrix inverse where the code uses a Cholesky factor,
# lm() where the code uses lm.fit(), and the formula of Piironen and Vehtari
# (2017) written out where the code evaluates it.

# The text of the error a call raises, or NA when it raises none.
error_text <- function(expr) {
  tryCatch({
    force(expr)
    NA_character_
  }, error = conditionMessage)
}

# Run `code` and put the session's random number state back afterwards,
# including the case where there was none to begin with. The tests below
# deliberately delete .Random.seed and switch generator kind, and neither may
# be left changed for the tests that follow.
with_rng_restored <- function(code) {
  had  <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  seed <- if (had) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  kind <- RNGkind()
  on.exit({
    RNGkind(kind[1L], kind[2L], kind[3L])
    if (had) {
      assign(".Random.seed", seed, envir = globalenv())
    } else {
      suppressWarnings(rm(list = ".Random.seed", envir = globalenv()))
    }
  }, add = TRUE)
  force(code)
}


test_that("validate_prior_precision accepts a scalar, a vector and a matrix", {
  # NULL takes the default, which is now one precision per coefficient: the
  # scale of a coefficient depends on the units of its column, so a single
  # number for all of them is what shrank a treatment effect of 110 by 3.
  per_coefficient <- c(1e-6, 4, 9)
  got <- validate_prior_precision(NULL, 3L, default = per_coefficient)
  expect_identical(dim(got), c(3L, 3L))
  expect_identical(diag(got), per_coefficient)
  expect_identical(sum(got[upper.tri(got) | lower.tri(got)]), 0)

  # A scalar default still means the same precision for every coefficient.
  expect_identical(validate_prior_precision(NULL, 2L, default = 1 / 100),
                   diag(0.01, 2L))
  expect_identical(validate_prior_precision(NULL, 2L), diag(1 / 100, 2L))

  # A scalar prior is a precision, so the prior standard deviation it implies
  # is 1/sqrt of it. Checking the variance rather than the number returned is
  # what catches a precision quietly used as a variance.
  scalar <- validate_prior_precision(1 / 25, 3L)
  expect_equal(sqrt(diag(solve(scalar))), rep(5, 3L))

  # A matrix is returned as given, correlations and all.
  full <- matrix(c(4, 1, 1, 2), 2L, 2L)
  expect_identical(validate_prior_precision(full, 2L), full)

  # A one by one matrix is a matrix, not a scalar, so it is checked for
  # dimension against d rather than broadcast to a diagonal.
  expect_identical(validate_prior_precision(matrix(4, 1L, 1L), 1L),
                   matrix(4, 1L, 1L))
  expect_error(validate_prior_precision(matrix(4, 1L, 1L), 2L),
               "theta_prior must be a 2 by 2 matrix, but is 1 by 1",
               fixed = TRUE)
})


test_that("validate_prior_precision refuses a matrix that is not a prior", {
  # Symmetry first: an asymmetric matrix is not a covariance structure at all,
  # and the message has to carry the argument name because the samplers pass
  # more than one prior through this function.
  asymmetric <- matrix(c(2, 1, 0, 2), 2L, 2L)
  expect_error(validate_prior_precision(asymmetric, 2L),
               "theta_prior must be symmetric", fixed = TRUE)
  expect_error(validate_prior_precision(asymmetric, 2L, arg = "beta0_prior"),
               "beta0_prior must be symmetric", fixed = TRUE)

  expect_error(validate_prior_precision(matrix(c(1, NA, NA, 1), 2L, 2L), 2L),
               "theta_prior must be finite", fixed = TRUE)
  expect_error(validate_prior_precision(matrix(c(1, Inf, Inf, 1), 2L, 2L), 2L),
               "theta_prior must be finite", fixed = TRUE)

  # An indefinite precision defines no proper prior. Refusing it here rather
  # than letting chol() fail mid-sweep is the point of the check, so the
  # message has to name the eigenvalue that made it indefinite.
  indefinite <- matrix(c(1, 2, 2, 1), 2L, 2L)
  said <- error_text(validate_prior_precision(indefinite, 2L))
  expect_true(grepl("theta_prior must be positive definite", said,
                    fixed = TRUE))
  expect_true(grepl("smallest eigenvalue is -1", said, fixed = TRUE))
  expect_true(grepl("does not define a proper prior", said, fixed = TRUE))

  # Positive semi-definite is not enough either: a zero eigenvalue is an
  # improper flat direction, and chol() would fail on it later.
  singular <- matrix(c(1, 1, 1, 1), 2L, 2L)
  expect_error(validate_prior_precision(singular, 2L, arg = "hs_prior"),
               "hs_prior must be positive definite", fixed = TRUE)
})


test_that("validate_prior_precision refuses a scalar that is not positive", {
  # The argument is a precision, and a user who passes a variance by mistake
  # gets a prior a hundred times too wide with no complaint, so the message
  # has to say which of the two the function wants.
  said <- error_text(validate_prior_precision(0, 2L))
  expect_true(grepl("is a prior PRECISION", said, fixed = TRUE))
  expect_true(grepl("Pass 1/variance", said, fixed = TRUE))
  expect_true(grepl("0 would give an improper prior", said, fixed = TRUE))
  expect_error(validate_prior_precision(-3, 2L, arg = "sigma_prior"),
               "sigma_prior is a prior PRECISION", fixed = TRUE)

  # The smallest legal scalar is any positive number.
  tiny <- validate_prior_precision(.Machine$double.xmin, 1L)
  expect_identical(tiny, diag(.Machine$double.xmin, 1L))

  # Anything that is neither a scalar nor a d by d matrix names both shapes,
  # since a length-d vector of precisions is a natural thing to try.
  for (bad in list(c(1, 2), character(1), NA_real_, Inf, list(1))) {
    expect_error(
      validate_prior_precision(bad, 2L),
      "must be NULL, a single positive number, or a 2 by 2 positive definite",
      fixed = TRUE)
  }
})


test_that("prior_sd_from_precision inverts the precision matrix", {
  # A correlated precision is where the two are different: the marginal
  # standard deviation is the square root of a diagonal entry of the INVERSE,
  # not the reciprocal square root of a diagonal entry of the precision.
  Q <- matrix(c(4, 1.5, 1.5, 2), 2L, 2L)
  expect_equal(prior_sd_from_precision(Q), sqrt(diag(solve(Q))))
  expect_false(isTRUE(all.equal(prior_sd_from_precision(Q),
                                1 / sqrt(diag(Q)))))

  # On a diagonal precision the two agree, which is why the mistake survives
  # every test that only uses a diagonal prior.
  D <- diag(c(4, 25))
  expect_equal(prior_sd_from_precision(D), c(0.5, 0.2))
})


test_that("columns_involving finds every column a variable contributes", {
  set.seed(3)
  dat <- data.frame(Y = stats::rnorm(30), A = stats::rbinom(30, 1L, 0.5),
                    X1 = stats::rnorm(30),
                    G = factor(sample(letters[1:3], 30, replace = TRUE)))
  mf <- stats::model.frame(Y ~ A * X1 + A:G + G, dat)
  tt <- stats::terms(mf)
  Z  <- stats::model.matrix(tt, mf)

  # A treatment that modifies its effect through an interaction contributes
  # more than one column. Missing the interaction columns is what collapses a
  # heterogeneous treatment effect onto a single coefficient, so the answer is
  # compared against the column names rather than against a hard-coded count.
  wanted <- which(colnames(Z) %in% c("A", "A:X1", "A:Gb", "A:Gc"))
  expect_identical(columns_involving(tt, Z, "A"), wanted)
  expect_gt(length(wanted), 1L)

  # A factor is spread over one column per non-reference level, and all of
  # them belong to the variable, including the ones inside the interaction.
  expect_identical(columns_involving(tt, Z, "G"),
                   which(colnames(Z) %in% c("Gb", "Gc", "A:Gb", "A:Gc")))

  # The intercept belongs to no variable.
  expect_false(1L %in% columns_involving(tt, Z, "A"))
  expect_false(1L %in% columns_involving(tt, Z, "X1"))

  # A variable that is not in the model contributes nothing, rather than
  # matching by accident on a name that is a substring of another.
  expect_identical(columns_involving(tt, Z, "W"), integer(0))
  expect_identical(columns_involving(tt, Z, "X"), integer(0))
})


test_that("columns_involving returns nothing when there is nothing to match", {
  set.seed(3)
  dat <- data.frame(Y = stats::rnorm(10))
  mf  <- stats::model.frame(Y ~ 1, dat)
  tt  <- stats::terms(mf)
  Z   <- stats::model.matrix(tt, mf)
  # An intercept-only formula has an empty factors table, which must not be
  # indexed into.
  expect_identical(columns_involving(tt, Z, "A"), integer(0))

  # A design matrix that lost its assign attribute, as one built by hand has,
  # cannot be mapped back to terms. Returning nothing leaves the shrinkage
  # prior shrinking everything, which is wrong but safe; guessing would not
  # be.
  dat2 <- data.frame(Y = stats::rnorm(10), A = stats::rbinom(10, 1L, 0.5))
  mf2  <- stats::model.frame(Y ~ A, dat2)
  tt2  <- stats::terms(mf2)
  Z2   <- stats::model.matrix(tt2, mf2)
  expect_identical(columns_involving(tt2, Z2, "A"), 2L)
  attr(Z2, "assign") <- NULL
  expect_identical(columns_involving(tt2, Z2, "A"), integer(0))
})


test_that("draw_storage allocates a named double draws_array", {
  X <- matrix(0, 5L, 2L, dimnames = list(NULL, c("age", "dose")))
  arr <- draw_storage(7L, 3L, X)

  # The posterior package's draws_array layout is iterations by chains by
  # parameters, and the diagnostics are handed this array unrearranged.
  expect_identical(dim(arr), c(7L, 3L, 3L))
  expect_identical(dimnames(arr)[[3]], c("(Intercept)", "age", "dose"))

  # Allocating logical NA would make the first assignment copy and coerce the
  # whole array, so the storage mode is part of the contract.
  expect_identical(typeof(arr), "double")
  expect_true(all(is.na(arr)))

  # Columns with no names, or with any name blank, are numbered rather than
  # left for array() to invent names of its own.
  expect_identical(dimnames(draw_storage(1L, 1L, matrix(0, 2L, 2L)))[[3]],
                   c("(Intercept)", "X1", "X2"))
  partly <- matrix(0, 2L, 2L, dimnames = list(NULL, c("age", "")))
  expect_identical(dimnames(draw_storage(1L, 1L, partly))[[3]],
                   c("(Intercept)", "X1", "X2"))

  # ps.formula = A ~ 1 gives a design with no columns at all. paste0("X",
  # integer(0)) is "X" rather than character(0), so the intercept-only case
  # would otherwise be handed one name too many and array() would refuse it.
  empty <- draw_storage(4L, 2L, matrix(numeric(0), nrow = 5L, ncol = 0L))
  expect_identical(dim(empty), c(4L, 2L, 1L))
  expect_identical(dimnames(empty)[[3]], "(Intercept)")
})


test_that("validate_chains takes a count and says why four are wanted", {
  expect_identical(validate_chains(1), 1L)
  expect_identical(validate_chains(4L), 4L)
  expect_identical(validate_chains(1e6), 1000000L)

  for (bad in list(0, -1, 2.5, NA_real_, NA, Inf, c(2, 2), "4", NULL)) {
    expect_error(validate_chains(bad),
                 "chains must be a single positive integer", fixed = TRUE)
  }
  # One chain is legal but cannot support R-hat, and the message has to say
  # so, since a user who sets chains = 1 to save time is about to lose the
  # diagnostic that decides whether the fit is reported as converged.
  expect_error(validate_chains(0), "R-hat compares chains against each other",
               fixed = TRUE)
})


test_that("validate_sampler_inputs refuses what would hang or mislead", {
  Y <- c(0, 1, 0, 1)
  X <- matrix(stats::rnorm(8), 4L, 2L)

  expect_error(validate_sampler_inputs(letters[1:4], X, 10L),
               "Y must be a numeric vector", fixed = TRUE)
  expect_error(validate_sampler_inputs(c(0, 1, 2, 1), X, 10L, binary = TRUE),
               "Y must be binary (0 or 1)", fixed = TRUE)
  # The binary check is skipped where it does not apply, so a continuous
  # response is not turned away from the linear sampler.
  expect_identical(dim(validate_sampler_inputs(c(0.5, 1.5, 2.5, 3.5), X, 10L)),
                   c(4L, 2L))

  expect_error(validate_sampler_inputs(Y, list(1, 2), 10L),
               "X must be a matrix or data.frame", fixed = TRUE)
  # A factor column survives as.matrix() as text, and the advice has to name
  # the fix rather than leave the user guessing at the coercion rule.
  factors <- data.frame(g = factor(c("a", "b", "a", "b")), z = 1:4)
  expect_error(validate_sampler_inputs(Y, factors, 10L),
               "X must be numeric", fixed = TRUE)
  expect_error(validate_sampler_inputs(Y, factors, 10L),
               "model.matrix()", fixed = TRUE)

  # A single NaN or Inf reaching pgdraw() makes its rejection sampler loop
  # forever, so the logistic samplers would hang rather than fail. The count
  # in the message is what tells a user how much of their data is affected.
  bad_X <- X
  bad_X[1L, 1L] <- NA
  bad_X[2L, 2L] <- Inf
  expect_error(validate_sampler_inputs(Y, bad_X, 10L),
               "X must be finite; found 2 missing or infinite value(s)",
               fixed = TRUE)
  expect_error(validate_sampler_inputs(c(0, 1, Inf, 1), X, 10L),
               "Y must be finite; found 1 missing or infinite value(s)",
               fixed = TRUE)

  expect_error(validate_sampler_inputs(Y, matrix(0, 5L, 2L), 10L),
               "X has 5 rows but Y has 4 elements", fixed = TRUE)

  for (bad in list(0, -1, 2.5, "10", c(10, 10))) {
    expect_error(validate_sampler_inputs(Y, X, bad),
                 "mc must be a single positive integer", fixed = TRUE)
  }
  # One iteration is legal here; drbayes_control() is where a number of
  # iterations too small to survive the burn-in is refused.
  expect_identical(dim(validate_sampler_inputs(Y, X, 1L)), c(4L, 2L))
})


test_that("validate_sampler_inputs returns X as a numeric matrix", {
  Y <- c(0.1, 0.2, 0.3, 0.4)
  frame <- data.frame(age = c(1, 2, 3, 4), dose = c(0.5, 0.5, 1, 1))
  got <- validate_sampler_inputs(Y, frame, 10L)
  # The samplers index the result as a matrix and name their draws from its
  # columns, so both the type and the column names have to survive.
  expect_true(is.matrix(got))
  expect_type(got, "double")
  expect_identical(colnames(got), c("age", "dose"))
  expect_equal(got[, "age"], frame$age, ignore_attr = TRUE)

  # ps.formula = A ~ 1 leaves a design with no columns, which is legitimate
  # and must not be mistaken for an empty matrix of the wrong shape.
  none <- validate_sampler_inputs(Y, matrix(numeric(0), nrow = 4L, ncol = 0L),
                                  10L)
  expect_identical(dim(none), c(4L, 0L))
})


test_that("validate_unshrunk names the indices it could not use", {
  expect_identical(validate_unshrunk(NULL, 5L), integer(0))
  expect_identical(validate_unshrunk(integer(0), 5L), integer(0))

  # Duplicates and any order are accepted and normalised, since the caller
  # builds this from columns_involving() and may pass the same column twice.
  expect_identical(validate_unshrunk(c(3, 1, 1), 5L), c(1L, 3L))
  expect_identical(validate_unshrunk(2, 5L), 2L)
  expect_type(validate_unshrunk(c(2, 4), 5L), "integer")

  # Names are what a user reaches for first, so the refusal has to say the
  # function wants indices rather than merely that the type is wrong.
  expect_error(validate_unshrunk("A", 5L),
               "unshrunk must be given as column indices, not names",
               fixed = TRUE)
  expect_error(validate_unshrunk(c("A", "X1"), 5L),
               "not names", fixed = TRUE)

  expect_error(validate_unshrunk(1.5, 5L),
               "unshrunk must be whole numbers indexing the columns of X",
               fixed = TRUE)
  expect_error(validate_unshrunk(TRUE, 5L), "whole numbers", fixed = TRUE)

  # Off-by-one is the mistake this catches: unshrunk indexes the columns of X,
  # not of the design matrix that carries the intercept in front of them, so
  # the message lists exactly which values fell outside.
  expect_error(validate_unshrunk(c(0, 2, 6), 5L),
               "unshrunk indexes columns outside X, which has 5 columns: 0, 6",
               fixed = TRUE)
  expect_identical(validate_unshrunk(c(1, 5), 5L), c(1L, 5L))
  expect_error(validate_unshrunk(6, 5L), "columns: 6", fixed = TRUE)
  expect_error(validate_unshrunk(0, 5L), "columns: 0", fixed = TRUE)
})


test_that("resolve_tau_prior follows Piironen and Vehtari (2017)", {
  # tau_0 = p0 / (p - p0) * sigma / sqrt(n), written out here rather than
  # taken from the code, so a rearranged formula is a failure and not a new
  # expected value.
  expect_equal(resolve_tau_prior(NULL, p = 50, n = 200, p0 = 5, sigma = 2),
               5 / (50 - 5) * 2 / sqrt(200))
  # With no error scale the dimensional part stands on its own.
  expect_equal(resolve_tau_prior(NULL, p = 50, n = 200, p0 = 5, sigma = NULL),
               5 / (50 - 5) * 1 / sqrt(200))

  # The prior scale shrinks as the sample grows and as the guessed number of
  # signals falls, which is the whole point of the recommendation.
  expect_lt(resolve_tau_prior(NULL, 50, 2000, 5, 1),
            resolve_tau_prior(NULL, 50, 200, 5, 1))
  expect_lt(resolve_tau_prior(NULL, 50, 200, 2, 1),
            resolve_tau_prior(NULL, 50, 200, 5, 1))

  # Fewer coefficients than guessed signals leaves nothing to shrink towards,
  # and p0 / (p - p0) would be negative or undefined there, so the unit-scale
  # half-Cauchy is used instead. The boundary is at p == p0.
  expect_identical(resolve_tau_prior(NULL, p = 5, n = 200, p0 = 5, sigma = 1),
                   1)
  expect_identical(resolve_tau_prior(NULL, p = 3, n = 200, p0 = 5, sigma = 1),
                   1)
  expect_lt(resolve_tau_prior(NULL, p = 6, n = 200, p0 = 5, sigma = 1), 1)

  # A supplied scale is returned untouched, whatever p, n and p0 say.
  expect_identical(resolve_tau_prior(0.05, p = 5, n = 200, p0 = 5, sigma = 1),
                   0.05)
  for (bad in list(0, -1, Inf, NA_real_, c(0.1, 0.1), "0.1")) {
    expect_error(resolve_tau_prior(bad, 50, 200, 5, 1),
                 "tau_prior must be a single positive number", fixed = TRUE)
  }

  # p0 counts coefficients, so it has to be a whole number of them, and the
  # message has to say what it counts or the name means nothing.
  for (bad in list(0, -1, 2.5, c(5, 5), "5")) {
    expect_error(resolve_tau_prior(NULL, 50, 200, bad, 1),
                 "p0 must be a single positive integer", fixed = TRUE)
  }
  expect_error(resolve_tau_prior(NULL, 50, 200, 0, 1),
               "a guess at how many coefficients are non-zero", fixed = TRUE)
  expect_identical(resolve_tau_prior(NULL, p = 2, n = 200, p0 = 1, sigma = 1),
                   1 / (2 - 1) * 1 / sqrt(200))
})


test_that("draw_normal_precision samples the distribution it claims to", {
  Q <- matrix(c(4, 1.5, 1.5, 3), 2L, 2L)
  b <- c(2, -1)

  # The Cholesky route is compared against an explicit inverse: the mean is
  # Q^-1 b and the noise is U^-1 z. A transposed solve would still return
  # numbers of a plausible size, so the draw is reproduced exactly rather
  # than only checked for shape.
  set.seed(21)
  z <- stats::rnorm(2L)
  wanted <- drop(solve(Q) %*% b + solve(chol(Q)) %*% z)
  set.seed(21)
  expect_equal(draw_normal_precision(Q, b, iter = 1L), wanted)

  # And the second moment, which is where using Q where Q^-1 belongs shows
  # up: the covariance of the draws is the INVERSE of the precision.
  set.seed(22)
  draws <- t(replicate(6000L, draw_normal_precision(Q, b, iter = 1L)))
  expect_equal(colMeans(draws), drop(solve(Q) %*% b), tolerance = 0.02)
  expect_equal(stats::cov(draws), solve(Q), tolerance = 0.06)
})


test_that("draw_normal_precision explains a precision that lost rank", {
  # Collinear columns and an underflowed shrinkage parameter both end here,
  # and chol()'s own message says nothing about either. The iteration number
  # is the part that lets a user find where the chain went wrong.
  singular <- matrix(c(1, 1, 1, 1), 2L, 2L)
  said <- error_text(draw_normal_precision(singular, c(1, 1), iter = 137L))
  expect_true(grepl("stopped being positive definite at iteration 137", said,
                    fixed = TRUE))
  expect_true(grepl("shrinkage parameter", said, fixed = TRUE))
  expect_true(grepl("collinear", said, fixed = TRUE))
  # chol()'s own message is kept, not thrown away, so the diagnosis can be
  # checked against it.
  expect_true(grepl("Original message:", said, fixed = TRUE))

  # The error is raised with call. = FALSE, so what a user sees is the
  # explanation and not an unreadable internal call.
  condition <- tryCatch(draw_normal_precision(singular, c(1, 1), 3L),
                        error = function(e) e)
  expect_null(conditionCall(condition))
})


test_that("with_preserved_rng leaves the caller's stream where it found it", {
  with_rng_restored({
    # CRAN policy forbids a function changing the user's random number state.
    # A script that seeds itself once at the top has to keep control of every
    # draw after this call.
    set.seed(101)
    wanted <- stats::rnorm(3L)
    set.seed(101)
    invisible(with_preserved_rng(7, stats::rnorm(50L)))
    expect_identical(stats::rnorm(3L), wanted)

    # The state is restored byte for byte, not merely made to look plausible.
    set.seed(101)
    before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    invisible(with_preserved_rng(7, stats::rnorm(50L)))
    expect_identical(
      get(".Random.seed", envir = globalenv(), inherits = FALSE), before)

    # Inside the call the seed is in force, so the value is exactly what the
    # caller would get from set.seed() themselves.
    set.seed(7)
    expect_identical(with_preserved_rng(7, stats::rnorm(4L)), stats::rnorm(4L))
  })
})


test_that("with_preserved_rng restores the generator kind as well", {
  with_rng_restored({
    # Switching kind re-seeds the generator, so the kind has to go back
    # before the seed does. Restoring them the other way round would leave a
    # correct kind beside a random state, which no test of the seed alone
    # would notice.
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    # Read rather than written out: R-devel reports a fourth component, for the
    # binomial method, and a test naming three would pass while the fourth was
    # being lost.
    before_kind <- RNGkind()
    set.seed(3)
    before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    invisible(with_preserved_rng(7, {
      RNGkind("L'Ecuyer-CMRG")
      stats::rnorm(5L)
    }))
    expect_identical(RNGkind(), before_kind)
    expect_identical(
      get(".Random.seed", envir = globalenv(), inherits = FALSE), before)
  })
})


test_that("with_preserved_rng leaves an unseeded session unseeded", {
  with_rng_restored({
    # A fresh session has no .Random.seed at all. Creating one is a visible
    # change to the caller's state, and it also silently fixes the stream for
    # everything that follows.
    suppressWarnings(rm(list = ".Random.seed", envir = globalenv()))
    expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    invisible(with_preserved_rng(11, stats::rnorm(3L)))
    expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))

    # Same again with the kind changed inside, since restoring the kind is
    # itself a way to create a .Random.seed that then has to be removed.
    suppressWarnings(rm(list = ".Random.seed", envir = globalenv()))
    invisible(with_preserved_rng(11, {
      RNGkind("L'Ecuyer-CMRG")
      stats::rnorm(3L)
    }))
    expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    expect_identical(RNGkind()[1L], "Mersenne-Twister")
  })
})


test_that("quick_coef_fit reproduces the frequentist fit it stands in for", {
  set.seed(5)
  x <- stats::rnorm(80L)
  Y <- 1 + 2 * x + stats::rnorm(80L)
  Z <- cbind(1, x)

  # The starting values are placed at this fit, so it has to be the fit and
  # not something close to it: compared against lm(), which arrives by a
  # different route and applies the residual scale itself.
  reference <- summary(stats::lm(Y ~ x))$coefficients
  got <- quick_coef_fit(Y, Z, "gaussian", "logit")
  expect_equal(got$coefficients, unname(reference[, 1L]))
  expect_equal(got$se, unname(reference[, 2L]))

  A <- stats::rbinom(80L, 1L, stats::plogis(0.3 + 0.8 * x))
  logit_fit <- stats::glm(A ~ x, family = stats::binomial())
  logit_reference <- summary(logit_fit)$coefficients
  logit <- quick_coef_fit(A, Z, "binomial", "logit")
  expect_equal(logit$coefficients, unname(logit_reference[, 1L]))
  expect_equal(logit$se, unname(logit_reference[, 2L]))

  # The link is honoured. A probit fit passed off as a logit one is a defect
  # this package has had, and the two give visibly different slopes.
  probit <- quick_coef_fit(A, Z, "binomial", "probit")
  expect_false(isTRUE(all.equal(probit$coefficients, logit$coefficients)))
  probit_reference <- stats::glm(A ~ x,
                                 family = stats::binomial("probit"))
  expect_equal(probit$coefficients, unname(stats::coef(probit_reference)))
})


test_that("quick_coef_fit gives up rather than start a chain from nonsense", {
  # An unusable fit must not stop a sampler that would have run perfectly
  # well from a poorer start, so every one of these is a NULL and not an
  # error. The caller falls back on the prior scale.
  collinear <- cbind(1, c(1, 2, 3, 4), c(2, 4, 6, 8))
  expect_null(quick_coef_fit(c(1, 2, 3, 4), collinear, "gaussian", "logit"))
  expect_null(quick_coef_fit(c(0, 1, 0, 1), collinear, "binomial", "logit"))

  # Complete separation: the logistic fit runs out of iterations without
  # converging, and its coefficients are on their way to infinity.
  x  <- seq(-1, 1, length.out = 40L)
  Zs <- cbind(1, x)
  A  <- as.numeric(x > 0)
  expect_false(suppressWarnings(
    stats::glm.fit(Zs, A, family = stats::binomial())$converged))
  expect_null(quick_coef_fit(A, Zs, "binomial", "logit"))
  expect_null(quick_coef_fit(A, Zs, "binomial", "probit"))

  # A full rank design can still give standard errors so large that a
  # starting value drawn on their scale is pure noise.
  set.seed(6)
  flat <- 1 + stats::rnorm(60L, sd = 1e-5)
  wild <- stats::rnorm(60L, sd = 1e5)
  Zf <- cbind(1, flat)
  expect_identical(suppressWarnings(stats::lm.fit(Zf, wild)$rank), 2L)
  expect_gt(max(summary(stats::lm(wild ~ flat))$coefficients[, 2L]), 1e6)
  expect_null(quick_coef_fit(wild, Zf, "gaussian", "logit"))
})


test_that("resolve_inits overdisperses the chains around the fit", {
  set.seed(8)
  x <- stats::rnorm(120L)
  Y <- 1 + 2 * x + stats::rnorm(120L)
  Z <- cbind(1, x)
  fit <- quick_coef_fit(Y, Z, "gaussian", "logit")

  # Chain one sits at the estimate so that at least one chain starts where
  # the data support, and the rest are drawn at `spread` standard errors, the
  # dispersion Vehtari et al. (2021) section 2 asks for. Reproduced here
  # from the caller's stream, which is also the claim that no seed is set
  # inside.
  set.seed(31)
  got <- resolve_inits(NULL, 3L, Y, Z, "gaussian", "logit", spread = 4)
  set.seed(31)
  wanted <- list(fit$coefficients,
                 fit$coefficients + 4 * fit$se * stats::rnorm(2L),
                 fit$coefficients + 4 * fit$se * stats::rnorm(2L))
  expect_equal(got, wanted)
  expect_length(got, 3L)

  # An unusable fit leaves the chains centred on zero and spread on the prior
  # scale, one entry per coefficient rather than one for all of them.
  collinear <- cbind(1, c(1, 2, 3, 4), c(2, 4, 6, 8))
  set.seed(32)
  fallback <- resolve_inits(NULL, 2L, c(1, 2, 3, 4), collinear, "gaussian",
                            "logit", prior_sd = c(5, 6, 7))
  set.seed(32)
  expect_equal(fallback[[1L]], c(0, 0, 0))
  expect_equal(fallback[[2L]], c(5, 6, 7) * stats::rnorm(3L))
})


test_that("resolve_inits checks supplied starting values against the design", {
  Z <- cbind(1, stats::rnorm(20L), stats::rnorm(20L))
  Y <- stats::rnorm(20L)
  ok <- list(c(0, 0, 0), c(1, 1, 1))

  expect_identical(resolve_inits(ok, 2L, Y, Z), ok)
  # Names and attributes are stripped, so a named vector cannot leak into the
  # draws array and rename its parameters.
  named <- list(c(a = 0, b = 0, c = 0), c(1, 1, 1))
  expect_identical(resolve_inits(named, 2L, Y, Z), ok)

  # A single starting value still has to be wrapped in a list, and the
  # message has to show the wrapping rather than only refuse.
  expect_error(resolve_inits(c(0, 0, 0), 1L, Y, Z),
               "still has to be wrapped in a list, as in list(rep(0, 3))",
               fixed = TRUE)
  expect_error(resolve_inits(ok, 4L, Y, Z),
               "init has 2 element(s) but chains is 4", fixed = TRUE)

  # The length is checked against the design matrix, which is the check
  # drbayes_control() cannot make because it has not built one yet. The
  # message spells out that the intercept comes first, since an init of
  # length ncol(X) rather than ncol(X) + 1 is the usual mistake.
  expect_error(resolve_inits(list(c(0, 0, 0), c(1, 1)), 2L, Y, Z),
               "init[[2]] must be a numeric vector of length 3", fixed = TRUE)
  expect_error(resolve_inits(list(c(0, 0, 0), c(1, 1)), 2L, Y, Z),
               "the intercept followed by the 2 column(s) of X", fixed = TRUE)
  expect_error(resolve_inits(list(c(0, 0, 0), letters[1:3]), 2L, Y, Z),
               "of class character and length 3", fixed = TRUE)
  expect_error(resolve_inits(list(c(0, 0, 0), c(1, NA, Inf)), 2L, Y, Z),
               "init[[2]] must be finite; found 2 missing or infinite",
               fixed = TRUE)
})


test_that("%||% supplies a default only for NULL", {
  # NA, 0, FALSE and the empty vector are all values a caller may have meant,
  # so only NULL may be replaced.
  expect_identical(NULL %||% 5, 5)
  expect_identical(3 %||% 5, 3)
  expect_identical(NA %||% 5, NA)
  expect_identical(0 %||% 5, 0)
  expect_identical(FALSE %||% 5, FALSE)
  expect_identical(integer(0) %||% 5, integer(0))
})


test_that("a precision matrix carrying dimnames is accepted", {
  # isSymmetric() compares attributes, and transposing swaps dimnames, so a
  # matrix labelled on one margin only is not symmetric to it. Dropping the
  # unname() guard turns a precision read from a file into "not a prior".
  prior <- matrix(c(4, 1, 1, 2), 2L, 2L)
  colnames(prior) <- c("age", "dose")
  expect_identical(dim(validate_prior_precision(prior, 2L)), c(2L, 2L))

  both <- prior
  rownames(both) <- colnames(both)
  expect_identical(dim(validate_prior_precision(both, 2L)), c(2L, 2L))
})


test_that("the caller's generator kind is restored, not reset to the default", {
  # Restoring a hard-coded Mersenne-Twister would look correct in a session
  # that already used it, which is every session the rest of the suite runs in.
  # future and furrr both set L'Ecuyer-CMRG, and both are Imports.
  old <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old)), add = TRUE)

  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  invisible(with_preserved_rng(7, stats::rnorm(3)))
  expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")

  RNGkind("Wichmann-Hill", "Inversion", "Rejection")
  invisible(with_preserved_rng(7, stats::rnorm(3)))
  expect_identical(RNGkind()[1], "Wichmann-Hill")
})

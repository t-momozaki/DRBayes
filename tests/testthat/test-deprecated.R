# The functions renamed for the first CRAN release still have to work under
# their old names. A shim that warns but computes something else is worse than
# no shim at all, because the warning is what persuades the reader that the
# answer below it is the same one. Each shim is therefore checked against the
# replacement its own warning names, called with the same arguments from the
# same seed.

old_names <- c("DRBayes.PC", "DRBayes.BB", "DRBayes.BB.nleqslv", "B.LM",
               "B.Logit", "B.Probit", "HS.LM", "HS.Logit", "HS.Probit")

sampler_shims <- list(B.LM      = "bayes_lm",
                      B.Logit   = "bayes_logit",
                      B.Probit  = "bayes_probit",
                      HS.LM     = "bayes_lm_hs",
                      HS.Logit  = "bayes_logit_hs",
                      HS.Probit = "bayes_probit_hs")

binary_outcome <- c("B.Logit", "B.Probit", "HS.Logit", "HS.Probit")

deprecated_test_data <- function(n = 120L, seed = 4L) {
  set.seed(seed)
  X <- cbind(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  eta <- 0.4 + 0.8 * X[, 1L] - 0.5 * X[, 2L]
  list(X = X,
       gaussian = as.numeric(eta + stats::rnorm(n)),
       binary   = stats::rbinom(n, 1L, stats::plogis(eta)),
       frame    = data.frame(YY = as.numeric(eta + stats::rnorm(n)),
                             AA = stats::rbinom(n, 1L, stats::plogis(eta)),
                             W1 = X[, 1L], W2 = X[, 2L]))
}

run_with_warnings <- function(expr) {
  seen <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(cnd) {
      seen <<- c(seen, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  list(value = value, warnings = seen)
}

deprecation_warnings <- function(warnings) {
  grep("is deprecated", warnings, value = TRUE)
}

# The replacement named in the warning, which is the promise the shim makes to
# the reader and so the function its answer is checked against.
replacement_named_in <- function(message) {
  sub("^.*Use '([^']+)' instead.*$", "\\1", gsub("\n", " ", message))
}


test_that("every renamed function is still reachable under its old name", {
  for (name in old_names) {
    expect_true(is.function(get0(name, mode = "function")), info = name)
  }
})


test_that("each renamed sampler forwards to the function its warning names", {
  d <- deprecated_test_data()

  for (old in names(sampler_shims)) {
    args <- list(Y = if (old %in% binary_outcome) d$binary else d$gaussian,
                 X = d$X, mc = 200L, chains = 2L)

    set.seed(11)
    run <- run_with_warnings(do.call(old, args))
    deprecations <- deprecation_warnings(run$warnings)

    # One warning, no more: a shim that warned twice, or that let the
    # replacement warn again on its way out, trains the reader to ignore it.
    expect_length(deprecations, 1L)
    expect_match(deprecations, paste0("'", old, "' is deprecated"),
                 fixed = TRUE, info = old)

    # The name in the warning has to be the name of the function that actually
    # produced the answer, which catches both a shim pointing at the wrong
    # replacement and one that forwards somewhere its message never mentions.
    named <- replacement_named_in(deprecations)
    expect_identical(named, sampler_shims[[old]], info = old)

    set.seed(11)
    expect_identical(run$value, do.call(named, args), info = old)

    # mc and chains reached the sampler rather than falling back to its
    # defaults of 5000 and 4, which would otherwise be invisible: the answer
    # would still be a plausible set of posterior draws.
    expect_identical(dim(run$value)[1:2], c(200L, 2L), info = old)

    # The replacement itself must be free of the deprecation warning, or every
    # caller of the new name is told to stop using it.
    set.seed(11)
    direct <- run_with_warnings(do.call(named, args))
    expect_length(deprecation_warnings(direct$warnings), 0L)
  }
})


test_that("arguments beyond the data reach the replacement", {
  # A dropped argument is the failure that matters here: it leaves the old name
  # quietly computing under a different prior from the new one, and a prior is
  # nowhere visible in the array of draws that comes back.
  d <- deprecated_test_data()

  set.seed(12)
  shimmed <- suppressWarnings(B.LM(d$gaussian, d$X, mc = 200L, chains = 2L,
                                   sigma_prior = c(4, 7)))
  set.seed(12)
  direct <- bayes_lm(d$gaussian, d$X, mc = 200L, chains = 2L,
                     sigma_prior = c(4, 7))
  expect_identical(shimmed, direct)

  set.seed(12)
  default_prior <- bayes_lm(d$gaussian, d$X, mc = 200L, chains = 2L)
  expect_false(identical(shimmed, default_prior))

  # The same, for an argument that changes which coefficients are shrunk.
  set.seed(13)
  shimmed_hs <- suppressWarnings(HS.LM(d$gaussian, d$X, mc = 200L,
                                       chains = 2L, unshrunk = 1L))
  set.seed(13)
  direct_hs <- bayes_lm_hs(d$gaussian, d$X, mc = 200L, chains = 2L,
                           unshrunk = 1L)
  expect_identical(shimmed_hs, direct_hs)

  set.seed(13)
  all_shrunk <- bayes_lm_hs(d$gaussian, d$X, mc = 200L, chains = 2L)
  expect_false(identical(shimmed_hs, all_shrunk))

  # An argument the replacement does not have has to fail rather than be
  # swallowed by the shim's own dots.
  expect_error(suppressWarnings(B.LM(d$gaussian, d$X, mc = 200L,
                                     shrinkage = 0.5)),
               "unused argument")
})


test_that("both Bayesian bootstrap names give the answer drbayes_bb gives", {
  d <- deprecated_test_data()
  args <- list(YY ~ AA + W1 + W2, AA ~ W1 + W2, d$frame,
               num_iterations = 200L, verbose = FALSE)

  set.seed(21)
  bb <- run_with_warnings(do.call("DRBayes.BB", args))
  expect_length(deprecation_warnings(bb$warnings), 1L)
  expect_match(deprecation_warnings(bb$warnings), "'DRBayes.BB' is deprecated",
               fixed = TRUE)
  expect_identical(replacement_named_in(deprecation_warnings(bb$warnings)),
                   "drbayes_bb")

  set.seed(21)
  expect_identical(bb$value, do.call("drbayes_bb", args))

  set.seed(21)
  nleqslv <- run_with_warnings(do.call("DRBayes.BB.nleqslv", args))
  expect_length(deprecation_warnings(nleqslv$warnings), 1L)
  expect_match(deprecation_warnings(nleqslv$warnings),
               "'DRBayes.BB.nleqslv' is deprecated", fixed = TRUE)
  expect_identical(replacement_named_in(deprecation_warnings(
    nleqslv$warnings)), "drbayes_bb")

  # The root finder was only ever another way to solve the same two weighted
  # regressions, which is why it was not reimplemented. If the two old names
  # ever stop agreeing, one of them is no longer the estimator it documents.
  expect_identical(nleqslv$value, bb$value)
})


test_that("DRBayes.PC gives the answer drbayes_pc gives", {
  skip_on_cran()

  d <- deprecated_test_data()
  args <- list(YY ~ AA + W1 + W2, AA ~ W1 + W2, d$frame,
               outcome.model = bayes_lm, ps.model = bayes_logit,
               mc = 600L, bn = 200L, chains = 2L, verbose = FALSE,
               control = drbayes_control(n_steps = 50))

  set.seed(31)
  utils::capture.output(
    shimmed <- run_with_warnings(do.call("DRBayes.PC", args)),
    file = nullfile())
  deprecations <- deprecation_warnings(shimmed$warnings)
  expect_length(deprecations, 1L)
  expect_match(deprecations, "'DRBayes.PC' is deprecated", fixed = TRUE)
  expect_identical(replacement_named_in(deprecations), "drbayes_pc")

  set.seed(31)
  utils::capture.output(
    direct <- suppressWarnings(do.call("drbayes_pc", args)),
    file = nullfile())

  # The call each records differs, because one of them was reached through the
  # shim's dots. Everything the estimator produced must not.
  expect_s3_class(shimmed$value, "DRBayes")
  expect_identical(shimmed$value$pc, direct$pc)
  expect_identical(shimmed$value$g.comp, direct$g.comp)
  expect_identical(shimmed$value$smc, direct$smc)
  expect_identical(shimmed$value$diagnostics, direct$diagnostics)
})

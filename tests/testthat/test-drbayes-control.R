# drbayes_control() is the package's contract with its users: a value it
# accepts is a value the sweep will act on, and a value it refuses has to say
# which setting was wrong and why. The tests below reach every refusal and
# assert on what it blames, pin the boundary on each side of every numeric
# setting, and pin the settings that are legal but INERT at a boundary. A
# setting that looks switched on while doing nothing is the failure this
# package has actually shipped, so those cases are pinned rather than left to
# be rediscovered.

# One observation, an intercept-only propensity score and an intercept-only
# outcome model, which is the smallest thing the tilting sweep will run on.
sweep_fixture <- function(S = 40, seed = 1) {
  set.seed(seed)
  betas.ps  <- matrix(stats::rnorm(S, sd = 0.2), ncol = 1L)
  betas.otc <- matrix(stats::rnorm(S, mean = 1, sd = 0.5), ncol = 1L)
  d <- tilting_data(Z.lm = matrix(1, 1L, 1L), Z.ps = matrix(1, 1L, 1L),
                    A = 1, Y = 0,
                    inverse_link = function(x) x,
                    ps_inverse_link = stats::plogis,
                    ps_formula_text = "A ~ 1")
  list(betas.ps = betas.ps, betas.otc = betas.otc, d = d,
       BB = moment_ipw(betas.ps, betas.otc, d))
}

# Four steps of Algorithm 2 on that fixture, with a tolerance no cloud will
# meet so that the sweep always takes them all.
sweep_particles <- function(fx, ..., seed = 4) {
  ctrl <- drbayes_control(n_steps = 4, lambda_max = 0.4, ...)
  set.seed(seed)
  out <- tilt_smc(fx$betas.ps, fx$betas.otc, fx$BB, fx$d, ctrl,
                  tol = 1e-12, direction = 1, moment_fn = moment_ipw)
  cbind(out$betas.ps, out$betas.otc)
}

# TRUE when every row of `after` is a row of `before`, which is what pure
# resampling leaves behind and what any added jitter destroys.
all_rows_from <- function(after, before) {
  all(apply(after, 1L, function(row) {
    any(apply(before, 1L, function(old) isTRUE(all.equal(row, old))))
  }))
}

# Did drbayes_control() accept this combination of settings?
accepts <- function(...) {
  !inherits(try(drbayes_control(...), silent = TRUE), "try-error")
}

# The text of the error a call raises, or NA when it raises none.
error_text <- function(expr) {
  tryCatch({
    force(expr)
    NA_character_
  }, error = conditionMessage)
}


test_that("every setting drbayes_control() accepts is returned", {
  ctrl <- drbayes_control()
  expect_s3_class(ctrl, "drbayes_control")

  # A setting missing from the returned list is one drbayes_pc() can never
  # read, so a user who names it gets silence instead of the behaviour asked
  # for. The two lists have to match name for name and in order.
  expect_identical(names(ctrl), names(formals(drbayes_control)))

  # The settings that default to NULL have to survive as elements rather than
  # disappear from the list, since drbayes_pc() reads them by name.
  expect_true(all(c("tol", "init") %in% names(ctrl)))
  expect_null(ctrl$tol)
  expect_null(ctrl$init)

  # Every setting given has to arrive in the returned list carrying the value
  # it was given, so one that is validated and then dropped is caught.
  given <- list(mc = 200, bn = 50, thin = 3, chains = 2L, lambda_max = 5,
                n_steps = 250, newton_steps = 7, smoothing = 0.95,
                tol = 0.02, ridge = 1e-8, ess_frac = 0.25,
                moment = "subclass", n_subclass = 4, pruning = TRUE,
                prune.rule = "quantile", prune.w = 0.4, prune.q = 0.2,
                rhat_max = 1.05, ess_min = 100, keep_particles = TRUE)
  set  <- do.call(drbayes_control, given)
  for (setting in names(given)) {
    expect_equal(set[[setting]], given[[setting]], label = setting)
  }
})


test_that("the print method shows the pruning rule actually in force", {
  # The printed summary is how a user checks what they asked for. Showing the
  # quantile rule's setting while the weight rule is the one running tells
  # them the sweep is doing something it is not.
  weight <- capture.output(print(drbayes_control(pruning = TRUE,
                                                 prune.w = 0.25)))
  expect_true(any(grepl("pruning = on (weight, prune.w = 0.25/S)", weight,
                        fixed = TRUE)))
  expect_false(any(grepl("prune.q", weight, fixed = TRUE)))

  quant <- capture.output(print(drbayes_control(pruning = TRUE,
                                                prune.rule = "quantile",
                                                prune.q = 0.2)))
  expect_true(any(grepl("pruning = on (quantile, prune.q = 0.2)", quant,
                        fixed = TRUE)))
  expect_false(any(grepl("prune.w", quant, fixed = TRUE)))

  # With pruning off neither threshold is in force, so neither is reported.
  off <- capture.output(print(drbayes_control(prune.w = 0.25, prune.q = 0.2)))
  expect_true(any(grepl("pruning = off", off, fixed = TRUE)))
  expect_false(any(grepl("prune.w", off, fixed = TRUE)))
  expect_false(any(grepl("prune.q", off, fixed = TRUE)))
})


test_that("the print method reports the tolerance and moment in force", {
  # tol = NULL means the Monte Carlo standard error, which is a rule and not a
  # number, so printing a number there would misreport the stopping rule.
  auto <- capture.output(print(drbayes_control()))
  expect_true(any(grepl("tol = Monte Carlo standard error", auto,
                        fixed = TRUE)))
  fixed_tol <- capture.output(print(drbayes_control(tol = 0.001)))
  expect_true(any(grepl("tol = 0.001", fixed_tol, fixed = TRUE)))

  # n_subclass only means anything under the subclassification moment, so it
  # is shown there and nowhere else.
  ipw <- capture.output(print(drbayes_control(n_subclass = 8)))
  expect_true(any(grepl("moment = ipw", ipw, fixed = TRUE)))
  expect_false(any(grepl("n_subclass", ipw, fixed = TRUE)))
  sub <- capture.output(print(drbayes_control(moment = "subclass",
                                              n_subclass = 8)))
  expect_true(any(grepl("subclass (n_subclass = 8)", sub, fixed = TRUE)))

  # The step size is lambda_max / n_steps, the quantity the sweep actually
  # takes, rather than either setting on its own.
  step <- capture.output(print(drbayes_control(lambda_max = 5,
                                               n_steps = 250)))
  expect_true(any(grepl("(step 0.02)", step, fixed = TRUE)))

  # The method returns its input invisibly, as print methods must, so that
  # ctrl <- print(ctrl) is not a way to lose a control object.
  shown <- capture.output(seen <- withVisible(print(drbayes_control())))
  expect_length(shown, 9L)
  expect_false(seen$visible)
  expect_identical(seen$value, drbayes_control())
})


test_that("the print method says which model the starting values are for", {
  two <- list(c(0, 0), c(1, 1))
  both <- capture.output(print(drbayes_control(
    chains = 2L, init = list(outcome = two, ps = two))))
  expect_true(any(grepl("init = supplied per model", both, fixed = TRUE)))

  one <- capture.output(print(drbayes_control(chains = 2L,
                                              init = list(ps = two))))
  expect_true(any(grepl("init = supplied for the ps model", one,
                        fixed = TRUE)))

  shared <- capture.output(print(drbayes_control(chains = 2L, init = two)))
  expect_true(any(grepl("init = supplied", shared, fixed = TRUE)))
  expect_false(any(grepl("per model", shared, fixed = TRUE)))

  auto <- capture.output(print(drbayes_control()))
  expect_true(any(grepl("init = automatic", auto, fixed = TRUE)))
})


test_that("mc, bn and thin are refused exactly when under two draws remain", {
  # drbayes_pc() retains seq(bn + 1, mc, by = thin). The control object has to
  # refuse precisely the combinations that rule leaves fewer than two draws
  # in, otherwise a fit either dies deep inside the sampler or runs on a
  # posterior of one point.
  grid <- expand.grid(mc = 2:10, bn = 0:4, thin = 1:3)
  wanted <- mapply(function(mc, bn, thin) {
    mc > bn + 1L && length(seq(bn + 1L, mc, by = thin)) >= 2L
  }, grid$mc, grid$bn, grid$thin)
  accepted <- mapply(function(mc, bn, thin) accepts(mc = mc, bn = bn,
                                                    thin = thin),
                     grid$mc, grid$bn, grid$thin)
  expect_identical(accepted, wanted)

  # The two refusals say different things and must not be confused: one is
  # about the burn-in swallowing everything, the other about the thinning.
  expect_error(drbayes_control(mc = 1001, bn = 1000),
               "mc must exceed bn + 1", fixed = TRUE)
  expect_error(drbayes_control(mc = 1001, bn = 1000),
               "mc = 1001, bn = 1000", fixed = TRUE)
  expect_error(drbayes_control(mc = 1002, bn = 1000, thin = 2),
               "leave only 1 draw per chain", fixed = TRUE)
  expect_error(drbayes_control(mc = 1002, bn = 1000, thin = 2),
               "Raise mc, or lower bn or thin", fixed = TRUE)

  # The smallest legal setting of each, and the first illegal one below it.
  expect_identical(drbayes_control(mc = 2, bn = 0, thin = 1)$mc, 2L)
  expect_identical(drbayes_control(bn = 0)$bn, 0L)
  expect_identical(drbayes_control(thin = 1)$thin, 1L)
  expect_error(drbayes_control(mc = 0),
               "mc must be a single positive integer", fixed = TRUE)
  expect_error(drbayes_control(bn = -1),
               "bn must be a single non-negative integer", fixed = TRUE)
  expect_error(drbayes_control(thin = 0),
               "thin must be a single positive integer", fixed = TRUE)
})


test_that("count settings refuse non-integers and report the value seen", {
  # The message has to name the setting and show what arrived, since a count
  # given as 2.5 or as a vector is nearly always a mis-typed argument name.
  expect_error(drbayes_control(thin = 2.5), "Received 2.5", fixed = TRUE)
  expect_error(drbayes_control(mc = c(10, 20)),
               "Received an object of class numeric and length 2",
               fixed = TRUE)
  expect_error(drbayes_control(n_steps = Inf), "Received Inf", fixed = TRUE)
  expect_error(drbayes_control(mc = NA_integer_), "Received NA", fixed = TRUE)
  expect_error(drbayes_control(mc = list(5000)),
               "Received an object of class list and length 1", fixed = TRUE)
  expect_error(drbayes_control(mc = "5000"),
               "mc must be a single positive integer", fixed = TRUE)

  # describe_value() shows a scalar and falls back on class and length for
  # anything that would print badly.
  expect_identical(describe_value(2.5), "2.5")
  expect_identical(describe_value(c(1, 2)),
                   "an object of class numeric and length 2")
  expect_identical(describe_value(list(1, 2)),
                   "an object of class list and length 2")
})


test_that("chains reaches the sampler as a whole number of chains", {
  # A count kept as a double would be handed to seq_len() and to the draws
  # array unchanged, so the coercion is part of the setting and not a detail.
  expect_identical(drbayes_control(chains = 4)$chains, 4L)
  expect_identical(drbayes_control(chains = 1L)$chains, 1L)
  for (bad in list(0, 2.5, NA_real_, "4")) {
    expect_error(drbayes_control(chains = bad),
                 "chains must be a single positive integer", fixed = TRUE)
  }
})


test_that("init is checked for shape and reports which chain is wrong", {
  two <- list(c(0, 0), c(1, 1))

  expect_null(drbayes_control()$init)
  expect_identical(drbayes_control(chains = 2L, init = two)$init, two)
  expect_error(drbayes_control(init = c(0, 0)),
               "init must be NULL or a list, but is of class numeric",
               fixed = TRUE)

  # Naming the models is all or nothing, and the message has to list the names
  # that were not understood rather than merely refuse.
  expect_error(
    drbayes_control(chains = 2L, init = list(outcome = two, otc = two)),
    "every element must be called", fixed = TRUE)
  expect_error(
    drbayes_control(chains = 2L, init = list(outcome = two, otc = two)),
    "Found: otc.", fixed = TRUE)

  # One starting value per chain, with the offending count in the message.
  expect_error(drbayes_control(chains = 4L, init = two),
               "init has 2 element(s) but chains is 4", fixed = TRUE)
  expect_error(
    drbayes_control(chains = 4L, init = list(ps = two)),
    "init$ps has 2 element(s) but chains is 4", fixed = TRUE)
  expect_error(drbayes_control(chains = 2L, init = list(ps = c(0, 0))),
               "init$ps must be a list holding one starting value per chain",
               fixed = TRUE)

  # Which element is wrong, not just that one is.
  expect_error(
    drbayes_control(chains = 2L, init = list(c(0, 0), "a")),
    "init[[2]] must be a non-empty numeric vector", fixed = TRUE)
  expect_error(drbayes_control(chains = 2L, init = list(c(0, 0), numeric(0))),
               "init[[2]] must be a non-empty numeric vector", fixed = TRUE)
  expect_error(drbayes_control(chains = 2L, init = list(c(0, NA), c(1, 1))),
               "init[[1]] must be finite; found 1 missing", fixed = TRUE)
  expect_error(
    drbayes_control(chains = 2L,
                    init = list(outcome = list(c(0, Inf), c(1, 1)),
                                ps = two)),
    "init$outcome[[1]] must be finite; found 1 missing", fixed = TRUE)
})


test_that("control_init_for hands each model its own starting values", {
  two   <- list(c(0, 0), c(1, 1))
  other <- list(c(2, 2), c(3, 3))

  # Two shapes of init mean different things, and reading the wrong one would
  # start the propensity score sampler from the outcome model's values.
  named <- drbayes_control(chains = 2L,
                           init = list(outcome = two, ps = other))
  expect_identical(control_init_for(named, "outcome"), two)
  expect_identical(control_init_for(named, "ps"), other)

  shared <- drbayes_control(chains = 2L, init = two)
  expect_identical(control_init_for(shared, "outcome"), two)
  expect_identical(control_init_for(shared, "ps"), two)

  # A control naming only one model leaves the other on the automatic values.
  one <- drbayes_control(chains = 2L, init = list(ps = other))
  expect_null(control_init_for(one, "outcome"))
  expect_identical(control_init_for(one, "ps"), other)

  expect_null(control_init_for(drbayes_control(), "outcome"))
})


test_that("lambda_max, n_steps and newton_steps are bounded", {
  # lambda_max and n_steps are the range and the resolution of one grid, so
  # the step size is what a user has changed when either moves.
  ctrl <- drbayes_control(lambda_max = 5, n_steps = 250)
  expect_identical(ctrl$lambda_max / ctrl$n_steps, 0.02)
  expect_identical(drbayes_control()$lambda_max /
                     drbayes_control()$n_steps, 0.01)

  expect_true(accepts(lambda_max = .Machine$double.eps))
  expect_true(accepts(lambda_max = 1e12))
  for (bad in list(0, -1, Inf, NA_real_, c(1, 2), "10")) {
    expect_error(drbayes_control(lambda_max = bad),
                 "lambda_max must be a single positive number", fixed = TRUE)
  }

  expect_identical(drbayes_control(n_steps = 1)$n_steps, 1L)
  expect_error(drbayes_control(n_steps = 0), "n_steps must be a single",
               fixed = TRUE)
  expect_error(drbayes_control(n_steps = 0),
               "the number of steps the sweep takes to reach lambda_max",
               fixed = TRUE)

  # newton_steps caps Algorithm 1 and n_steps resolves Algorithm 2's grid.
  # They are separate settings and their messages must not point at the same
  # algorithm.
  expect_identical(drbayes_control(newton_steps = 1)$newton_steps, 1L)
  expect_error(drbayes_control(newton_steps = 0),
               "the largest number of Newton iterations Algorithm 1 may take",
               fixed = TRUE)
})


test_that("smoothing accepts 1, where the Liu and West kernel does nothing", {
  expect_identical(drbayes_control(smoothing = 1)$smoothing, 1)
  expect_true(accepts(smoothing = .Machine$double.eps))
  for (bad in list(0, -0.5, 1 + 1e-8, 2, NA_real_, Inf, c(0.9, 0.9), "0.9")) {
    expect_error(drbayes_control(smoothing = bad),
                 "smoothing must be a single number in (0, 1]", fixed = TRUE)
  }

  skip_on_cran()
  # At a = 1 the kernel shrinks by nothing and adds noise with variance
  # (1 - a^2) = 0, so a step of the sweep is plain resampling: every particle
  # returned is a copy of one it started from and the cloud loses distinct
  # values. The setting is legal, and this is what it costs.
  fx     <- sweep_fixture()
  before <- cbind(fx$betas.ps, fx$betas.otc)
  sharp  <- sweep_particles(fx, smoothing = 1)
  expect_true(all_rows_from(sharp, before))
  expect_lt(nrow(unique(sharp)), nrow(before))

  # Just below 1 the kernel jitters every particle, so no row survives intact.
  soft <- sweep_particles(fx, smoothing = 0.99)
  expect_false(all_rows_from(soft, before))
  expect_identical(nrow(unique(soft)), nrow(before))
})


test_that("tol and ridge are positive or, for tol, deliberately absent", {
  expect_null(drbayes_control()$tol)
  expect_identical(drbayes_control(tol = 1e-6)$tol, 1e-6)
  for (bad in list(0, -1, Inf, NA_real_, c(1, 2), "1e-6")) {
    expect_error(drbayes_control(tol = bad),
                 "tol must be NULL or a single positive number", fixed = TRUE)
  }
  # NULL is a rule, not a missing value, and the message has to say which.
  expect_error(drbayes_control(tol = 0), "Monte Carlo standard error",
               fixed = TRUE)

  expect_identical(drbayes_control(ridge = 1e-12)$ridge, 1e-12)
  for (bad in list(0, -1e-5, Inf, NA_real_, c(1, 2))) {
    expect_error(drbayes_control(ridge = bad),
                 "ridge must be a single positive number", fixed = TRUE)
  }
  # Zero is the value a user reaches for when they want the ridge gone, so the
  # refusal has to offer the thing they actually want.
  expect_error(drbayes_control(ridge = 0), "1e-10", fixed = TRUE)
})


test_that("ess_frac = 1 is legal and is a floor no sweep can meet", {
  expect_identical(drbayes_control(ess_frac = 1)$ess_frac, 1)
  expect_true(accepts(ess_frac = 1e-8))
  for (bad in list(0, -0.1, 1 + 1e-8, NA_real_, Inf, c(0.1, 0.1))) {
    expect_error(drbayes_control(ess_frac = bad),
                 "ess_frac must be a single number in (0, 1]", fixed = TRUE)
  }

  # Resampling S particles with replacement leaves about 63% of them
  # distinct, so a floor of the whole cloud reports every sweep as degenerate
  # however well the moment condition was met. Legal, but never satisfied.
  expect_true(tilting_degeneracy(NA_real_, 999L, 1000L, 1)$degenerate)
  expect_false(tilting_degeneracy(NA_real_, 1000L, 1000L, 1)$degenerate)
  expect_false(tilting_degeneracy(NA_real_, 999L, 1000L, 0.1)$degenerate)

  # The report has to name the floor that was applied, otherwise a user cannot
  # tell which of ess_frac and ess_min they should be looking at.
  flagged <- tilting_degeneracy(NA_real_, 12L, 1000L, 0.1)
  expect_true(grepl("ess_frac = 0.1", flagged$message, fixed = TRUE))
  expect_true(grepl("12 of the 1000 particles distinct", flagged$message,
                    fixed = TRUE))
})


test_that("moment picks one of two named conditions", {
  expect_identical(drbayes_control()$moment, "ipw")
  expect_identical(drbayes_control(moment = "subclass")$moment, "subclass")
  # The formals hold both, so the unmodified default has to resolve to the
  # first rather than being stored as a length-two vector.
  expect_identical(drbayes_control(moment = c("subclass", "ipw"))$moment,
                   "subclass")

  # Case matters, and the message has to show what arrived so that "IPW" is
  # recognisable as the mistake it is.
  expect_error(drbayes_control(moment = "IPW"), "Received IPW", fixed = TRUE)
  expect_error(drbayes_control(moment = "ips"), "equation (3.4)", fixed = TRUE)
  expect_error(drbayes_control(moment = character(0)),
               "an object of class character and length 0", fixed = TRUE)
  expect_error(drbayes_control(moment = 1), "Received 1", fixed = TRUE)
})


test_that("n_subclass is bounded below at two and is not bounded above", {
  expect_identical(drbayes_control(n_subclass = 2)$n_subclass, 2L)
  expect_error(drbayes_control(n_subclass = 1),
               "the number of propensity score strata, which cannot be ",
               fixed = TRUE)
  expect_error(drbayes_control(n_subclass = 0), "n_subclass", fixed = TRUE)

  # Whether a stratification is usable depends on the data, which the control
  # object has not seen, so a stratum per observation is accepted here and
  # refused where equation (F.2) would divide by zero.
  n <- 40L
  expect_identical(drbayes_control(n_subclass = n)$n_subclass, n)

  set.seed(2)
  X1 <- stats::rnorm(n)
  A  <- stats::rbinom(n, 1L, stats::plogis(0.4 * X1))
  d  <- tilting_data(Z.lm = cbind(1, X1), Z.ps = cbind(1, X1), A = A,
                     Y = stats::rnorm(n),
                     inverse_link = function(x) x,
                     ps_inverse_link = stats::plogis,
                     ps_formula_text = "A ~ X1")
  betas.ps <- matrix(c(0.1, 0.4), nrow = 3L, ncol = 2L, byrow = TRUE)

  expect_length(subclass_weights(betas.ps, d, 5L), n)
  # One observation per stratum leaves every stratum all treated or all
  # control, so the refusal has to name the stratum and the reason.
  expect_error(subclass_weights(betas.ps, d, n),
               "Propensity score stratum 1 of 40", fixed = TRUE)
  expect_error(subclass_weights(betas.ps, d, n),
               "use a smaller n_subclass", fixed = TRUE)
})


test_that("pruning thresholds that cannot act are inert, not refused", {
  # prune.q is a share of the cloud, so any value below 1/S rounds down to
  # discarding nothing: pruning = TRUE is then bit-identical to pruning =
  # FALSE. This is the shape of a defect the package has shipped before, a
  # threshold outside the range where the rule can act, so it is pinned.
  ww <- normalised_weights(seq(0, 3, length.out = 5))
  quantile_at <- function(q) {
    drbayes_control(pruning = TRUE, prune.rule = "quantile", prune.q = q)
  }
  expect_identical(floor(0.1 * length(ww)), 0)
  expect_identical(prune_weights(ww, ww, quantile_at(0.1)), ww)
  expect_false(identical(prune_weights(ww, ww, quantile_at(0.2)), ww))

  # The weight rule thresholds at prune.w / S, so with even weights nothing
  # falls below it until prune.w passes 1, and then everything does.
  even <- rep(1 / 20, 20)
  weight_at <- function(w) drbayes_control(pruning = TRUE, prune.w = w)
  expect_identical(prune_weights(even, even, weight_at(1)), even)
  expect_error(prune_weights(even, even, weight_at(1 + 1e-8)),
               "Sample pruning left 0 particle(s)", fixed = TRUE)

  # The refusal has to blame the rule in force. Naming prune.q while the
  # weight rule is running would send a user to tune a setting that is not
  # being read.
  said <- error_text(prune_weights(even, even, weight_at(1 + 1e-8)))
  expect_true(grepl("Lower prune.w, which is 1", said, fixed = TRUE))
  expect_false(grepl("prune.q", said, fixed = TRUE))

  collapsed <- normalised_weights(c(0, -800, -800, -800))
  said_q <- error_text(prune_weights(collapsed, collapsed, quantile_at(0.5)))
  expect_true(grepl("Lower prune.q, which is 0.5", said_q, fixed = TRUE))
  expect_false(grepl("prune.w", said_q, fixed = TRUE))
})


test_that("the pruning settings are inert while pruning is off", {
  skip_on_cran()
  # prune.rule, prune.w and prune.q are stored whatever pruning says, so a
  # sweep with pruning off has to ignore them completely. If it did not, a
  # user who left an old threshold in place would get a different fit from
  # one who never set it.
  fx <- sweep_fixture()
  plain  <- sweep_particles(fx)
  unused <- sweep_particles(fx, prune.rule = "quantile", prune.q = 0.9,
                            prune.w = 0.9)
  expect_identical(plain, unused)

  # Switching pruning on with the same thresholds does change the answer, so
  # the comparison above is not vacuous.
  pruned <- sweep_particles(fx, pruning = TRUE, prune.rule = "quantile",
                            prune.q = 0.5)
  expect_false(identical(plain, pruned))
})


test_that("pruning, prune.rule, prune.w and prune.q are validated", {
  expect_false(drbayes_control()$pruning)
  for (bad in list(NA, "TRUE", 1, c(TRUE, TRUE), NULL)) {
    expect_error(drbayes_control(pruning = bad),
                 "pruning must be TRUE or FALSE", fixed = TRUE)
  }

  expect_identical(drbayes_control()$prune.rule, "weight")
  # match.arg() abbreviates, so a shortened rule has to resolve rather than
  # silently fall back on the default.
  expect_identical(drbayes_control(prune.rule = "q")$prune.rule, "quantile")
  expect_error(drbayes_control(prune.rule = "smallest"), "quantile")

  expect_true(accepts(prune.w = .Machine$double.eps))
  expect_true(accepts(prune.w = 1e6))
  for (bad in list(0, -0.3, Inf, NA_real_, c(0.3, 0.3), "0.3")) {
    expect_error(drbayes_control(prune.w = bad),
                 "prune.w must be a single positive number", fixed = TRUE)
  }

  # prune.q = 0 is legal and documented as leaving the sweep untouched; 1
  # would delete the whole cloud and is not.
  expect_identical(drbayes_control(prune.q = 0)$prune.q, 0)
  expect_true(accepts(prune.q = 1 - 1e-12))
  for (bad in list(-1e-8, 1, 1.5, Inf, NA_real_, c(0.1, 0.1))) {
    expect_error(drbayes_control(prune.q = bad),
                 "prune.q must be a single number in [0, 1)", fixed = TRUE)
  }
  expect_error(drbayes_control(prune.q = 1), "Received 1", fixed = TRUE)
})


test_that("rhat_max and ess_min keep the diagnostic thresholds meetable", {
  # R-hat is at least 1 by construction, so a threshold of 1 can never be met
  # and would report every fit as unconverged.
  expect_error(drbayes_control(rhat_max = 1),
               "a threshold of 1 or less can never be met", fixed = TRUE)
  expect_error(drbayes_control(rhat_max = 0.99), "rhat_max", fixed = TRUE)
  expect_true(accepts(rhat_max = 1 + 1e-12))
  expect_identical(drbayes_control()$rhat_max, 1.01)
  for (bad in list(Inf, NA_real_, c(1.01, 1.01), "1.01")) {
    expect_error(drbayes_control(rhat_max = bad),
                 "rhat_max must be a single number greater than 1",
                 fixed = TRUE)
  }

  expect_true(accepts(ess_min = .Machine$double.eps))
  expect_identical(drbayes_control()$ess_min, 400)
  for (bad in list(0, -1, Inf, NA_real_, c(400, 400))) {
    expect_error(drbayes_control(ess_min = bad),
                 "ess_min must be a single positive number", fixed = TRUE)
  }
})


test_that("keep_particles is a flag and nothing else", {
  expect_false(drbayes_control()$keep_particles)
  expect_true(drbayes_control(keep_particles = TRUE)$keep_particles)
  for (bad in list(NA, "TRUE", 1, c(TRUE, FALSE), NULL)) {
    expect_error(drbayes_control(keep_particles = bad),
                 "keep_particles must be TRUE or FALSE", fixed = TRUE)
  }
})


test_that("resolve_control revalidates whatever it is handed", {
  expect_identical(resolve_control(NULL), drbayes_control())
  expect_identical(resolve_control(list()), drbayes_control())
  expect_identical(resolve_control(list(n_steps = 200)),
                   drbayes_control(n_steps = 200))
  expect_identical(resolve_control(drbayes_control(n_steps = 200)),
                   drbayes_control(n_steps = 200))

  # An object edited by hand after construction has to be checked again, not
  # trusted because of its class, or an impossible setting reaches the
  # sampler with the class vouching for it.
  tampered <- drbayes_control()
  tampered$mc <- -5
  expect_error(resolve_control(tampered),
               "mc must be a single positive integer", fixed = TRUE)

  expect_error(resolve_control(1:3),
               "but is of class integer", fixed = TRUE)
  expect_error(resolve_control(list(2000)),
               "must be a NAMED list", fixed = TRUE)
  expect_error(resolve_control(list(mc = 10, 2000)),
               "1 of its 2 element(s) have no name", fixed = TRUE)

  # A misspelt setting must be refused rather than ignored, and the message
  # has to list what was available so the right name can be found.
  expect_error(resolve_control(list(mcmc = 10)),
               "unknown setting(s): mcmc", fixed = TRUE)
  expect_error(resolve_control(list(mcmc = 10)),
               "The available settings are: mc, bn, thin", fixed = TRUE)
  expect_error(resolve_control(list(mc = 10, mc = 20)),
               "sets the same setting more than once: mc", fixed = TRUE)

  # The argument name travels into the message, since drbayes_pc() resolves
  # more than one control-shaped argument.
  expect_error(resolve_control(list(mcmc = 10), arg = "ps.control"),
               "ps.control contains unknown setting", fixed = TRUE)
})


test_that("every setting prints the value it holds", {
  # Nine of the thirteen printed settings were covered only by a line count, so
  # a print statement reading the wrong field, or one field printed in another's
  # place, went unnoticed. The printed block is how a user checks what they
  # asked for, so it is pinned as a whole against distinctive values.
  ctl <- drbayes_control(mc = 4321, bn = 765, thin = 3, chains = 5L,
                         lambda_max = 7.5, n_steps = 250, newton_steps = 33,
                         smoothing = 0.95, tol = 0.002, ridge = 3e-4,
                         ess_frac = 0.25, n_subclass = 6L, rhat_max = 1.05,
                         ess_min = 250, keep_particles = TRUE)
  shown <- paste(utils::capture.output(print(ctl)), collapse = " ")

  for (pair in list(c("mc", "4321"), c("bn", "765"), c("thin", "3"),
                    c("chains", "5"), c("lambda_max", "7.5"),
                    c("n_steps", "250"), c("newton_steps", "33"),
                    c("smoothing", "0.95"), c("ridge", "3e-04"),
                    c("ess_frac", "0.25"),
                    c("rhat_max", "1.05"), c("ess_min", "250"),
                    c("keep_particles", "TRUE"))) {
    expect_match(shown, paste0(pair[1], " = ", pair[2]), fixed = TRUE,
                 info = pair[1])
  }

  # The stratum count belongs to the subclassification moment, so it is shown
  # only when that moment is the one in force.
  expect_no_match(shown, "n_subclass", fixed = TRUE)
  sub <- utils::capture.output(
    print(drbayes_control(moment = "subclass", n_subclass = 6L)))
  expect_match(paste(sub, collapse = " "), "n_subclass = 6", fixed = TRUE)
})


test_that("a count past the integer bound is refused rather than made NA", {
  # as.integer() turns anything past the bound into NA with only a warning, so
  # a setting that had passed every guard came back missing, and mc and bn died
  # later in a base R comparison naming neither the setting nor the package.
  for (arg in c("mc", "bn", "thin", "chains", "n_steps", "newton_steps",
                "n_subclass")) {
    expect_error(do.call(drbayes_control, stats::setNames(list(3e9), arg)),
                 arg, fixed = TRUE, info = arg)
  }
  expect_identical(drbayes_control(n_subclass = 6L)$n_subclass, 6L)
})

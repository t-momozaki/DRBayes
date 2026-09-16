# Algorithm 4: which covariates drbayes_select() keeps, which coefficients the
# sweep is allowed to move, and what the resulting estimand is.

# Regularization-induced confounding, the setting of Hahn et al. (2018) that
# section 7.2 motivates the algorithm with. W1 to W6 are weak confounders:
# each predicts treatment, and each has an outcome coefficient small enough
# relative to the residual scale that the horseshoe shrinks it away. Z1 to Zk
# are noise. The unit level treatment effect is exactly `ate`.
weak_confounder_data <- function(seed, n = 250, n_weak = 6, n_noise = 24,
                                 alpha = 0.8, beta = 0.5, sigma = 3,
                                 ate = 2) {
  set.seed(seed)
  X <- matrix(rnorm(n * (n_weak + n_noise)), n, n_weak + n_noise)
  colnames(X) <- c(paste0("W", seq_len(n_weak)), paste0("Z", seq_len(n_noise)))
  W <- X[, seq_len(n_weak), drop = FALSE]
  d <- as.data.frame(X)
  d$A <- rbinom(n, 1, plogis(drop(W %*% rep(alpha, n_weak))))
  d$Y <- ate * d$A + drop(W %*% rep(beta, n_weak)) + rnorm(n, sd = sigma)
  d
}

covariate_formulas <- function(d) {
  covariates <- paste(setdiff(names(d), c("A", "Y")), collapse = " + ")
  list(outcome = stats::as.formula(paste("Y ~ A +", covariates)),
       ps      = stats::as.formula(paste("A ~", covariates)))
}

fit_quietly <- function(...) {
  utils::capture.output(
    res <- suppressWarnings(drbayes_select(..., verbose = FALSE)),
    file = nullfile())
  res
}

# The outcome columns the sweep may move, recomputed from the fit the way
# drbayes_select() computes them.
updated_columns <- function(fit, md) {
  ps_columns <- 1L + match(fit$selected, colnames(md$Z.ps)[-1L])
  sort(unique(c(1L, md$treatment_columns + 1L,
                selection_columns(ps_columns, md)$columns)))
}

# How many distinct values each coefficient takes across the particles. A
# coefficient the kernel refreshed is distinct at every particle; a frozen one
# takes as many values as the cloud has ancestors.
distinct_per_column <- function(particles) {
  apply(particles, 2L, function(x) length(unique(x)))
}


test_that("the fit carries the selection alongside the usual components", {
  skip_on_cran()
  d <- weak_confounder_data(1, n = 150, n_weak = 3, n_noise = 5)
  f <- covariate_formulas(d)

  fit <- fit_quietly(f$outcome, f$ps, d, threshold = 0.05, seed = 1,
                     mc = 600, bn = 200, thin = 2, chains = 2L,
                     control = drbayes_control(n_steps = 50),
                     diagnostics = "none")

  expect_s3_class(fit, "DRBayes_select")
  expect_s3_class(fit, "DRBayes")
  expect_identical(fit$threshold, 0.05)
  expect_identical(fit$n_candidates, 8L)
  expect_identical(fit$rejuvenate, "all")
  expect_length(fit$posterior_means, 8L)
  expect_named(fit$posterior_means, c(paste0("W", 1:3), paste0("Z", 1:5)))
  expect_true(all(fit$selected %in% names(fit$posterior_means)))
  # The selected set is exactly the covariates above the threshold.
  above <- abs(fit$posterior_means) >= 0.05
  expect_identical(fit$selected, names(fit$posterior_means)[above])

  # Inherited methods still work, which is the point of the class vector.
  expect_length(fit$pc, length(fit$g.comp))
  expect_s3_class(summary(fit), "summary.DRBayes")
  # The diagnostics run_tilting() records, so print() and summary() read a
  # selection fit exactly as they read a drbayes_pc() one, plus the ancestry
  # count that only a partly frozen particle needs.
  expect_named(fit$smc, c("lambda", "B_mean", "B_mean_weighted", "tol",
                          "n_steps", "max_steps", "converged", "ess",
                          "n_ancestors_min", "n_distinct",
                          "n_ancestors_cumulative"))
})


test_that("print reports the selection and what the sweep moved", {
  skip_on_cran()
  d <- weak_confounder_data(2, n = 150, n_weak = 3, n_noise = 5)
  f <- covariate_formulas(d)

  fit <- fit_quietly(f$outcome, f$ps, d, threshold = 0.05, seed = 2,
                     mc = 600, bn = 200, thin = 2, chains = 2L,
                     control = drbayes_control(n_steps = 50),
                     diagnostics = "none")

  out <- paste(utils::capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Confounder selection")
  expect_match(out, paste0("Selected ", length(fit$selected), " of 8"))
  expect_match(out, "Average treatment effect", fixed = TRUE)
  # Which reading of Step 2 produced these draws is part of the answer.
  expect_match(out, "rejuvenated as a\\s+whole")
})


test_that("threshold = 0 is Algorithm 2 on the whole parameter vector", {
  skip_on_cran()
  d <- weak_confounder_data(3, n = 150, n_weak = 3, n_noise = 5)
  f <- covariate_formulas(d)
  # A tolerance these draws do not already meet, so that a sweep really runs:
  # the equivalence is a claim about the sweep and says nothing if Lemma 1
  # skips it, which at the default tolerance is what this fixture does.
  args <- list(f$outcome, f$ps, d, seed = 3, mc = 600, bn = 200, thin = 2,
               chains = 2L,
               control = drbayes_control(n_steps = 50, tol = 1e-8),
               diagnostics = "none")

  sel <- do.call(fit_quietly, c(args, list(threshold = 0)))
  expect_length(sel$selected, sel$n_candidates)
  expect_gt(sel$smc$n_steps, 0L)

  utils::capture.output(
    pc <- suppressWarnings(do.call(
      drbayes_pc, c(args, list(outcome.prior = "horseshoe",
                               ps.prior = "horseshoe", verbose = FALSE)))),
    file = nullfile())

  # With nothing excluded the restricted moment condition is the unrestricted
  # one, so drbayes_select() is driving the sweep of R/tilting.R with the same
  # arguments drbayes_pc() drives it with, and must agree draw for draw.
  expect_equal(sel$pc, pc$pc)
  expect_equal(sel$g.comp, pc$g.comp)
  expect_equal(sel$smc$lambda, pc$smc$lambda)

  # And again through select_smc(), the separate loop the literal reading of
  # Step 2 uses. With every coefficient in the updated block it has nothing to
  # freeze, so it has to reproduce tilt_smc() step for step; nothing else in
  # the suite holds the two loops together.
  literal <- do.call(fit_quietly, c(args, list(threshold = 0,
                                               rejuvenate = "selected")))
  expect_equal(literal$pc, pc$pc)
  expect_equal(literal$smc$lambda, pc$smc$lambda)
  # Nothing was carried, so there is no ancestry for the count to describe,
  # and print() says that rather than reporting NA distinct draws.
  expect_true(is.na(literal$smc$n_ancestors_cumulative))
  expect_match(gsub("\\s+", " ",
                    paste(utils::capture.output(print(literal)),
                          collapse = " ")),
               "nothing outside that block to carry")
})


test_that("the tilting settings on the control object are honoured", {
  skip_on_cran()
  d <- weak_confounder_data(6, n = 200, n_weak = 3, n_noise = 5)
  f <- covariate_formulas(d)
  fit_with <- function(control, ...) {
    fit_quietly(f$outcome, f$ps, d, threshold = 0.05, seed = 6, mc = 600,
                bn = 200, thin = 2, chains = 2L, control = control,
                diagnostics = "none", ...)
  }

  plain <- fit_with(drbayes_control(n_steps = 50))
  # Sample pruning and the subclassification moment condition are settings of
  # the sweep, so they have to reach it here as they do through drbayes_pc().
  # prune.q is read by the quantile rule; the default weight rule reads
  # prune.w instead and discards nothing at these tilting parameters.
  pruning <- drbayes_control(n_steps = 50, pruning = TRUE,
                             prune.rule = "quantile", prune.q = 0.5)
  pruned <- fit_with(pruning)
  strata <- fit_with(drbayes_control(n_steps = 50, moment = "subclass",
                                     n_subclass = 4L))

  expect_false(isTRUE(all.equal(plain$pc, pruned$pc)))
  expect_false(isTRUE(all.equal(plain$pc, strata$pc)))
  # Half the particles are given zero weight before every resampling step, so
  # at most half of them can survive one however the draws fall.
  survivors <- length(plain$pc) - floor(0.5 * length(plain$pc))
  expect_lte(pruned$smc$n_ancestors_min, survivors)
  expect_gt(plain$smc$n_ancestors_min, pruned$smc$n_ancestors_min)

  # The literal reading of Step 2 runs a sweep of its own, so pruning has to
  # reach that one on the same terms rather than by another rule or not at all.
  literal <- fit_with(pruning, rejuvenate = "selected")
  expect_lte(literal$smc$n_ancestors_min, survivors)
  expect_false(isTRUE(all.equal(literal$pc,
                                fit_with(drbayes_control(n_steps = 50),
                                         rejuvenate = "selected")$pc)))
})


test_that("a threshold no covariate reaches names the largest one", {
  skip_on_cran()
  d <- weak_confounder_data(4, n = 150, n_weak = 3, n_noise = 5)
  f <- covariate_formulas(d)

  expect_error(
    fit_quietly(f$outcome, f$ps, d, threshold = 50, seed = 4, mc = 600,
                bn = 200, thin = 2, chains = 2L,
                control = drbayes_control(n_steps = 50),
                diagnostics = "none"),
    "Lower threshold")
})


# Everything below refuses the call before Step 1 is sampled, so it costs
# milliseconds and runs on CRAN as well.
test_that("arguments that cannot be honoured say what to do about it", {
  d <- weak_confounder_data(4, n = 150, n_weak = 3, n_noise = 5)
  f <- covariate_formulas(d)

  expect_error(drbayes_select(f$outcome, f$ps, d, threshold = -1),
               "must be a single non-negative number")
  expect_error(drbayes_select(f$outcome, f$ps, d, threshold = c(0.1, 0.2)),
               "must be a single non-negative number")
  expect_error(drbayes_select(f$outcome, A ~ 1, d),
               "no covariates, so there is nothing to select")
  expect_error(drbayes_select(Y ~ A, f$ps, d),
               "nothing to shrink")
  expect_error(drbayes_select(f$outcome, f$ps, d, rejuvenate = "some"),
               "should be one of")
  expect_error(drbayes_select(f$outcome, f$ps, d, method = "mcmc"),
               "should be one of")
  # A sampler drawing on a scale the moment condition would not evaluate it on
  # is refused before it is run, not silently misread.
  expect_error(drbayes_select(f$outcome, f$ps, d, ps.model = bayes_probit,
                              verbose = FALSE),
               "ps.link = \"probit\"")
})


test_that("selection is on the propensity score, not the outcome model", {
  skip_on_cran()
  set.seed(11)
  n <- 250
  d <- data.frame(Ptreat = rnorm(n), Pout = rnorm(n), Pnone = rnorm(n))
  # Ptreat drives treatment only, Pout drives the outcome only.
  d$A <- rbinom(n, 1, plogis(1.5 * d$Ptreat))
  d$Y <- 2 * d$A + 3 * d$Pout + rnorm(n)

  fit <- fit_quietly(Y ~ A + Ptreat + Pout + Pnone,
                     A ~ Ptreat + Pout + Pnone, d, threshold = 0.25,
                     seed = 11, mc = 800, bn = 200, thin = 2, chains = 2L,
                     control = drbayes_control(n_steps = 50),
                     diagnostics = "none")

  # The point of section 7.2: a covariate is kept for predicting TREATMENT.
  # A large outcome coefficient earns Pout nothing here.
  expect_true("Ptreat" %in% fit$selected)
  expect_false("Pout" %in% fit$selected)
  expect_false("Pnone" %in% fit$selected)
})


test_that("the treatment is never a candidate and is never frozen", {
  skip_on_cran()
  set.seed(12)
  n <- 200
  d <- data.frame(X1 = rnorm(n), X2 = rnorm(n), X3 = rnorm(n))
  d$A <- rbinom(n, 1, plogis(1.2 * d$X1))
  d$Y <- 2 * d$A + d$X1 + 0.7 * d$A * d$X2 + rnorm(n)

  fit <- fit_quietly(Y ~ A + X1 + X2 + X3 + A:X2, A ~ X1 + X2 + X3,
                     d, threshold = 0.4, seed = 12, mc = 800, bn = 200,
                     thin = 2, chains = 2L, rejuvenate = "selected",
                     control = drbayes_control(n_steps = 50,
                                               keep_particles = TRUE),
                     diagnostics = "none")

  # The treatment is the response of ps.formula, so it is not among the
  # candidates at all.
  expect_false("A" %in% names(fit$posterior_means))
  expect_identical(fit$n_candidates, 3L)
  expect_false("X2" %in% fit$selected)

  # A and A:X2 must stay in the block the sweep updates, even though X2 was
  # not selected: freezing them would leave the tilting no way to move the
  # treatment effect however far the moment condition is from zero.
  md <- prepare_model_data(Y ~ A + X1 + X2 + X3 + A:X2, A ~ X1 + X2 + X3,
                           d, "na.omit")
  treatment <- colnames(md$Z.lm)[md$treatment_columns + 1L]
  expect_setequal(treatment, c("A", "A:X2"))

  updated <- updated_columns(fit, md)
  expect_true(all(md$treatment_columns + 1L %in% updated))
  expect_false(match("X2", colnames(md$Z.lm)) %in% updated)

  # And the particles agree: both treatment coefficients were refreshed.
  distinct <- distinct_per_column(fit$particles$outcome)
  expect_equal(unname(distinct[treatment]),
               rep(nrow(fit$particles$outcome), 2L))
})


test_that("a selected covariate is matched through the terms, not the name", {
  skip_on_cran()
  set.seed(14)
  n <- 400
  d <- data.frame(X1 = runif(n, 1, 5), X2 = rnorm(n), X3 = rnorm(n))
  d$A <- rbinom(n, 1, plogis(1.5 * (log(d$X1) - 1) + 0.8 * d$X2))
  d$Y <- 2 * d$A + 3 * log(d$X1) + 0.5 * d$X2 + 0.5 * d$X3 + rnorm(n)

  # The confounder is written X1 in the treatment model and log(X1) in the
  # outcome model, which is the case that matching design matrix column names
  # gets wrong: it freezes the coefficient Algorithm 4 exists to revive.
  fit <- fit_quietly(Y ~ A + log(X1) + X2 + X3, A ~ X1 + X2 + X3, d,
                     threshold = 0.2, seed = 14, mc = 800, bn = 200,
                     thin = 2, chains = 2L, rejuvenate = "selected",
                     control = drbayes_control(n_steps = 100,
                                               keep_particles = TRUE),
                     diagnostics = "none")
  expect_setequal(fit$selected, c("X1", "X2"))

  # Read off the particles rather than from the mapping itself: a refreshed
  # coefficient is distinct at every particle, a frozen one is not.
  distinct <- distinct_per_column(fit$particles$outcome)
  particles <- nrow(fit$particles$outcome)
  expect_identical(unname(distinct[["log(X1)"]]), particles)
  expect_identical(unname(distinct[["X2"]]), particles)
  # X3 was not selected, so it is the one the sweep leaves alone.
  expect_lt(distinct[["X3"]], particles)
  expect_identical(unname(distinct[["X3"]]), fit$smc$n_ancestors_cumulative)

  # The same mapping through poly() and through a factor.
  md <- prepare_model_data(Y ~ A + poly(X1, 2) + X2, A ~ X1 + X2, d,
                           "na.omit")
  reached <- selection_columns(2L, md)$columns
  expect_setequal(colnames(md$Z.lm)[reached],
                  c("poly(X1, 2)1", "poly(X1, 2)2"))

  d$f <- factor(rep(letters[1:3], length.out = n))
  md <- prepare_model_data(Y ~ A + f + X2, A ~ f + X2, d, "na.omit")
  # One dummy column of a factor is selected; the term it belongs to owns
  # every dummy column of that factor in the outcome model.
  expect_setequal(colnames(md$Z.lm)[selection_columns(2L, md)$columns],
                  c("fb", "fc"))
})


test_that("a selected covariate the outcome model lacks is reported", {
  skip_on_cran()
  set.seed(15)
  n <- 300
  d <- data.frame(V1 = rnorm(n), V2 = rnorm(n), Q1 = rnorm(n), Q2 = rnorm(n))
  d$A <- rbinom(n, 1, plogis(1.5 * d$V1))
  d$Y <- 2 * d$A + 2 * d$Q1 + rnorm(n)

  # The two formulas share no covariate at all, so the coupling has no outcome
  # coefficient of a confounder to revive and, read literally, the sweep can
  # move nothing but the intercept and the treatment coefficient.
  utils::capture.output(
    expect_warning(
      drbayes_select(Y ~ A + Q1 + Q2, A ~ V1 + V2, d, threshold = 0.01,
                     seed = 15, mc = 600, bn = 200, thin = 2, chains = 2L,
                     rejuvenate = "selected",
                     control = drbayes_control(n_steps = 40),
                     diagnostics = "none", verbose = FALSE),
      "contribute no column to the outcome model"),
    file = nullfile())
})


test_that("the carried block is rejuvenated by default, frozen on request", {
  skip_on_cran()
  set.seed(41)
  n <- 300
  p <- 8
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", seq_len(p))
  d <- as.data.frame(X)
  d$A <- rbinom(n, 1, plogis(0.8 * d$X1 + 0.6 * d$X2))
  # A binary outcome is where freezing does damage: the risk difference is
  # non-linear in beta, so the coefficients outside the selected set reach the
  # credible interval.
  d$Y <- rbinom(n, 1, plogis(0.8 * d$A + 0.4 * d$X1 + 1.0 * d$X4))
  covariates <- paste(colnames(X), collapse = " + ")
  outcome <- stats::as.formula(paste("Y ~ A +", covariates))
  ps      <- stats::as.formula(paste("A ~", covariates))

  # The warnings are collected rather than suppressed: whether the collapse is
  # reported is part of what is being tested.
  fit_with <- function(rejuvenate) {
    seen <- character(0)
    utils::capture.output(
      fit <- withCallingHandlers(
        drbayes_select(outcome, ps, d, threshold = 0.05, seed = 41, mc = 700,
                       bn = 200, thin = 2, chains = 2L,
                       rejuvenate = rejuvenate,
                       control = drbayes_control(n_steps = 150,
                                                 keep_particles = TRUE),
                       diagnostics = "none", verbose = FALSE),
        warning = function(cond) {
          seen <<- c(seen, conditionMessage(cond))
          invokeRestart("muffleWarning")
        }),
      file = nullfile())
    list(fit = fit, warnings = seen)
  }
  refreshed <- fit_with("all")
  frozen    <- fit_with("selected")

  md      <- prepare_model_data(outcome, ps, d, "na.omit")
  carried <- setdiff(seq_len(ncol(md$Z.lm)),
                     updated_columns(refreshed$fit, md))
  expect_gt(length(carried), 0L)

  block <- function(fit) {
    apply(fit$particles$outcome[, carried, drop = FALSE], 1L, paste,
          collapse = "|")
  }
  particles <- nrow(refreshed$fit$particles$outcome)

  # The default leaves every particle its own draw of the carried block, and
  # has nothing to warn about.
  expect_identical(length(unique(block(refreshed$fit))), particles)
  expect_true(is.na(refreshed$fit$smc$n_ancestors_cumulative))
  expect_false(any(grepl("point mass", refreshed$warnings)))

  # The literal reading does not, and n_ancestors_cumulative counts exactly
  # how far it has collapsed.
  expect_identical(length(unique(block(frozen$fit))),
                   frozen$fit$smc$n_ancestors_cumulative)
  expect_lt(frozen$fit$smc$n_ancestors_cumulative, 0.1 * particles)
  # The per-step count is the statistic that cannot see this: most particles
  # survive any single resampling step while the block is a handful of points.
  expect_gt(frozen$fit$smc$n_ancestors_min,
            10 * frozen$fit$smc$n_ancestors_cumulative)
  # So the collapse is reported, and the fit is not called converged.
  expect_true(any(grepl("point mass", frozen$warnings)))
  expect_false(frozen$fit$smc$converged)

  # The propensity score coefficients outside the selected set ride the same
  # particle and collapse with it, so $particles$ps is as much a point mass as
  # $particles$outcome and the warning has to say so.
  ps_columns <- 1L + match(frozen$fit$selected, colnames(md$Z.ps)[-1L])
  carried_ps <- setdiff(seq_len(ncol(md$Z.ps)), c(1L, ps_columns))
  expect_gt(length(carried_ps), 0L)
  expect_equal(max(distinct_per_column(
                     frozen$fit$particles$ps[, carried_ps, drop = FALSE])),
               frozen$fit$smc$n_ancestors_cumulative)
  expect_true(any(grepl("propensity score coefficient", frozen$warnings)))

  # Read back with the line breaks flattened, since the text is wrapped to
  # whatever the console width happens to be.
  printed <- function(fit) {
    gsub("\\s+", " ", paste(utils::capture.output(print(fit)), collapse = " "))
  }

  # The convergence line print() inherits reads every fit that is not
  # converged as one that missed the moment condition, so a fit stopped by the
  # collapse alone has to say which of the two it was. Which one stops a given
  # sweep is a property of the data, so the moment condition is put inside the
  # tolerance here rather than hunted for with a fixture.
  met <- frozen$fit
  met$smc$B_mean <- 0.5 * met$smc$tol
  expect_match(printed(met), "THE MOMENT CONDITION WAS NOT MET")
  expect_match(printed(met), "not the moment condition")
  # A fit that really did miss it keeps the inherited explanation.
  expect_gt(abs(frozen$fit$smc$B_mean), frozen$fit$smc$tol)
  expect_no_match(printed(frozen$fit), "not the moment condition")
})


test_that("frozen coefficients are carried whole, not dropped or mixed", {
  skip_on_cran()
  d <- weak_confounder_data(5, n = 200, n_weak = 3, n_noise = 7)
  f <- covariate_formulas(d)

  fit <- fit_quietly(f$outcome, f$ps, d, threshold = 0.3, seed = 5,
                     mc = 800, bn = 200, thin = 2, chains = 2L,
                     rejuvenate = "selected",
                     control = drbayes_control(n_steps = 60,
                                               keep_particles = TRUE),
                     diagnostics = "none")
  expect_lt(length(fit$selected), fit$n_candidates)

  md      <- prepare_model_data(f$outcome, f$ps, d, "na.omit")
  updated <- updated_columns(fit, md)
  frozen  <- setdiff(seq_len(ncol(md$Z.lm)), updated)
  expect_gt(length(frozen), 0L)

  # The Step 1 draws the sweep started from: drbayes_select() sets the seed
  # once and then samples the outcome model before the propensity score one.
  set.seed(5)
  step1 <- bayes_lm_hs(md$Y, md$X.lm, mc = 800, chains = 2L, init = NULL,
                       unshrunk = md$treatment_columns)
  start <- matrix(step1[seq(201L, 800L, by = 2L), , , drop = FALSE],
                  ncol = ncol(md$Z.lm),
                  dimnames = list(NULL, colnames(md$Z.lm)))

  # Every frozen coefficient of a tilted particle must be an exact copy of the
  # corresponding coefficient of ONE Step 1 draw: the sweep resamples the whole
  # particle and then leaves this block alone. A per-coefficient copy would
  # leave the particle incoherent, and dropping the block would lose it.
  ancestor <- function(x) apply(x[, frozen, drop = FALSE], 1L, paste,
                                collapse = "|")
  expect_true(all(ancestor(fit$particles$outcome) %in% ancestor(start)))

  index <- match(ancestor(fit$particles$outcome), ancestor(start))
  expect_false(anyNA(index))
  # The selected block did move, so the two blocks really were treated
  # differently rather than the sweep having done nothing at all.
  expect_gt(max(abs(fit$particles$outcome[, updated, drop = FALSE] -
                      start[index, updated, drop = FALSE])), 0)
})


test_that("the g-computation uses the full beta, not the selected block", {
  skip_on_cran()
  set.seed(13)
  n <- 300
  d <- data.frame(Ptreat = rnorm(n), Pout = rnorm(n))
  d$A <- rbinom(n, 1, plogis(1.5 * d$Ptreat))
  # A binary outcome makes the estimand non-linear in beta, so a coefficient
  # outside the selected set genuinely changes the risk difference.
  d$Y <- rbinom(n, 1, plogis(0.8 * d$A + 2 * d$Pout))

  fit <- fit_quietly(Y ~ A + Ptreat + Pout, A ~ Ptreat + Pout,
                     d, threshold = 0.3, seed = 13, mc = 800, bn = 200,
                     thin = 2, chains = 2L,
                     control = drbayes_control(n_steps = 50,
                                               keep_particles = TRUE),
                     diagnostics = "none")
  expect_false("Pout" %in% fit$selected)

  md <- prepare_model_data(Y ~ A + Ptreat + Pout, A ~ Ptreat + Pout, d,
                           "na.omit")
  risk_difference <- function(betas) {
    rowMeans(stats::plogis(tcrossprod(betas, md$Z.lm1))) -
      rowMeans(stats::plogis(tcrossprod(betas, md$Z.lm0)))
  }

  expect_equal(risk_difference(fit$particles$outcome), fit$pc)

  # Dropping the carried coefficients would answer a different question.
  without <- fit$particles$outcome
  without[, "Pout"] <- 0
  expect_gt(abs(mean(risk_difference(without)) - mean(fit$pc)), 1e-6)
})


test_that("Step 2 can be run by either algorithm, and reports which", {
  skip_on_cran()
  set.seed(9)
  n <- 200
  d <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  d$A <- rbinom(n, 1, plogis(0.8 * d$X1))
  d$Y <- 2 * d$A + d$X1 + rnorm(n)
  args <- list(Y ~ A + X1 + X2, A ~ X1 + X2, d, threshold = 0.05, seed = 9,
               mc = 400, bn = 100, thin = 2, chains = 2L,
               diagnostics = "none")

  # Algorithm 1, which the paper offers alongside Algorithm 2 in Step 2. It
  # reweights and resamples once, so it reports an effective sample size and
  # never sweeps.
  one <- do.call(fit_quietly,
                 c(args, list(method = "is",
                              control = drbayes_control(n_steps = 20))))
  expect_identical(one$method, "is")
  # Algorithm 1 applies no kernel, so the fit does not claim a block was
  # refreshed even though rejuvenate has a default value.
  expect_true(is.na(one$rejuvenate))
  expect_false(is.na(one$smc$ess))
  expect_true(is.na(one$smc$n_ancestors_min))
  expect_identical(one$smc$max_steps, drbayes_control()$newton_steps)
  expect_identical(one$smc$n_ancestors_cumulative, one$smc$n_distinct)
  expect_match(paste(utils::capture.output(print(one)), collapse = " "),
               "resampled them once")

  two <- do.call(fit_quietly,
                 c(args, list(control = drbayes_control(n_steps = 20))))
  expect_true(is.na(two$smc$ess))
  expect_false(is.na(two$smc$n_ancestors_min))
  expect_false(isTRUE(all.equal(one$pc, two$pc)))

  # Lemma 1: a tolerance the untilted draws already meet skips the sweep, and
  # then every draw is still in play. run_tilting() documents ess and the
  # distinct counts as the number of draws in exactly that case.
  skipped <- do.call(fit_quietly,
                     c(args, list(control = drbayes_control(n_steps = 20,
                                                            tol = 100))))
  expect_identical(skipped$smc$n_steps, 0L)
  expect_identical(skipped$smc$ess, as.numeric(length(skipped$pc)))
  expect_identical(skipped$smc$n_distinct, length(skipped$pc))
  expect_identical(skipped$smc$n_ancestors_cumulative, length(skipped$pc))
  expect_equal(skipped$pc, skipped$g.comp)
})


test_that("Step 1 takes either prior and a sampler of the caller's own", {
  skip_on_cran()
  set.seed(9)
  n <- 200
  d <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  d$A <- rbinom(n, 1, plogis(0.8 * d$X1))
  d$Y <- 2 * d$A + d$X1 + rnorm(n)
  args <- list(Y ~ A + X1 + X2, A ~ X1 + X2, d, threshold = 0.05, seed = 9,
               mc = 400, bn = 100, thin = 2, chains = 2L,
               control = drbayes_control(n_steps = 20), diagnostics = "none")

  # "Shrinkage priors, such as the horseshoe prior": the horseshoe is the
  # default, not the only option.
  shrunk <- do.call(fit_quietly, args)
  expect_identical(shrunk$outcome.prior, "horseshoe")
  flat <- do.call(fit_quietly, c(args, list(outcome.prior = "normal",
                                            ps.prior = "normal")))
  expect_identical(flat$outcome.prior, "normal")
  expect_identical(flat$ps.prior, "normal")
  expect_false(isTRUE(all.equal(shrunk$pc, flat$pc)))

  # A sampler handed in directly, which is how Step 1 draws from a backend the
  # package does not provide. bayes_lm here only has to reach the sweep.
  supplied <- do.call(fit_quietly, c(args, list(outcome.model = bayes_lm)))
  expect_equal(supplied$pc, flat$pc)
  expect_identical(supplied$outcome.prior, "supplied")
})


test_that("a positivity failure names the model that overflowed", {
  skip_on_cran()
  set.seed(42)
  n <- 1500
  d <- data.frame(X1 = rnorm(n), X2 = rnorm(n), Sep = rnorm(n))
  # Sep separates the treatment groups perfectly, so the propensity score
  # reaches 0 and 1 and the moment condition is not finite.
  d$A <- as.integer(d$Sep > 0)
  d$Y <- rnorm(n, 2 * d$A + d$X1)

  # X2 is not selected, so the score that overflowed is the one on Sep and X1
  # alone. Naming the formula the user wrote would send them to a different
  # model from the one that failed.
  expect_error(
    fit_quietly(Y ~ A + X1 + X2 + Sep, A ~ X1 + X2 + Sep, d,
                threshold = 0.2, seed = 42, mc = 400, bn = 100, thin = 2,
                chains = 2L, control = drbayes_control(n_steps = 20),
                diagnostics = "none"),
    "glm(A ~ Sep, binomial", fixed = TRUE)
})


test_that("the Step 1 convergence check reports how far off it is", {
  skip_on_cran()
  set.seed(9)
  n <- 200
  d <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  d$A <- rbinom(n, 1, plogis(0.8 * d$X1))
  d$Y <- 2 * d$A + d$X1 + rnorm(n)

  # Twenty draws per chain cannot support the diagnostics, and the message has
  # to say by how much, not only which parameters failed.
  expect_error(
    fit_quietly(Y ~ A + X1 + X2, A ~ X1 + X2, d, threshold = 0.05, seed = 9,
                mc = 60, bn = 20, thin = 2, chains = 2L,
                control = drbayes_control(n_steps = 20),
                diagnostics = "error"),
    "min ESS")
  expect_error(
    fit_quietly(Y ~ A + X1 + X2, A ~ X1 + X2, d, threshold = 0.05, seed = 9,
                mc = 60, bn = 20, thin = 2, chains = 2L,
                control = drbayes_control(n_steps = 20),
                diagnostics = "error"),
    "makes R-hat itself unreliable")
})


test_that("a weak confounder killed by shrinkage is kept and revived", {
  skip_on_cran()
  truth <- 2
  weak  <- paste0("W", 1:6)

  # Step 1 on its own: the premise of section 7.2. The horseshoe shrinks the
  # outcome coefficients of the weak confounders towards zero even though each
  # is 0.5 in the data generating process.
  d  <- weak_confounder_data(21, ate = truth)
  f  <- covariate_formulas(d)
  md <- prepare_model_data(f$outcome, f$ps, d, "na.omit")
  set.seed(200)
  step1 <- bayes_lm_hs(md$Y, md$X.lm, mc = 1200, chains = 2L,
                       unshrunk = md$treatment_columns)
  shrunk <- colMeans(apply(step1[401:1200, , , drop = FALSE], 3L, c))
  expect_lt(mean(abs(shrunk[weak])), 0.25)

  # One dataset is a coin flip, so the claim is made over several. Over seeds
  # 21 to 25 of this generator coupling is the closer of the two on four, and
  # loses on one, so the assertions below are made on the five together rather
  # than on any single fit. The margin is comfortable, but the individual
  # errors move with the sampler, so they are asserted and not written down.
  errors <- vapply(21:25, function(seed) {
    d   <- weak_confounder_data(seed, ate = truth)
    f   <- covariate_formulas(d)
    fit <- fit_quietly(f$outcome, f$ps, d, threshold = 0.01, seed = seed,
                       mc = 1200, bn = 400, thin = 2, chains = 2L,
                       control = drbayes_control(n_steps = 150),
                       diagnostics = "none")
    # Every weak confounder survives selection: it is kept for predicting
    # treatment, which is exactly what the shrunk outcome coefficient hid.
    expect_true(all(weak %in% fit$selected))
    # And the selection separates them from the noise rather than keeping
    # whatever the threshold happens to fall above: the smallest posterior
    # mean among the confounders is larger than the largest among the noise
    # covariates, by a factor of at least one and a half on these seeds.
    noise <- setdiff(names(fit$posterior_means), weak)
    expect_gt(min(abs(fit$posterior_means[weak])),
              max(abs(fit$posterior_means[noise])))
    c(g.comp = abs(mean(fit$g.comp) - truth),
      pc     = abs(mean(fit$pc) - truth))
  }, numeric(2))

  expect_gte(sum(errors["pc", ] < errors["g.comp", ]), 4L)
  expect_lt(mean(errors["pc", ]), 0.9 * mean(errors["g.comp", ]))
})


test_that("selection does not cost accuracy against unselected coupling", {
  skip_on_cran()
  truth <- 2
  d <- weak_confounder_data(21, ate = truth)
  f <- covariate_formulas(d)
  args <- list(f$outcome, f$ps, d, seed = 21, mc = 1200, bn = 400, thin = 2,
               chains = 2L, control = drbayes_control(n_steps = 150),
               diagnostics = "none")

  sel <- do.call(fit_quietly, c(args, list(threshold = 0.01)))
  utils::capture.output(
    nosel <- suppressWarnings(do.call(
      drbayes_pc, c(args, list(outcome.prior = "horseshoe",
                               ps.prior = "horseshoe", verbose = FALSE)))),
    file = nullfile())

  # Section 7.3 reports the selected fit as the more accurate one, 0.188
  # against 0.043 for a G-formula value of 0.229. On simulated data of this
  # kind the two are instead indistinguishable. What does reproduce is that
  # both are far better than the separate-shrinkage g-formula, so this test
  # asks only that restricting the moment condition to the selected block does
  # not undo the coupling.
  expect_lt(abs(mean(sel$pc) - truth), abs(mean(sel$g.comp) - truth))
  expect_lt(abs(mean(nosel$pc) - truth), abs(mean(nosel$g.comp) - truth))
  expect_lt(abs(mean(sel$pc) - mean(nosel$pc)), stats::sd(nosel$pc))
})


# The blocks below sample nothing and run in milliseconds, so they carry no
# skip_on_cran(): they are the coverage of Algorithm 4 that CRAN itself gets.

test_that("the threshold keeps the candidates whose mean reaches it", {
  alpha_bar <- c(W1 = 0.4, W2 = -0.02, Z1 = 0.004, Z2 = -0.5)
  expect_identical(select_confounders(alpha_bar, 0.01), c("W1", "W2", "Z2"))
  expect_identical(select_confounders(alpha_bar, 0), names(alpha_bar))
  # The comparison is |alpha_bar| >= threshold, so the boundary is selected.
  expect_identical(select_confounders(alpha_bar, 0.02), c("W1", "W2", "Z2"))
  # Selecting nothing is refused with the value that would select something.
  expect_error(select_confounders(alpha_bar, 1), "0.5")
  expect_error(select_confounders(alpha_bar, 1), "Z2")

  expect_identical(validate_threshold(0L), 0)
  expect_error(validate_threshold(NA_real_), "single non-negative number")
  expect_error(validate_threshold("0.1"), "single non-negative number")
})


test_that("a selected covariate reaches its outcome columns by term", {
  d <- data.frame(X1 = runif(20, 1, 5), X2 = rnorm(20), X3 = rnorm(20),
                  A = rep(0:1, 10), Y = rnorm(20))

  # Column 2 of the propensity score design matrix is X1, which the outcome
  # model writes as log(X1). Matching column names would find nothing.
  md    <- prepare_model_data(Y ~ A + log(X1) + X2 + X3, A ~ X1 + X2 + X3, d,
                              "na.omit")
  reach <- selection_columns(2L, md)
  expect_identical(colnames(md$Z.lm)[reach$columns], "log(X1)")
  expect_false(any(reach$unmatched))
  # The moment condition scores the treatment on the selected covariates only,
  # so that is the model a positivity failure inside the sweep belongs to.
  expect_identical(selected_formula_text(2L, md), "A ~ X1")

  md2 <- prepare_model_data(Y ~ A + X2, A ~ X1 + X2, d, "na.omit")
  expect_true(selection_columns(2L, md2)$unmatched)
  expect_warning(warn_unreached("X1", TRUE, "selected"),
                 "nothing but the intercept")
  expect_silent(warn_unreached(character(0), FALSE, "all"))
})


test_that("the carried block is held to the floor the sweep is held to", {
  # A fit with nothing carried has no count, and no collapse either.
  expect_false(carried_collapse(NA_integer_, 100L, 0.1))
  expect_false(carried_collapse(10L, 100L, 0.1))
  expect_true(carried_collapse(9L, 100L, 0.1))
})


test_that("the outcome link and the propensity score link are separate", {
  expect_identical(resolve_links(NULL, NULL, "logit", c(0, 1, 1))$family,
                   "binomial")
  expect_identical(resolve_links(NULL, NULL, "logit", c(0.4, 1, 2))$family,
                   "gaussian")
  # A probit propensity score alongside a logit outcome model: the moment
  # condition has to turn the draws into probabilities with its own link.
  both <- resolve_links("binomial", "logit", "probit", c(0, 1, 0))
  expect_identical(both$link, "logit")
  expect_equal(both$ps_inverse_link(0), 0.5)
  expect_equal(both$inverse_link(0), 0.5)

  expect_error(resolve_links("binomial", NULL, "logit", c(0.2, 1)),
               "needs a response that is 0 or 1")
  expect_error(resolve_links("gaussian", "logit", "logit", c(0.2, 1)),
               "not one of the links")

  expect_error(check_sampler_link(bayes_probit, "logit", "ps.model",
                                  "ps.link"),
               "Set ps.link = \"probit\"", fixed = TRUE)
  expect_null(check_sampler_link(NULL, "logit", "ps.model", "ps.link"))
  expect_null(check_sampler_link(bayes_logit, "logit", "ps.model", "ps.link"))

  expect_match(step1_description(NULL, "logit", "horseshoe"),
               "bayes_logit_hs (logit link, horseshoe prior)", fixed = TRUE)
  expect_match(step1_description(bayes_lm, "identity", "normal"),
               "the supplied sampler")
  # The treatment columns are exempted from shrinkage only by a sampler that
  # has somewhere to put them.
  expect_true(has_formal(bayes_lm_hs, "unshrunk"))
  expect_false(has_formal(bayes_lm, "unshrunk"))
})

# Helpers shared by the simulation scripts of Orihara, Momozaki and Sugasawa
# (2025), "Bayesian Doubly Robust Causal Inference via Posterior Coupling":
# the two estimators that are not part of the package, the five performance
# measures of section 5.2, and the seeding and parallelism the studies use.
#
# It is sourced by scripts/simulation_table1.R, scripts/simulation_table2.R,
# scripts/simulation_appendix_d1.R and scripts/simulation_appendix_d2.R, and is
# not meant to be run on its own.
#
# References
#   Bang, H., & Robins, J. M. (2005). Doubly robust estimation in missing data
#     and causal inference models. Biometrics, 61(4), 962-973.
#   Luo, Y., Graham, D. J., & McCoy, E. J. (2023). Semiparametric Bayesian
#     doubly robust causal estimation. Journal of Statistical Planning and
#     Inference, 225, 171-187.
#   Schennach, S. M. (2005). Bayesian exponentially tilted empirical
#     likelihood. Biometrika, 92(1), 31-46.
#   Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
#     Causal Inference via Posterior Coupling. arXiv:2506.04868.


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

#' Evaluate an expression under a fixed random seed
#'
#' The samplers that take a \code{seed} argument put the caller's random number
#' stream back the way they found it, so a fit is a pure function of its data
#' and its seed. \code{drbayes_bb()} and \code{luo_estimate()} take no seed
#' argument, and a whole replication is easier to reason about as one seeded
#' block than as a list of separately seeded fits, so this gives any expression
#' the same property.
#'
#' The generator is pinned as well as the seed. \code{set.seed()} keeps
#' whatever generator is in force, and a furrr worker runs under
#' L'Ecuyer-CMRG rather than the default Mersenne-Twister, so without this a
#' replication would draw a different stream in a worker than in this process
#' and the study would depend on how many workers ran it.
#'
#' @param seed Integer seed.
#' @param expr Expression to evaluate.
#'
#' @return The value of \code{expr}.
with_seed <- function(seed, expr) {
  # Restoring .Random.seed restores the generator too: its first element
  # encodes which one produced it.
  restore <- if (exists(".Random.seed", envir = globalenv(),
                        inherits = FALSE)) {
    saved <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    function() assign(".Random.seed", saved, envir = globalenv())
  } else {
    saved_kind <- RNGkind()
    function() {
      rm(".Random.seed", envir = globalenv())
      do.call(RNGkind, as.list(saved_kind))
    }
  }
  on.exit(restore(), add = TRUE)

  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  expr
}


#' One seed per replication and specification
#'
#' @param study_seed The single seed the script is run with.
#' @param n_replications Number of replications.
#' @param n_specifications Number of specifications each replication is run
#'   under.
#'
#' @return A \code{n_replications} by \code{n_specifications} integer matrix of
#'   seeds.
replication_seeds <- function(study_seed, n_replications,
                              n_specifications = 1L) {
  with_seed(study_seed,
            matrix(sample.int(.Machine$integer.max,
                              n_replications * n_specifications),
                   nrow = n_replications, ncol = n_specifications))
}


#' Apply a function over replications, in parallel when asked to
#'
#' Every replication is seeded by the caller, so the study returns the same
#' numbers however many workers run it and in whatever order they finish.
#' furrr is given a seed as well because it refuses to let a future touch the
#' random number generator without one, not because the result depends on it.
#'
#' Note that furrr and future are needed whatever \code{n_workers} is:
#' \code{run_parallel_simulation()} generates the datasets through them.
#'
#' @param indices Replication indices.
#' @param fun Function of one index.
#' @param n_workers Number of parallel workers. One runs the replications in
#'   this process.
#' @param seed Seed handed to furrr.
#'
#' @return A list, one element per index.
map_replications <- function(indices, fun, n_workers, seed) {
  if (n_workers <= 1L) {
    return(lapply(indices, fun))
  }
  old_plan <- future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(old_plan), add = TRUE)
  furrr::future_map(indices, fun,
                    .options = furrr::furrr_options(seed = seed))
}


# ---------------------------------------------------------------------------
# Estimators that are not part of the package
# ---------------------------------------------------------------------------

#' Summarise posterior draws the way Table 1 does
#'
#' @param draws Numeric vector of posterior draws of the treatment effect.
#' @param method Name of the estimator, as it should appear in the table.
#'
#' @return A one row data frame with the posterior mean and the bounds of an
#'   equal-tailed 95 percent credible interval.
summarise_draws <- function(draws, method) {
  bounds <- stats::quantile(draws, c(0.025, 0.975), names = FALSE)
  data.frame(method   = method,
             estimate = mean(draws),
             lower    = bounds[1],
             upper    = bounds[2],
             stringsAsFactors = FALSE)
}


#' Rows both models can be fitted to
#'
#' \code{lm()} and \code{glm()} drop incomplete cases one model at a time, and
#' the estimators below also read the treatment straight out of \code{data}.
#' A row that only one of the two models dropped would leave those vectors
#' different lengths, which R would recycle rather than refuse. Restricting
#' both models to the rows that are complete for both is what keeps them
#' aligned.
#'
#' @param outcome.formula,ps.formula The two model formulas.
#' @param data Data frame holding every variable in both formulas.
#'
#' @return \code{data} with the incomplete rows removed.
complete_rows <- function(outcome.formula, ps.formula, data) {
  variables <- union(all.vars(outcome.formula), all.vars(ps.formula))
  absent <- setdiff(variables, names(data))
  if (length(absent) > 0) {
    stop("These variables are named in the formulas but are not columns of ",
         "data: ", toString(absent), ". Add them, or take them out of the ",
         "formulas.", call. = FALSE)
  }

  complete <- stats::complete.cases(data[, variables, drop = FALSE])
  if (!all(complete)) {
    message(sum(!complete), " of ", nrow(data), " rows are missing a value ",
            "used by one of the two models and are dropped, so that both are ",
            "fitted to the same rows.")
  }
  data[complete, , drop = FALSE]
}


#' Augmented inverse probability weighted estimator
#'
#' The "DR" row of Table 1: the ordinary non-Bayesian doubly robust estimator
#' that is asymptotically equivalent to Bang and Robins (2005),
#' \deqn{n^{-1} \sum_i \left\{ m_1(X_i) - m_0(X_i) +
#'   \frac{A_i - e(X_i)}{e(X_i)(1 - e(X_i))} (Y_i - m_{A_i}(X_i)) \right\}.}
#' It is not part of the package, which is Bayesian throughout, so it lives
#' here with the rest of the comparison.
#'
#' The interval is a Wald interval built from the empirical variance of the
#' summand, which is the estimator's influence function when both models are
#' correctly specified. Nothing is subtracted for the estimation of the two
#' nuisance models, and that omission is not what shapes the coverage in
#' Table 1: over 2000 replications of the section 5.1 mechanism the average
#' standard error tracks the empirical one, 0.108 against 0.110 at n = 500 and
#' 0.062 against 0.062 at n = 1500, with coverage 94.5 and 95.4 percent. Where
#' the outcome model is misspecified the summand is no longer the influence
#' function and the interval is far too wide rather than too short, which is
#' why Table 1's Correct / Incorrect block reports a coverage of 100 on an
#' interval of average length 10.
#'
#' @param outcome.formula Formula for the outcome model.
#' @param ps.formula Formula for the propensity score model, with the treatment
#'   on the left.
#' @param data Data frame holding every variable in both formulas.
#'
#' @return A one row data frame in the same shape as \code{summarise_draws()}.
#'
#' @references
#' Bang, H., & Robins, J. M. (2005). Doubly robust estimation in missing data
#' and causal inference models. Biometrics, 61(4), 962-973.
aipw_estimate <- function(outcome.formula, ps.formula, data) {
  data <- complete_rows(outcome.formula, ps.formula, data)

  ps_fit      <- stats::glm(ps.formula, family = stats::binomial, data = data)
  outcome_fit <- stats::lm(outcome.formula, data = data)

  treatment <- all.vars(ps.formula)[1]
  A <- data[[treatment]]
  Y <- stats::model.response(stats::model.frame(outcome.formula, data))

  # Set the treatment in the DATA and predict again, rather than editing a
  # column of the design matrix: predict() rebuilds the model frame through
  # the fitted terms, so interactions with the treatment and any transformed
  # term are recomputed instead of being left at their observed values.
  counterfactual <- function(value) {
    modified <- data
    modified[[treatment]] <- value
    stats::predict(outcome_fit, newdata = modified)
  }

  e  <- stats::fitted(ps_fit)
  mu <- stats::fitted(outcome_fit)
  influence <- counterfactual(1) - counterfactual(0) +
    (A - e) / (e * (1 - e)) * (Y - mu)

  estimate <- mean(influence)
  halfwidth <- stats::qnorm(0.975) *
    stats::sd(influence) / sqrt(length(influence))

  data.frame(method   = "DR",
             estimate = estimate,
             lower    = estimate - halfwidth,
             upper    = estimate + halfwidth,
             stringsAsFactors = FALSE)
}


#' Bayesian empirical likelihood estimator of Luo et al. (2023)
#'
#' The "Luo" row of Table 1. Section 4.3 of the paper writes the method out.
#' The propensity score is fitted first and enters the outcome working model
#' \deqn{E[Y_i \mid A_i, X_i] = \tau A_i + h_1(X_i; \beta) +
#'   \phi e_i(\hat{\alpha})}
#' as one more covariate; the estimating equation (4.2) built from that model
#' has one component per parameter,
#' \deqn{U_i(\beta, \phi, \tau) = (A_i,\; \partial h_1 / \partial \beta,\;
#'   e_i)^\top (Y_i - \tau A_i - h_1(X_i; \beta) - \phi e_i);}
#' and the exponentially tilted empirical likelihood
#' \eqn{p_i \propto \exp(\lambda^\top U_i)} built from it replaces the
#' likelihood, with \eqn{\lambda} chosen at each parameter value so that the
#' weighted moments vanish. That inner problem is the convex minimisation of
#' \eqn{\log \sum_j \exp(\lambda^\top U_j)}, whose gradient is the moment
#' vector averaged under the weights \eqn{p_i} themselves, so a stationary
#' point is a \eqn{\lambda} at which the weighted moments vanish. It is solved
#' by BFGS from zero. The plain sum has the same stationary points, its
#' gradient being that average times a positive constant, but it overflows for
#' moderate \eqn{\lambda^\top U_j} where the logarithm, evaluated by factoring
#' out the largest exponent, does not.
#'
#' The posterior is sampled by a random walk Metropolis chain over
#' \eqn{(\tau, \beta, \phi)} under a flat prior, started at the solution of the
#' estimating equation and proposing from the least squares covariance of that
#' solution, scaled by the usual \eqn{2.4^2 / d}.
#'
#' This is a reimplementation from the description in section 4.3, not the
#' authors' own code, and it is here rather than in the package because the
#' package implements posterior coupling and nothing else.
#'
#' @param outcome.formula Formula for the outcome model. The treatment must
#'   appear as a main effect; its coefficient is the estimand.
#' @param ps.formula Formula for the propensity score model.
#' @param data Data frame holding every variable in both formulas.
#' @param n_draws Number of posterior draws to keep.
#' @param n_burn Number of leading draws to discard.
#'
#' @return Numeric vector of \code{n_draws} posterior draws of the treatment
#'   effect, carrying the Metropolis acceptance rate as the attribute
#'   \code{"acceptance"}.
#'
#' @references
#' Luo, Y., Graham, D. J., & McCoy, E. J. (2023). Semiparametric Bayesian
#' doubly robust causal estimation. Journal of Statistical Planning and
#' Inference, 225, 171-187.
#'
#' Schennach, S. M. (2005). Bayesian exponentially tilted empirical likelihood.
#' Biometrika, 92(1), 31-46.
luo_estimate <- function(outcome.formula, ps.formula, data,
                         n_draws = 4000L, n_burn = 1000L) {
  data <- complete_rows(outcome.formula, ps.formula, data)

  treatment <- all.vars(ps.formula)[1]
  ps_fit <- stats::glm(ps.formula, family = stats::binomial, data = data)
  e <- unname(stats::fitted(ps_fit))

  A <- data[[treatment]]
  Y <- stats::model.response(stats::model.frame(outcome.formula, data))
  X <- stats::model.matrix(outcome.formula, data)
  if (!treatment %in% colnames(X)) {
    stop("The treatment ", treatment, " is not a main effect of ",
         "outcome.formula, so this estimator has no coefficient to report. ",
         "Write it as ", deparse(outcome.formula[[2]]), " ~ ", treatment,
         " + ...", call. = FALSE)
  }
  # The design matrix of h_1, that is everything in the outcome model except
  # the treatment itself.
  H <- X[, setdiff(colnames(X), treatment), drop = FALSE]

  n <- length(Y)
  # The moment vector is the estimating equation (4.2) evaluated observation by
  # observation, so it has one column per parameter and the empirical
  # likelihood is just identified.
  moments <- function(theta) {
    residual <- Y - theta[1] * A - drop(H %*% theta[-c(1L, length(theta))]) -
      theta[length(theta)] * e
    cbind(A, H, e) * residual
  }

  log_posterior <- function(theta) {
    U <- moments(theta)
    objective <- function(lambda) {
      v <- drop(U %*% lambda)
      max(v) + log(sum(exp(v - max(v))))
    }
    gradient <- function(lambda) {
      v <- drop(U %*% lambda)
      w <- exp(v - max(v))
      drop(crossprod(U, w / sum(w)))
    }
    solved <- stats::optim(rep(0, ncol(U)), objective, gradient,
                           method = "BFGS",
                           control = list(maxit = 200L, reltol = 1e-10))
    v <- drop(U %*% solved$par)
    value <- sum(v) - n * (max(v) + log(sum(exp(v - max(v)))))
    if (is.finite(value)) value else -Inf
  }

  # The estimating equation is linear in the parameters, so its solution and
  # the covariance that scales the proposal both come from one least squares
  # fit.
  design <- cbind(A, H, e)
  start  <- stats::lm.fit(design, Y)
  theta  <- unname(start$coefficients)
  if (anyNA(theta)) {
    stop("The working model of section 4.3 is rank deficient: the treatment, ",
         "the outcome covariates and the fitted propensity score are ",
         "collinear. Drop a covariate from outcome.formula, or leave this ",
         "estimator out of the comparison.", call. = FALSE)
  }
  residual_variance <- sum(start$residuals^2) / (n - length(theta))
  proposal <- chol(residual_variance * chol2inv(qr.R(qr(design))) *
                     2.4^2 / length(theta))

  current  <- log_posterior(theta)
  kept     <- numeric(n_draws)
  accepted <- 0L

  for (iteration in seq_len(n_burn + n_draws)) {
    candidate <- theta + drop(stats::rnorm(length(theta)) %*% proposal)
    proposed  <- log_posterior(candidate)
    if (log(stats::runif(1)) < proposed - current) {
      theta    <- candidate
      current  <- proposed
      accepted <- accepted + 1L
    }
    if (iteration > n_burn) kept[iteration - n_burn] <- theta[1]
  }

  attr(kept, "acceptance") <- accepted / (n_burn + n_draws)
  kept
}


# ---------------------------------------------------------------------------
# Performance measures
# ---------------------------------------------------------------------------

#' The five performance measures of section 5.2
#'
#' ABias is the absolute difference between the average estimate and the truth,
#' ESE the standard deviation of the estimates across replications, RMSE the
#' root mean squared error about the truth, CP the percentage of intervals
#' containing it, and AvL their average length. All five are computed over
#' replications, not averaged over per-replication values. CP and AvL are NA
#' for an estimator reported without an interval, as in Table 2.
#'
#' @param results Data frame with one row per replication, block and method,
#'   holding \code{n}, \code{block}, \code{method}, \code{estimate},
#'   \code{lower} and \code{upper}.
#' @param true_ate The true average treatment effect.
#' @param block_order,method_order The order the rows are reported in.
#'
#' @return A data frame with one row per sample size, block and method.
performance_table <- function(results, true_ate, block_order, method_order) {
  parts <- split(results, list(results$n, results$block, results$method),
                 drop = TRUE)

  rows <- lapply(parts, function(part) {
    error <- part$estimate - true_ate
    data.frame(
      n      = part$n[1],
      block  = part$block[1],
      method = part$method[1],
      ABias  = abs(mean(error)),
      ESE    = stats::sd(part$estimate),
      RMSE   = sqrt(mean(error^2)),
      CP     = 100 * mean(part$lower <= true_ate & true_ate <= part$upper),
      AvL    = mean(part$upper - part$lower),
      stringsAsFactors = FALSE)
  })

  table <- do.call(rbind, rows)
  table <- table[order(table$n,
                       match(table$block, block_order),
                       match(table$method, method_order)), ]
  rownames(table) <- NULL
  table
}


#' Print a performance table the way the paper lays one out
#'
#' @param table Output of \code{performance_table()}.
#' @param measures Columns to show.
#'
#' @return \code{table}, invisibly.
print_performance_table <- function(table,
                                    measures = c("ABias", "ESE", "RMSE",
                                                 "CP", "AvL")) {
  digits <- c(ABias = 3L, ESE = 3L, RMSE = 3L, CP = 1L, AvL = 3L)

  for (n in unique(table$n)) {
    part <- table[table$n == n, ]
    cat("\nn = ", n, "\n", sep = "")

    shown <- data.frame(Block = part$block, Method = part$method,
                        check.names = FALSE, stringsAsFactors = FALSE)
    for (measure in measures) {
      shown[[measure]] <- formatC(part[[measure]], format = "f",
                                  digits = digits[[measure]])
    }

    # The paper leaves the block cell empty on the continuation rows, so that
    # the methods sharing a specification read as one block instead of several
    # unrelated rows.
    shown$Block[duplicated(part$block)] <- ""
    print(shown, row.names = FALSE)
  }

  invisible(table)
}

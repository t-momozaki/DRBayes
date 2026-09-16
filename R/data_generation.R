#' Generate Synthetic Dataset for Causal Inference Simulation
#'
#' Generates a synthetic dataset with confounders, treatment assignment, and outcomes
#' following a specific data generating process designed for testing causal inference methods.
#' Only the observed outcome is returned; the potential outcomes are used to build it and then discarded.
#'
#' @param nn Integer. Number of observations to generate. Default is 200.
#' @param pp Integer. Number of additional confounders to include beyond the base 4
#'   confounders. Default is 40. Can be 0 for only the 4 primary confounders.
#' @param seed Integer or NULL. Random seed for reproducibility. If NULL, no seed is set.
#'
#' @return A data.frame with \code{nn} rows and \code{6 + pp} columns:
#' \describe{
#'   \item{AA}{Binary treatment assignment (0 or 1)}
#'   \item{YY}{Observed outcome (continuous)}
#'   \item{WW1, WW2, WW3, WW4}{Primary confounding variables}
#'   \item{WW5, ..., WW(4+pp)}{Irrelevant covariates (only if pp > 0)}
#' }
#'
#' @details
#' The data generating process follows these steps:
#'
#' **Covariates:**
#' \itemize{
#'   \item 4 confounders: \eqn{X_1, X_2, X_3, X_4 \sim N(0, 1)}
#'   \item \code{pp} further covariates \eqn{X_{4+j} \sim N(\mu_j, 1)} with
#'     \eqn{\mu_j \sim \mathrm{Unif}(-2, 2)}. These enter neither the treatment
#'     nor the outcome model, so they are irrelevant covariates rather than
#'     confounders, matching the high-dimensional setting of the paper's
#'     Appendix D.1. Note that the paper draws \eqn{\mu_j \sim \mathrm{Unif}(-1, 1)},
#'     so this generator is twice as dispersed as the published one.
#' }
#'
#' **Treatment assignment:**
#' \deqn{P(A = 1 | X) = \mathrm{logit}^{-1}(X_1 - 0.5 X_2 + 0.25 X_3 + 0.1 X_4)}
#'
#' There is no intercept, so treatment prevalence is one half.
#'
#' **Potential outcomes:**
#' \deqn{Y^{(0)} = 100 + 27.4 X_1 + 13.7(X_2 + X_3 + X_4) + \epsilon}
#' \deqn{Y^{(1)} = 210 + 27.4 X_1 + 13.7(X_2 + X_3 + X_4) + \epsilon}
#'
#' where \eqn{\epsilon \sim N(0, 1)}. The SAME error appears in both potential
#' outcomes, as the paper's single \eqn{\epsilon_i} implies, so the unit level
#' treatment effect is exactly 110 for every observation.
#'
#' **Observed outcome:**
#' \deqn{Y = A \cdot Y^{(1)} + (1-A) \cdot Y^{(0)}}
#'
#' The true average treatment effect (ATE) is \eqn{E[Y^{(1)} - Y^{(0)}] = 210 - 100 = 110}.
#'
#' @examples
#' # Generate a dataset with only primary confounders
#' data_minimal <- generate_dataset(nn = 100, pp = 0, seed = 123)
#' head(data_minimal)
#' ncol(data_minimal)  # 6: AA, YY, WW1-WW4
#'
#' # Generate a dataset with additional, unrelated covariates
#' data_extended <- generate_dataset(nn = 100, pp = 20, seed = 123)
#' ncol(data_extended)  # 26: AA, YY, WW1-WW24
#'
#' # Check treatment assignment
#' table(data_minimal$AA)
#'
#' # The unit level treatment effect is exactly 110, so the difference in
#' # means below is that plus the confounding the covariates carry.
#' aggregate(YY ~ AA, data = data_minimal, FUN = mean)
#'
#' # Estimate propensity scores
#' ps_model <- glm(AA ~ WW1 + WW2 + WW3 + WW4, data = data_minimal,
#'                 family = binomial)
#' coef(ps_model)
#'
#' @seealso \code{\link{run_parallel_simulation}}, \code{\link{drbayes_pc}}
#'
#' @importFrom stats rnorm runif plogis
#'
#' @export
generate_dataset <- function(nn = 200, pp = 40, seed = NULL) {

  if (!is.numeric(nn) || length(nn) != 1L || nn <= 0 || nn != round(nn)) {
    stop("nn must be a single positive integer")
  }
  if (!is.numeric(pp) || length(pp) != 1L || pp < 0 || pp != round(pp)) {
    stop("pp must be a single non-negative integer")
  }

  # A seed here is LOCAL. Setting it and walking away would replace the
  # caller's stream, so a script that seeds itself at the top would lose
  # control of everything after this call.
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("seed must be NULL or a single non-missing number")
    }
    return(with_preserved_rng(seed,
                              generate_dataset(nn = nn, pp = pp, seed = NULL)))
  }

  # Four confounders, standard normal, as in section 5.1 of the paper.
  XX <- matrix(stats::rnorm(nn * 4), nn, 4)

  # Covariates that enter neither model. rnorm recycles the length-pp mean
  # vector over the nn * pp draws, so draw k has mean mu[((k - 1) %% pp) + 1];
  # byrow = TRUE then places draw (i - 1) * pp + j at [i, j], and since
  # (i - 1) * pp is a multiple of pp that draw carries mean mu[j]. Every column
  # therefore gets its own mean. Dropping the byrow = TRUE would leave every
  # column with roughly the same mean, so it is load bearing.
  noise <- if (pp > 0) {
    matrix(stats::rnorm(nn * pp, mean = stats::runif(pp, -2, 2)), nn, pp, TRUE)
  } else {
    NULL
  }

  # e(X) = expit(X1 - 0.5 X2 + 0.25 X3 + 0.1 X4), no intercept.
  AA_ln <- drop(tcrossprod(c(1, -0.5, 0.25, 0.1), XX))
  AA    <- as.integer(stats::runif(nn) < stats::plogis(AA_ln))

  # One error per unit, shared by both potential outcomes, as the paper's
  # Y_ai = 100 + 110a + 13.7(2 X1 + X2 + X3 + X4) + eps_i implies. The unit
  # level treatment effect is therefore exactly 110 for everyone.
  ee  <- stats::rnorm(nn, sd = 1)
  mu  <- 100 + 27.4 * XX[, 1] + 13.7 * (XX[, 2] + XX[, 3] + XX[, 4]) + ee
  YY0 <- mu
  YY1 <- mu + 110

  YY <- AA * YY1 + (1 - AA) * YY0

  DD <- data.frame(cbind(AA, YY, XX, noise))
  names(DD) <- c("AA", "YY", paste0("WW", seq_len(4 + pp)))
  DD
}


#' Run Parallel Simulation Study
#'
#' Executes multiple simulation runs in parallel, generating synthetic datasets
#' using \code{\link{generate_dataset}}. This function is designed for conducting
#' large-scale simulation studies to evaluate causal inference methods.
#'
#' @param num_simulations Integer. Number of simulation replications to run. Default is 1000.
#' @param nn Integer. Number of observations per dataset. Default is 200.
#' @param pp Integer. Number of irrelevant covariates per dataset. Default is 40.
#'   Can be 0 for only the 4 primary confounders.
#' @param seed Integer or NULL. Controls the whole study. Each replication is
#'   given its own L'Ecuyer-CMRG substream, derived in the parent process from
#'   this one seed, so the datasets are identical however many workers run.
#'   When NULL a seed is drawn from the caller's stream, which keeps a script
#'   that has already called set.seed() reproducible while letting two
#'   successive calls give different studies. The seed actually used is
#'   recorded in attr(result, "seed"), so any run can be replayed exactly.
#'   Either way the caller's own generator is left as it was found, in both its
#'   state and its kind: the parallel machinery switches to L'Ecuyer-CMRG while
#'   it derives the substreams and is switched back here, so a script that
#'   seeds itself keeps control of its own stream.
#' @param n_cores Integer or NULL. When given, a future plan with this many workers is
#'   installed for the duration of the call and the caller's plan is restored afterwards.
#'   When NULL the caller's existing plan is used and left alone, which is the recommended
#'   division of labour: the end user owns the plan, not the package.
#' @param verbose Logical. Whether to report progress through message(). Defaults to
#'   interactive(). Unlike printing, this can be silenced with suppressMessages().
#' @param strategy Character. Parallel processing strategy. Options are "multisession"
#'   (cross-platform) or "multicore" (Unix-like systems only). Default is "multisession".

#' @return A list of length \code{num_simulations}, where each element is a data.frame
#'   generated by \code{\link{generate_dataset}} with the specified parameters.
#'
#' @details
#' This function uses the \code{future} and \code{furrr} packages to enable efficient
#' parallel processing. It automatically handles:
#' \itemize{
#'   \item Parallel execution across multiple CPU cores
#'   \item Proper random number generation in parallel environments
#'   \item Progress reporting during execution
#'   \item Execution time measurement
#'   \item Automatic cleanup of parallel workers
#' }
#'
#' The function sets unique seeds for each simulation to ensure reproducibility
#' while maintaining independence between parallel workers.
#'
#' **Memory considerations:**
#' With large \code{num_simulations}, \code{nn}, or \code{pp}, the resulting list
#' can consume significant memory. Consider processing results in batches for
#' very large simulation studies.
#'
#' **Platform compatibility:**
#' \itemize{
#'   \item "multisession": Works on all platforms (Windows, macOS, Linux)
#'   \item "multicore": Unix-like systems only, generally faster but uses forking
#' }
#'
#' @examples
#' # With n_cores left at NULL the caller's future plan is used, which is a
#' # sequential one until the caller installs another.
#' results_minimal <- run_parallel_simulation(
#'   num_simulations = 10,
#'   nn = 100,
#'   pp = 0,  # Only primary confounders
#'   seed = 1
#' )
#'
#' length(results_minimal)  # 10
#' dim(results_minimal[[1]])  # 100 x 6: AA, YY, WW1-WW4
#' attr(results_minimal, "seed")  # replays the study exactly
#'
#' # Summary statistics across simulations
#' exposure_rates <- sapply(results_minimal, function(df) mean(df$AA))
#' outcome_means <- sapply(results_minimal, function(df) mean(df$YY))
#' mean(exposure_rates)  # Average treatment assignment rate
#' mean(outcome_means)   # Average outcome across simulations
#'
#' \donttest{
#' # Two worker processes, which the function starts and shuts down again.
#' # The datasets are the same as above: the substreams are derived from the
#' # seed in the parent, not from the workers.
#' results_parallel <- run_parallel_simulation(
#'   num_simulations = 10,
#'   nn = 100,
#'   pp = 0,
#'   seed = 1,
#'   n_cores = 2,
#'   strategy = "multisession"
#' )
#' identical(results_parallel[[1]], results_minimal[[1]])
#' }
#'
#' @seealso \code{\link{generate_dataset}}, \code{\link{drbayes_pc}}
#'
#' @importFrom future plan
#' @importFrom furrr future_map furrr_options
#'
#' @export
run_parallel_simulation <- function(num_simulations = 1000,
                                    nn = 200,
                                    pp = 40,
                                    seed = NULL,
                                    n_cores = NULL,
                                    strategy = c("multisession", "multicore",
                                                 "sequential"),
                                    verbose = interactive()) {

  if (!is.numeric(num_simulations) || length(num_simulations) != 1L ||
      num_simulations <= 0 || num_simulations != round(num_simulations)) {
    stop("num_simulations must be a single positive integer")
  }
  if (!is.numeric(nn) || length(nn) != 1L || nn <= 0 || nn != round(nn)) {
    stop("nn must be a single positive integer")
  }
  if (!is.numeric(pp) || length(pp) != 1L || pp < 0 || pp != round(pp)) {
    stop("pp must be a single non-negative integer")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE")
  }

  # Resolve the seed BEFORE saving the caller's state. Capturing the state first
  # and restoring it on exit would undo this draw, and every seed = NULL call
  # would then return the same study.
  if (is.null(seed)) {
    # Draw one integer from the caller's stream, so a script that has already
    # called set.seed() is reproducible and two successive calls differ.
    seed <- sample.int(.Machine$integer.max, 1L)
  } else {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("seed must be NULL or a single non-missing integer")
    }
    seed <- as.integer(seed)
  }

  # furrr switches the generator to L'Ecuyer-CMRG to derive its substreams, and
  # puts back whatever .Random.seed it found. In a session that has never drawn
  # a random number there is nothing to put back, so the switch would otherwise
  # be permanent and every later set.seed() in the caller's script would produce
  # a different stream than it did before this call. CRAN policy forbids a
  # package leaving the kind or the seed of the user's generator changed, so
  # both are restored here, including the case where the caller had no
  # .Random.seed at all. The effect of this function on the generator is then
  # exactly the one draw the seed = NULL branch above takes, and none at all
  # when a seed is given.
  old_kind <- RNGkind()
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
            add = TRUE)
  }
  # Setting the kind re-seeds the generator, so it has to be put back before the
  # seed is rather than after it. on.exit() runs its expressions in the order
  # they were registered, and after = FALSE puts this one first.
  on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]), add = TRUE,
          after = FALSE)

  # Only touch the plan when asked to, and put back what was there. Resetting
  # to sequential would shut down worker pools the caller set up for the rest
  # of their script, which shows up as a silent loss of parallelism.
  if (!is.null(n_cores)) {
    if (!is.numeric(n_cores) || length(n_cores) != 1L || n_cores < 1 ||
        n_cores != round(n_cores)) {
      stop("n_cores must be NULL or a single positive integer")
    }
    strategy <- match.arg(strategy)
    n_cores  <- max(1L, min(as.integer(n_cores), as.integer(num_simulations)))
    old_plan <- future::plan(strategy, workers = n_cores)
    on.exit(future::plan(old_plan), add = TRUE)
  }

  if (verbose) {
    message("Generating ", num_simulations, " datasets of ", nn,
            " observations with ", 4 + pp, " covariates (",
            pp, " of them unrelated to treatment and outcome), seed ", seed, ".")
  }

  # furrr assigns each element its own L'Ecuyer-CMRG substream, derived in the
  # parent from this one seed. The substreams are provably disjoint and the
  # assignment does not depend on how many workers run, so the study is
  # reproducible on any machine and at any level of parallelism.
  results <- furrr::future_map(
    .x = seq_len(num_simulations),
    .f = function(i) generate_dataset(nn = nn, pp = pp),
    .options  = furrr::furrr_options(seed = seed, scheduling = 1),
    .progress = verbose
  )

  # Record the seed so a run started without one can be replayed exactly.
  attr(results, "seed") <- seed
  results
}


#' Convert Simulation Results to Array Format
#'
#' Converts a list of simulation results (data.frames) to a 3-dimensional array
#' format for easier analysis and manipulation. This is a utility function for
#' post-processing simulation results from \code{\link{run_parallel_simulation}}.
#'
#' @param simulation_results List of data.frames, typically the output from
#'   \code{\link{run_parallel_simulation}}.
#'
#' @return A 3-dimensional array with dimensions:
#' \describe{
#'   \item{Dimension 1}{Observations (rows of each data.frame)}
#'   \item{Dimension 2}{Variables (columns of each data.frame)}
#'   \item{Dimension 3}{Simulations (list elements)}
#' }
#'
#' The array includes appropriate dimension names for easy indexing.
#'
#' @details
#' This function is useful when you need to:
#' \itemize{
#'   \item Perform array-based operations across simulations
#'   \item Calculate statistics across the simulation dimension
#'   \item Export results to other software that expects array format
#'   \item Use vectorized operations for faster computation
#' }
#'
#' **Memory warning:** The resulting array can be very large for many simulations
#' or high-dimensional datasets. Consider the memory requirements before conversion.
#'
#' @examples
#' # Generate simulation results with only primary confounders
#' results <- run_parallel_simulation(num_simulations = 5, nn = 100, pp = 0,
#'                                    seed = 1)
#'
#' # Convert to array format
#' result_array <- convert_to_array(results)
#'
#' # Dimensions: observations, variables, simulations
#' dim(result_array)  # 100 x 6 x 5 for pp = 0
#' dimnames(result_array)[[2]]
#'
#' # Access a specific simulation
#' sim_1 <- result_array[, , 1]  # First simulation as a matrix
#'
#' # The unadjusted difference in means for each simulation. The unit level
#' # treatment effect is 110, so the gap is the confounding.
#' ate_by_sim <- apply(result_array, 3, function(sim) {
#'   treated <- sim[sim[, "AA"] == 1, "YY"]
#'   control <- sim[sim[, "AA"] == 0, "YY"]
#'   mean(treated) - mean(control)
#' })
#' round(ate_by_sim, 1)
#'
#' @seealso \code{\link{run_parallel_simulation}}, \code{\link{generate_dataset}}
#'
#' @export
convert_to_array <- function(simulation_results) {
  if (!is.list(simulation_results) || length(simulation_results) == 0L) {
    stop("simulation_results must be a non-empty list of data frames")
  }
  if (!all(vapply(simulation_results, is.data.frame, logical(1)))) {
    stop("every element of simulation_results must be a data.frame")
  }

  n_sim     <- length(simulation_results)
  n_obs     <- nrow(simulation_results[[1]])
  var_names <- names(simulation_results[[1]])
  n_vars    <- length(var_names)

  consistent <- vapply(simulation_results, function(df) {
    nrow(df) == n_obs && identical(names(df), var_names)
  }, logical(1))
  if (!all(consistent)) {
    stop("every data frame in simulation_results must have the same dimensions ",
         "and the same variable names; ", sum(!consistent), " of ", n_sim,
         " differ from the first")
  }

  result_array <- array(
    NA_real_,
    dim = c(n_obs, n_vars, n_sim),
    dimnames = list(paste0("obs_", seq_len(n_obs)), var_names,
                    paste0("sim_", seq_len(n_sim))))

  for (i in seq_len(n_sim)) {
    result_array[, , i] <- as.matrix(simulation_results[[i]])
  }
  result_array
}

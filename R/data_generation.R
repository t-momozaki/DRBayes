#' Generate Synthetic Dataset for Causal Inference Simulation
#'
#' Generates a synthetic dataset with confounders, treatment assignment, and outcomes
#' following a specific data generating process designed for testing causal inference methods.
#' The dataset includes potential outcomes under both treatment and control conditions.
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
#'   \item{WW5, ..., WW(4+pp)}{Additional confounding variables (only if pp > 0)}
#' }
#'
#' @details
#' The data generating process follows these steps:
#'
#' **Confounders:**
#' \itemize{
#'   \item 4 primary confounders: \eqn{X_1, X_2, X_3, X_4 \sim N(0, 1)}
#'   \item \code{pp} additional confounders: \eqn{X_{5+j} \sim N(\mu_j, 1)} where \eqn{\mu_j \sim \text{Uniform}(-2, 2)} (only if pp > 0)
#' }
#'
#' **Treatment assignment:**
#' \deqn{P(A = 1 | X) = \text{logit}^{-1}(1 - 0.5 X_1 + 0.25 X_2 + 0.1 X_3)}
#'
#' **Potential outcomes:**
#' \deqn{Y^{(0)} = 100 + 27.4 X_1 + 13.7(X_2 + X_3 + X_4) + \epsilon_0}
#' \deqn{Y^{(1)} = 210 + 27.4 X_1 + 13.7(X_2 + X_3 + X_4) + \epsilon_1}
#'
#' where \eqn{\epsilon_0, \epsilon_1 \sim N(0, 1)} independently.
#'
#' **Observed outcome:**
#' \deqn{Y = A \cdot Y^{(1)} + (1-A) \cdot Y^{(0)}}
#'
#' The true average treatment effect (ATE) is \eqn{E[Y^{(1)} - Y^{(0)}] = 210 - 100 = 110}.
#'
#' @examples
#' \dontrun{
#' # Generate a dataset with only primary confounders
#' data_minimal <- generate_dataset(nn = 100, pp = 0, seed = 123)
#' head(data_minimal)
#' ncol(data_minimal)  # Should be 6 (AA, YY, WW1-WW4)
#'
#' # Generate a dataset with additional confounders
#' data_extended <- generate_dataset(nn = 100, pp = 20, seed = 123)
#' head(data_extended)
#' ncol(data_extended)  # Should be 26 (AA, YY, WW1-WW24)
#'
#' # Check treatment assignment
#' table(data_minimal$AA)
#'
#' # Check outcome distribution by treatment
#' aggregate(YY ~ AA, data = data_minimal, FUN = function(x) c(mean = mean(x), sd = sd(x)))
#'
#' # Estimate propensity scores
#' ps_model <- glm(AA ~ WW1 + WW2 + WW3 + WW4, data = data_minimal, family = binomial)
#' summary(ps_model)
#' }
#'
#' @seealso \code{\link{run_parallel_simulation}}, \code{\link{DRBayes.PC}}
#'
#' @importFrom stats rnorm rbinom runif plogis
#'
#' @export
generate_dataset <- function(nn = 200, pp = 40, seed = NULL) {
  # Input validation
  if (!is.numeric(nn) || nn <= 0 || nn != round(nn)) {
    stop("nn must be a positive integer")
  }

  if (!is.numeric(pp) || pp < 0 || pp != round(pp)) {
    stop("pp must be a non-negative integer")
  }

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1)) {
    stop("seed must be NULL or a single numeric value")
  }

  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }

  ## Data generating process
  ### Primary confounders (4 variables)
  XX <- matrix(stats::rnorm(nn * 4), nn, 4)

  ### Additional confounders with random means (only if pp > 0)
  if (pp > 0) {
    XX.add <- matrix(stats::rnorm(nn * pp, mean = stats::runif(pp, -2, 2)), nn, pp, TRUE)
  } else {
    XX.add <- NULL
  }

  ### Treatment assignment (logistic model)
  AA_ln <- 1 * drop(tcrossprod(c(1, -0.5, 0.25, 0.1), XX))
  # AA    <- stats::rbinom(nn, 1, stats::plogis(AA_ln))
  AA    <- (runif(nn) < stats::plogis(AA_ln))

  ### Potential outcomes
  ee1 <- ee0 <- stats::rnorm(nn, sd = 1)
  YY0 <- 100 + 27.4 * XX[,1] + 13.7 * (XX[,2] + XX[,3] + XX[,4]) + ee0
  YY1 <- 210 + 27.4 * XX[,1] + 13.7 * (XX[,2] + XX[,3] + XX[,4]) + ee1

  ### Observed outcome (switching equation)
  YY <- AA * YY1 + (1 - AA) * YY0

  ### Create dataset
  if (pp > 0) {
    # Include additional confounders
    DD <- data.frame(cbind(AA, YY, XX, XX.add))
    names(DD) <- c("AA", "YY",
                   paste0("WW", 1:4),
                   paste0("WW", 5:(4 + pp)))
  } else {
    # Only primary confounders
    DD <- data.frame(cbind(AA, YY, XX))
    names(DD) <- c("AA", "YY", paste0("WW", 1:4))
  }

  return(DD)
}

#' Run Parallel Simulation Study
#'
#' Executes multiple simulation runs in parallel, generating synthetic datasets
#' using \code{\link{generate_dataset}}. This function is designed for conducting
#' large-scale simulation studies to evaluate causal inference methods.
#'
#' @param num_simulations Integer. Number of simulation replications to run. Default is 1000.
#' @param nn Integer. Number of observations per dataset. Default is 200.
#' @param pp Integer. Number of additional confounders per dataset. Default is 40.
#'   Can be 0 for only the 4 primary confounders.
#' @param n_cores Integer or NULL. Number of CPU cores to use for parallel computation.
#'   If NULL, automatically detects available cores minus 1. Default is NULL.
#' @param strategy Character. Parallel processing strategy. Options are "multisession"
#'   (cross-platform) or "multicore" (Unix-like systems only). Default is "multisession".
#' @param use_lecuyer_rng Logical. Whether to use L'Ecuyer-CMRG random number generator
#'   for better reproducibility in parallel computing. Default is TRUE.
#'
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
#' \dontrun{
#' # Small simulation study with only primary confounders
#' results_minimal <- run_parallel_simulation(
#'   num_simulations = 10,
#'   nn = 100,
#'   pp = 0,  # Only primary confounders
#'   n_cores = 2,
#'   strategy = "multisession"
#' )
#'
#' # Check results
#' length(results_minimal)  # Should be 10
#' dim(results_minimal[[1]])  # Should be 100 x 6 (AA, YY, WW1-WW4)
#'
#' # Small simulation study with additional confounders
#' results_extended <- run_parallel_simulation(
#'   num_simulations = 10,
#'   nn = 100,
#'   pp = 20,
#'   n_cores = 2,
#'   strategy = "multisession"
#' )
#'
#' # Check results
#' length(results_extended)  # Should be 10
#' dim(results_extended[[1]])  # Should be 100 x 26 (AA, YY, WW1-WW24)
#'
#' # Calculate summary statistics across simulations
#' exposure_rates <- sapply(results_minimal, function(df) mean(df$AA))
#' outcome_means <- sapply(results_minimal, function(df) mean(df$YY))
#'
#' mean(exposure_rates)  # Average treatment assignment rate
#' mean(outcome_means)   # Average outcome across simulations
#'
#' # Large-scale simulation (be careful with memory)
#' # results_large <- run_parallel_simulation(
#' #   num_simulations = 1000,
#' #   nn = 500,
#' #   pp = 50,
#' #   n_cores = NULL  # Auto-detect cores
#' # )
#' }
#'
#' @seealso \code{\link{generate_dataset}}, \code{\link{DRBayes.PC}}
#'
#' @importFrom future plan availableCores sequential
#' @importFrom furrr future_map furrr_options
#' @importFrom parallelly availableCores
#'
#' @export
run_parallel_simulation <- function(num_simulations = 1000,
                                    nn = 200,
                                    pp = 40,
                                    n_cores = NULL,
                                    strategy = "multisession",
                                    use_lecuyer_rng = TRUE) {

  # Input validation
  if (!is.numeric(num_simulations) || num_simulations <= 0 || num_simulations != round(num_simulations)) {
    stop("num_simulations must be a positive integer")
  }

  if (!is.numeric(nn) || nn <= 0 || nn != round(nn)) {
    stop("nn must be a positive integer")
  }

  if (!is.numeric(pp) || pp < 0 || pp != round(pp)) {
    stop("pp must be a non-negative integer")
  }

  if (!is.null(n_cores) && (!is.numeric(n_cores) || n_cores <= 0 || n_cores != round(n_cores))) {
    stop("n_cores must be NULL or a positive integer")
  }

  if (!strategy %in% c("multisession", "multicore", "sequential")) {
    stop("strategy must be one of: 'multisession', 'multicore', 'sequential'")
  }

  if (!is.logical(use_lecuyer_rng)) {
    stop("use_lecuyer_rng must be logical (TRUE or FALSE)")
  }

  # Auto-set available cores using parallelly package
  if (is.null(n_cores)) {
    n_cores <- min(parallelly::availableCores() - 1, num_simulations)
  }

  # Set Future plan
  future::plan(strategy, workers = n_cores)

  # Configure RNG handling - use warning level for safety, restore after function
  old_rng_option <- base::getOption("future.rng.onMisuse", "warning")
  old_rng_kind <- NULL

  # Set RNG kind for parallel processing if requested
  if (use_lecuyer_rng) {
    old_rng_kind <- base::RNGkind()
    base::RNGkind("L'Ecuyer-CMRG")
  }

  # Ensure cleanup on function exit
  on.exit({
    options(future.rng.onMisuse = old_rng_option)
    if (!is.null(old_rng_kind)) {
      base::RNGkind(old_rng_kind[1])
    }
    future::plan(future::sequential)
  })
  options(future.rng.onMisuse = "warning")

  # Print simulation information
  cat("Starting parallel data generation...\n")
  cat("Number of simulations:", num_simulations, "\n")
  cat("Sample size:", nn, "\n")
  cat("Number of additional confounders:", pp, "\n")
  if (pp == 0) {
    cat("  (Using only 4 primary confounders)\n")
  }
  cat("Total number of confounders:", 4 + pp, "\n")
  cat("Number of parallel cores:", n_cores, "\n")
  cat("Strategy:", strategy, "\n")
  cat("RNG kind:", if(use_lecuyer_rng) "L'Ecuyer-CMRG" else "Default", "\n\n")

  # Measure execution time
  start_time <- Sys.time()

  # Generate simulation indices for better seed management
  sim_indices <- 1:num_simulations

  # Parallel execution using future_map with proper RNG
  simulation_results <- furrr::future_map(
    .x = sim_indices,
    .f = ~ {
      # Set unique seed for each simulation
      current_seed <- .x * 12345 + 67890  # Ensure seed diversity
      generate_dataset(nn = nn, pp = pp, seed = current_seed)
    },
    .options = furrr::furrr_options(
      seed = TRUE,  # Enable proper RNG handling
      scheduling = 1.0
    ),
    .progress = TRUE
  )

  end_time <- Sys.time()

  # Print completion information
  cat("\n")
  cat("Data generated!\n")
  cat("Execution time:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
  cat("Each dataset has dimensions:", nrow(simulation_results[[1]]), "x", ncol(simulation_results[[1]]), "\n\n")

  return(simulation_results)
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
#' \dontrun{
#' # Generate simulation results with only primary confounders
#' results <- run_parallel_simulation(num_simulations = 5, nn = 100, pp = 0)
#'
#' # Convert to array format
#' result_array <- convert_to_array(results)
#'
#' # Check dimensions: [observations, variables, simulations]
#' dim(result_array)  # Should be [100, 6, 5] for pp=0
#'
#' # Generate simulation results with additional confounders
#' results_extended <- run_parallel_simulation(num_simulations = 5, nn = 100, pp = 10)
#' result_array_extended <- convert_to_array(results_extended)
#' dim(result_array_extended)  # Should be [100, 16, 5] for pp=10
#'
#' # Access specific simulation
#' sim_1 <- result_array[, , 1]  # First simulation as matrix
#'
#' # Calculate mean across simulations for each variable
#' var_means <- apply(result_array, c(1, 2), mean)
#'
#' # Calculate treatment effect for each simulation
#' ate_by_sim <- apply(result_array, 3, function(sim) {
#'   treated <- sim[sim[, "AA"] == 1, "YY"]
#'   control <- sim[sim[, "AA"] == 0, "YY"]
#'   mean(treated) - mean(control)
#' })
#' }
#'
#' @seealso \code{\link{run_parallel_simulation}}, \code{\link{generate_dataset}}
#'
convert_to_array <- function(simulation_results) {
  # Input validation
  if (!is.list(simulation_results)) {
    stop("simulation_results must be a list")
  }

  if (length(simulation_results) == 0) {
    stop("simulation_results cannot be empty")
  }

  # Check that all elements are data.frames with same structure
  if (!all(sapply(simulation_results, is.data.frame))) {
    stop("All elements of simulation_results must be data.frames")
  }

  # Get dimensions
  n_sim <- length(simulation_results)
  n_obs <- nrow(simulation_results[[1]])
  n_vars <- ncol(simulation_results[[1]])
  var_names <- names(simulation_results[[1]])

  # Check consistency across simulations
  consistent_dims <- all(sapply(simulation_results, function(df) {
    nrow(df) == n_obs && ncol(df) == n_vars && identical(names(df), var_names)
  }))

  if (!consistent_dims) {
    stop("All data.frames in simulation_results must have the same dimensions and variable names")
  }

  # Initialize array with appropriate dimension names
  result_array <- base::array(
    dim = c(n_obs, n_vars, n_sim),
    dimnames = list(
      paste0("obs_", 1:n_obs),
      var_names,
      paste0("sim_", 1:n_sim)
    )
  )

  # Fill array with data
  for (i in 1:n_sim) {
    result_array[, , i] <- base::as.matrix(simulation_results[[i]])
  }

  return(result_array)
}

# Example usage and demonstration (not run during package build)
if (FALSE) {
  # Load required libraries
  library(future)
  library(future.apply)
  library(furrr)
  library(parallelly)

  # Test with pp = 0 (only primary confounders)
  cat("=== Testing with pp = 0 ===\n")
  test_data_minimal <- generate_dataset(nn = 50, pp = 0, seed = 123)
  cat("Dataset dimensions with pp=0:", dim(test_data_minimal), "\n")
  cat("Variable names:", names(test_data_minimal), "\n")
  print(head(test_data_minimal))

  # Test with pp > 0 (additional confounders)
  cat("\n=== Testing with pp = 5 ===\n")
  test_data_extended <- generate_dataset(nn = 50, pp = 5, seed = 123)
  cat("Dataset dimensions with pp=5:", dim(test_data_extended), "\n")
  cat("Variable names:", names(test_data_extended), "\n")
  print(head(test_data_extended))

  # Parameter settings for simulation
  num_simulations <- 10   # Number of simulations
  nn <- 100               # Sample size
  pp_values <- c(0, 20)   # Test both cases

  for (pp in pp_values) {
    cat("\n=== Simulation with pp =", pp, "===\n")

    # Execute simulation
    simulation_data <- run_parallel_simulation(
      num_simulations = num_simulations,
      nn = nn,
      pp = pp,
      n_cores = 2,  # Number of cores to use
      strategy = "multisession"
    )

    # Check results
    cat("Number of generated datasets:", length(simulation_data), "\n")
    cat("Structure of each dataset (with", 4 + pp, "confounders):\n")
    str(simulation_data[[1]])

    # Simple statistics of results
    exposure_rates <- base::sapply(simulation_data, function(df) base::mean(df$AA))
    outcome_means  <- base::sapply(simulation_data, function(df) base::mean(df$YY))

    cat("\n== Summary of simulation results ==\n")
    cat("Mean exposure rate:", round(base::mean(exposure_rates), 3), "\n")
    cat("SD of exposure rate:", round(stats::sd(exposure_rates), 3), "\n")
    cat("Mean outcome:", round(base::mean(outcome_means), 2), "\n")
    cat("SD of outcome:", round(stats::sd(outcome_means), 2), "\n")

    # Example of conversion to array format
    result_array <- convert_to_array(simulation_data)
    cat("Array dimensions:", dim(result_array), "\n")  # [observations, variables, simulations]
  }
}

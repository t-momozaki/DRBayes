#' Bayesian Doubly Robust Causal Inference via Bayesian Bootstrap
#'
#' Implements the Bayesian bootstrap-based doubly robust estimator
#' proposed by Saarela et al. (2016). This method combines Bayesian bootstrap
#' resampling with inverse probability of treatment weighting to provide
#' doubly robust causal effect estimation that protects against misspecification
#' of either the outcome model or the propensity score model.
#'
#' @param Y Numeric vector of observed outcomes. Supports both continuous
#'   and binary (0/1) outcomes depending on the \code{family} parameter.
#' @param A Binary vector (0/1) of treatment assignments.
#' @param X.lm Data frame or matrix of covariates for the outcome model.
#'   Should include treatment indicator and interaction terms if desired.
#' @param X.ps Data frame or matrix of covariates for the propensity score model.
#'   Should not include the treatment indicator.
#' @param num_iterations Integer. Number of MCMC iterations for the Bayesian
#'   bootstrap procedure. Default is 3000.
#' @param family Character string specifying the outcome model family.
#'   Options are "gaussian" for continuous outcomes or "binomial" for
#'   binary outcomes. Default is "gaussian".
#'
#' @return A numeric vector of length \code{num_iterations} containing posterior
#'   samples of the average treatment effect (ATE). Each element represents
#'   one draw from the posterior distribution of the causal effect.
#'
#' @details
#' This function implements the importance sampling-based doubly robust
#' estimator from equation (16) in Saarela et al. (2016). The method works
#' by:
#'
#' **Bayesian Bootstrap Procedure:**
#' \enumerate{
#'   \item Draw Dirichlet weights \eqn{\xi_i \sim \text{Dir}(1, ..., 1)}
#'   \item Fit weighted propensity score model: \eqn{\text{logit}(e_i) = \gamma_0 + \gamma^T X_{ps,i}}
#'   \item Fit weighted outcome model: \eqn{E[Y_i] = \phi_0 + \phi^T X_{lm,i}}
#'   \item Compute doubly robust estimate using importance sampling weights
#' }
#'
#' **Doubly Robust Formula:**
#' The estimator combines model-based and inverse probability weighted components:
#' \deqn{\hat{\tau} = \sum_{i=1}^n \xi_i \left[ (Y_i - \hat{\mu}_i) \left( \frac{A_i}{\hat{e}_i} - \frac{1-A_i}{1-\hat{e}_i} \right) + (\hat{\mu}_{1i} - \hat{\mu}_{0i}) \right]}
#'
#' where:
#' \itemize{
#'   \item \eqn{\hat{e}_i} is the estimated propensity score
#'   \item \eqn{\hat{\mu}_i} is the fitted outcome under observed treatment
#'   \item \eqn{\hat{\mu}_{1i}, \hat{\mu}_{0i}} are fitted outcomes under treatment/control
#'   \item \eqn{\xi_i} are Bayesian bootstrap weights
#' }
#'
#' **Double Robustness Property:**
#' The estimator is consistent if either:
#' \itemize{
#'   \item The outcome model is correctly specified, OR
#'   \item The propensity score model is correctly specified
#' }
#'
#' **Model Specifications:**
#' \itemize{
#'   \item Propensity score: Logistic regression with \code{quasibinomial} family
#'   \item Outcome model: Linear regression
#'   \item Both models use Bayesian bootstrap weights for parameter estimation
#' }
#'
#' @note
#' **Important considerations:**
#' \itemize{
#'   \item \code{X.lm} should include treatment indicator and any desired interaction terms
#'   \item \code{X.ps} should NOT include the treatment indicator
#'   \item Large \code{num_iterations} provides more stable posterior estimates
#'   \item Extreme propensity scores may lead to high variance; consider weight truncation
#'   \item The method assumes no unmeasured confounding (ignorability)
#'   \item For binary outcomes (\code{family="binomial"}), Y must contain only 0 and 1 values
#'   \item The ATE interpretation differs: mean difference for continuous, probability difference for binary
#' }
#'
#' @examples
#' \dontrun{
#' # Generate synthetic data
#' set.seed(123)
#' n <- 200
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' X3 <- rnorm(n)
#'
#' # Treatment assignment
#' ps_true <- plogis(0.5 * X1 + 0.3 * X2)
#' A <- rbinom(n, 1, ps_true)
#'
#' # Outcome with treatment effect = 2
#' Y <- 1 + 2 * A + X1 + X2 + X3 + rnorm(n, sd = 1)
#'
#' # Prepare covariate matrices
#' # Outcome model: include treatment and interactions
#' X.lm <- data.frame(
#'   A = A,
#'   X1 = X1, X2 = X2, X3 = X3,
#'   A_X1 = A * X1,  # Treatment-covariate interaction
#'   A_X2 = A * X2
#' )
#'
#' # Propensity score model: exclude treatment
#' X.ps <- data.frame(X1 = X1, X2 = X2, X3 = X3)
#'
#' # Estimate causal effect for continuous outcome
#' posterior_samples <- DRBayes.BB(
#'   Y = Y,
#'   A = A,
#'   X.lm = X.lm,
#'   X.ps = X.ps,
#'   num_iterations = 1000,
#'   family = "gaussian"  # For continuous outcomes
#' )
#'
#' # For binary outcomes
#' Y_binary <- rbinom(n, 1, plogis(0.5 + 1.5 * A + X1 + X2))
#' posterior_samples_binary <- DRBayes.BB(
#'   Y = Y_binary,
#'   A = A,
#'   X.lm = X.lm,
#'   X.ps = X.ps,
#'   num_iterations = 1000,
#'   family = "binomial"  # For binary outcomes
#' )
#'
#' # Posterior summary
#' cat("Posterior mean ATE:", round(mean(posterior_samples), 3), "\n")
#' cat("95% Credible interval:",
#'     round(quantile(posterior_samples, c(0.025, 0.975)), 3), "\n")
#'
#' # Posterior diagnostics
#' hist(posterior_samples, main = "Posterior Distribution of ATE")
#' plot(posterior_samples, type = "l", main = "MCMC Trace")
#'
#' # Compare with simple difference in means
#' naive_ate <- mean(Y[A == 1]) - mean(Y[A == 0])
#' cat("Naive ATE:", round(naive_ate, 3), "\n")
#' }
#'
#' @references
#' Saarela, O., Belzile, L. R., & Stephens, D. A. (2016). A Bayesian view of
#' doubly robust causal inference. \emph{Biometrika}, 103(3), 667-681.
#' \doi{10.1093/biomet/asw025}
#'
#' @seealso
#' \code{\link{DRBayes.PC}}, \code{\link{generate_dataset}},
#' \code{\link{run_parallel_simulation}}
#'
#' @importFrom MCMCpack rdirichlet
#' @importFrom stats glm lm plogis coef quasibinomial
#'
#' @export
DRBayes.BB <- function(Y, A, X.lm, X.ps, num_iterations = 3000, family = "gaussian") {
  # Input validation
  if (!is.numeric(Y)) {
    stop("Y must be a numeric vector")
  }

  if (!all(A %in% c(0, 1))) {
    stop("A must be a binary vector with values 0 and 1")
  }

  if (length(Y) != length(A)) {
    stop("Y and A must have the same length")
  }

  if (nrow(X.lm) != length(Y) || nrow(X.ps) != length(Y)) {
    stop("X.lm and X.ps must have the same number of rows as the length of Y and A")
  }

  if (!is.numeric(num_iterations) || num_iterations <= 0 || num_iterations != round(num_iterations)) {
    stop("num_iterations must be a positive integer")
  }

  if (!family %in% c("gaussian", "binomial")) {
    stop("family must be either 'gaussian' for continuous outcomes or 'binomial' for binary outcomes")
  }

  if (family == "binomial" && !all(Y %in% c(0, 1))) {
    stop("When family='binomial', Y must contain only 0 and 1 values")
  }

  if (any(is.na(Y)) || any(is.na(A)) || any(is.na(X.lm)) || any(is.na(X.ps))) {
    warning("Missing values detected in input data")
  }

  # Convert to data frames for model fitting
  df.lm <- data.frame(Y, X.lm)
  df.ps <- data.frame(A, X.ps)

  # Convert to matrices for linear algebra operations
  X.lm <- as.matrix(X.lm)
  Z.lm <- cbind(1, X.lm)  # Design matrix with intercept

  X.ps <- as.matrix(X.ps)
  Z.ps <- cbind(1, X.ps)  # Design matrix with intercept

  # Design matrices for potential outcomes
  # Assuming first column of X.lm is treatment indicator
  Z.lm1 <- cbind(1, 1, X.lm[, -1])  # Treatment = 1
  Z.lm0 <- cbind(1, 0, X.lm[, -1])  # Treatment = 0

  # Sample size and dimensions
  nn    <- nrow(X.lm)
  pp.lm <- ncol(X.lm)
  pp.ps <- ncol(X.ps)

  # Initialize posterior samples storage
  post.ate <- numeric(num_iterations)

  # Bayesian bootstrap MCMC loop
  for (iter in seq_len(num_iterations)) {

    # Draw Dirichlet weights (Bayesian bootstrap)
    # Each observation gets a random weight from Dir(1,...,1)
    xi <- drop(MCMCpack::rdirichlet(1, rep(1, nn)))

    # Fit weighted propensity score model using Bayesian bootstrap weights
    psmodel <- glm(A ~ .,
                   family  = quasibinomial(link = logit),
                   weights = nn * xi,
                   data    = df.ps)

    # Predict propensity scores
    ps <- plogis(drop(tcrossprod(coef(psmodel), Z.ps)))

    # Fit weighted outcome model using appropriate family
    if (family == "gaussian") {
      # Linear regression for continuous outcomes
      otcmodel <- lm(Y ~ .,
                     weights = nn * xi,
                     data    = df.lm)

      # Predict outcomes under observed and counterfactual treatments
      mu <- drop(tcrossprod(coef(otcmodel), Z.lm))    # Observed treatment
      y1 <- drop(tcrossprod(coef(otcmodel), Z.lm1))   # Treatment = 1
      y0 <- drop(tcrossprod(coef(otcmodel), Z.lm0))   # Treatment = 0

    } else if (family == "binomial") {
      # Logistic regression for binary outcomes
      otcmodel <- glm(Y ~ .,
                      family  = quasibinomial(link = logit),
                      weights = nn * xi,
                      data    = df.lm)

      # Predict probabilities under observed and counterfactual treatments
      mu <- plogis(drop(tcrossprod(coef(otcmodel), Z.lm)))    # Observed treatment
      y1 <- plogis(drop(tcrossprod(coef(otcmodel), Z.lm1)))   # Treatment = 1
      y0 <- plogis(drop(tcrossprod(coef(otcmodel), Z.lm0)))   # Treatment = 0
    }

    # Compute importance sampling weights (clever covariate) for IPW component
    cc <- A / ps - (1 - A) / (1 - ps)

    # Bayesian bootstrap doubly robust estimator:
    # Combines IPW component with regression imputation component
    post.ate[iter] <- sum(xi * (Y - mu) * cc) + sum(xi * (y1 - y0))
  }

  return(post.ate)
}

# Example usage and testing (not run during package build)
if (FALSE) {
  # Load required libraries for testing
  library(MCMCpack)

  # Example with simulated data from run_parallel_simulation
  if (FALSE) {
    simulation_data <- run_parallel_simulation(
      num_simulations = 100,
      nn = 200,
      pp = 40,
      n_cores = 4,
      strategy = "multisession"
    )
  }

  # Use first simulated dataset for demonstration
  data <- simulation_data[[1]]

  # Extract outcomes and treatment
  Y <- data$YY
  A <- data$AA

  # Prepare covariate matrices
  # Outcome model: include treatment, covariates, and interactions
  X.lm <- cbind(A, data[, -(1:2)], A * data[, -(1:2)])

  # Propensity score model: only covariates (no treatment)
  X.ps <- data[, -(1:2)]

  # Estimate causal effect using Bayesian Bootstrap method
  posterior_ate <- DRBayes.BB(Y, A, X.lm, X.ps, num_iterations = 1000, family = "gaussian")

  # Summary statistics
  cat("Posterior mean ATE:", round(mean(posterior_ate), 3), "\n")
  cat("Posterior SD:", round(sd(posterior_ate), 3), "\n")
  cat("95% Credible interval:",
      round(quantile(posterior_ate, c(0.025, 0.975)), 3), "\n")

  # Compare with naive estimator
  naive_ate <- mean(Y[A == 1]) - mean(Y[A == 0])
  cat("Naive ATE (no covariate adjustment):", round(naive_ate, 3), "\n")

  # True ATE for this data generating process is 110
  cat("True ATE: 110\n")

  # Posterior diagnostics
  hist(posterior_ate, main = "Posterior Distribution of ATE",
       xlab = "Average Treatment Effect")
  abline(v = mean(posterior_ate), col = "red", lwd = 2)
  abline(v = 110, col = "blue", lwd = 2, lty = 2)
  legend("topright", c("Posterior Mean", "True ATE"),
         col = c("red", "blue"), lty = c(1, 2), lwd = 2)
}

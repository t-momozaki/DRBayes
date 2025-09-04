#' Bayesian Doubly Robust Causal Inference via Bayesian Bootstrap (nleqslv implementation)
#'
#' Implements the Bayesian bootstrap-based doubly robust estimator
#' proposed by Saarela et al. (2016) using direct score function estimation
#' with nleqslv. This method combines Bayesian bootstrap resampling with
#' inverse probability of treatment weighting to provide doubly robust
#' causal effect estimation that protects against misspecification
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
#' estimator from equation (16) in Saarela et al. (2016) using direct
#' score function optimization via nleqslv. The method works by:
#'
#' **Bayesian Bootstrap Procedure:**
#' \enumerate{
#'   \item Draw Dirichlet weights \eqn{\xi_i \sim \text{Dir}(1, ..., 1)}
#'   \item Solve weighted propensity score equations: \eqn{\sum w_i X_i (A_i - \text{expit}(X_i^T\gamma)) = 0}
#'   \item Solve weighted outcome model equations: \eqn{\sum w_i X_i (Y_i - g(X_i^T\phi)) = 0}
#'   \item Compute doubly robust estimate using importance sampling weights
#' }
#'
#' **Score Function Approach:**
#' Instead of using glm/lm, this implementation directly solves the weighted
#' score equations using nleqslv for numerical root finding. This provides:
#' \itemize{
#'   \item Greater theoretical transparency
#'   \item Direct implementation of the mathematical formulation
#'   \item Flexibility for custom model specifications
#' }
#'
#' @note
#' **Implementation differences from standard glm/lm approach:**
#' \itemize{
#'   \item Uses nleqslv for numerical root finding instead of built-in regression functions
#'   \item May have different convergence properties in edge cases
#'   \item Provides more direct mathematical control over the estimation process
#'   \item Requires careful initialization and convergence monitoring
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
#' X.lm <- data.frame(A = A, X1 = X1, X2 = X2, X3 = X3)
#' X.ps <- data.frame(X1 = X1, X2 = X2, X3 = X3)
#'
#' # Estimate causal effect
#' posterior_samples <- DRBayes.BB.nleqslv(
#'   Y = Y, A = A, X.lm = X.lm, X.ps = X.ps,
#'   num_iterations = 1000, family = "gaussian"
#' )
#'
#' # Posterior summary
#' cat("Posterior mean ATE:", round(mean(posterior_samples), 3), "\n")
#' cat("95% Credible interval:",
#'     round(quantile(posterior_samples, c(0.025, 0.975)), 3), "\n")
#' }
#'
#' @references
#' Saarela, O., Belzile, L. R., & Stephens, D. A. (2016). A Bayesian view of
#' doubly robust causal inference. \emph{Biometrika}, 103(3), 667-681.
#' \doi{10.1093/biomet/asw025}
#'
#' @importFrom MCMCpack rdirichlet
#' @importFrom nleqslv nleqslv
#' @importFrom stats plogis
#'
#' @export
DRBayes.BB.nleqslv <- function(Y, A, X.lm, X.ps, num_iterations = 3000, family = "gaussian") {

  # Input validation (same as original)
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

  # Convert to matrices for linear algebra operations
  X.lm <- as.matrix(X.lm)
  X.ps <- as.matrix(X.ps)

  # Design matrices with intercepts
  Z.lm <- cbind(1, X.lm)  # Design matrix for outcome model
  Z.ps <- cbind(1, X.ps)  # Design matrix for propensity score model

  # Design matrices for potential outcomes
  # Assuming first column of X.lm is treatment indicator
  Z.lm1 <- cbind(1, 1, X.lm[, -1])  # Treatment = 1
  Z.lm0 <- cbind(1, 0, X.lm[, -1])  # Treatment = 0

  # Sample size and dimensions
  nn    <- length(Y)
  pp.lm <- ncol(Z.lm)
  pp.ps <- ncol(Z.ps)

  # Initialize posterior samples storage
  post.ate <- numeric(num_iterations)

  # Define score functions

  # Propensity score model score function
  ps_score <- function(alpha, weights) {
    linear_pred <- drop( base::tcrossprod(alpha, Z.ps) )
    prob        <- plogis(linear_pred)
    residuals   <- A - prob
    score       <- drop( base::crossprod(Z.ps, weights * residuals) )
    return(as.numeric(score))
  }

  # Outcome model score function (depends on family)
  if (family == "gaussian") {
    outcome_score <- function(beta, weights) {
      linear_pred <- drop( base::tcrossprod(beta, Z.lm) )
      residuals   <- Y - linear_pred
      score       <- drop( base::crossprod(Z.lm, weights * residuals) )
      return(as.numeric(score))
    }
  } else if (family == "binomial") {
    outcome_score <- function(beta, weights) {
      linear_pred <- drop( base::tcrossprod(beta, Z.lm) )
      prob        <- plogis(linear_pred)
      residuals   <- Y - prob
      score       <- drop( base::crossprod(Z.lm, weights * residuals) )
      return(as.numeric(score))
    }
  }

  # Bayesian bootstrap MCMC loop
  for (iter in seq_len(num_iterations)) {

    # Draw Dirichlet weights (Bayesian bootstrap)
    xi <- drop(MCMCpack::rdirichlet(1, rep(1, nn)))
    weights <- xi  # Scale weights to preserve effective sample size

    # Solve propensity score model using nleqslv
    ps_result <- nleqslv::nleqslv(
      x = rep(0, pp.ps),
      fn = function(alpha) ps_score(alpha, weights),
      control = list(ftol = 1e-8, xtol = 1e-8)
    )

    if (sqrt(sum(ps_result$fvec^2)) > 1e-6) {
      warning(paste("Propensity score model convergence issue at iteration", iter,
                    "- residual norm:", sqrt(sum(ps_result$fvec^2))))
    }

    alpha_hat <- ps_result$x

    # Solve outcome model using nleqslv
    outcome_result <- nleqslv::nleqslv(
      x = rep(0, pp.lm),
      fn = function(beta) outcome_score(beta, weights),
      control = list(ftol = 1e-8, xtol = 1e-8)
    )

    if (sqrt(sum(outcome_result$fvec^2)) > 1e-6) {
      warning(paste("Outcome model convergence issue at iteration", iter,
                    "- residual norm:", sqrt(sum(outcome_result$fvec^2))))
    }

    beta_hat <- outcome_result$x

    # Predict propensity scores
    ps <- plogis( drop( base::tcrossprod(alpha_hat, Z.ps) ) )

    # Predict outcomes under observed and counterfactual treatments
    if (family == "gaussian") {
      mu <- drop( base::tcrossprod(beta_hat, Z.lm) )           # Observed treatment
      y1 <- drop( base::tcrossprod(beta_hat, Z.lm1) )          # Treatment = 1
      y0 <- drop( base::tcrossprod(beta_hat, Z.lm0) )          # Treatment = 0
    } else if (family == "binomial") {
      mu <- plogis( drop( base::tcrossprod(beta_hat, Z.lm) ) )   # Observed treatment
      y1 <- plogis( drop( base::tcrossprod(beta_hat, Z.lm1) ) )  # Treatment = 1
      y0 <- plogis( drop( base::tcrossprod(beta_hat, Z.lm0) ) )  # Treatment = 0
    }

    # Compute importance sampling weights (clever covariate)
    cc <- A / ps - (1 - A) / (1 - ps)

    # Bayesian bootstrap doubly robust estimator
    post.ate[iter] <- sum(xi * (Y - mu) * cc) + sum(xi * (y1 - y0))
  }

  return(post.ate)
}

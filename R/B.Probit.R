#' Bayesian Probit Regression using Latent Variable Approach
#'
#' @title Bayesian Probit Regression with Latent Variables
#'
#' @description
#' Performs Bayesian probit regression using the latent variable data augmentation
#' scheme proposed by Albert and Chib (1993). This method provides an exact and
#' efficient MCMC algorithm for binary probit regression by introducing latent
#' Gaussian variables that render the likelihood conditionally Gaussian.
#'
#' @param Y A numeric vector of binary outcomes (0 or 1). Must contain only
#'   values 0 and 1.
#' @param X A matrix or data.frame of covariates. Each row corresponds to an
#'   observation and each column to a predictor variable. An intercept term
#'   will be automatically added.
#' @param mc A positive integer specifying the number of MCMC iterations
#'   (default: 5000).
#' @param theta_prior Either a scalar or a square matrix specifying the precision
#'   for the normal prior on regression coefficients. If scalar, assumes
#'   diagonal precision matrix with all elements equal to this value.
#'   If NULL (default), a weakly informative prior with precision
#'   1/100 * I is used, where I is the identity matrix.
#'
#' @return A matrix of posterior samples with dimensions \code{mc} x (\code{ncol(X)} + 1).
#'   Each row represents one MCMC draw of the regression coefficients theta.
#'   Columns are named with "(Intercept)" for the intercept term and
#'   "X1", "X2", etc. for the predictor coefficients.
#'
#' @details
#' The function implements the Albert and Chib (1993) latent variable approach
#' which introduces continuous latent variables W_i such that:
#' \itemize{
#'   \item Y_i = 1 if W_i > 0, Y_i = 0 if W_i ≤ 0
#'   \item W_i ~ N(X_i^T * theta, 1)
#' }
#'
#' The algorithm alternates between:
#' \enumerate{
#'   \item Sampling latent variables W_i from truncated normal distributions
#'   \item Sampling regression coefficients theta from a multivariate normal distribution
#' }
#'
#' For observations with Y_i = 1, W_i is sampled from N(X_i^T * theta, 1) truncated
#' to be positive. For Y_i = 0, W_i is sampled from the same distribution truncated
#' to be negative or zero.
#'
#' The prior for regression coefficients is theta ~ N(0, V), where V^(-1) is
#' specified by \code{theta_prior}. The default prior is relatively uninformative
#' with precision 0.01 for each coefficient.
#'
#' @note
#' The function requires the \code{RcppTN} package for sampling truncated
#' normal random variables.
#'
#' @examples
#' \dontrun{
#' # Simulate binary probit regression data
#' set.seed(123)
#' n <- 200
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' X <- cbind(X1, X2)
#'
#' # True coefficients: intercept = 0.5, X1 = 1.2, X2 = -0.8
#' eta <- 0.5 + 1.2 * X1 - 0.8 * X2
#' prob <- pnorm(eta)  # probit transformation
#' Y <- rbinom(n, 1, prob)
#'
#' # Fit Bayesian probit regression
#' posterior_samples <- B.Probit(Y, X, mc = 1000)
#'
#' # Summarize posterior
#' colMeans(posterior_samples)  # posterior means
#' apply(posterior_samples, 2, sd)  # posterior standard deviations
#'
#' # Compare with classical probit regression
#' classical_fit <- glm(Y ~ X, family = binomial(link = "probit"))
#' cbind(
#'   Bayesian = colMeans(posterior_samples),
#'   Classical = coef(classical_fit)
#' )
#' }
#'
#' @references
#' Albert, J. H., & Chib, S. (1993). Bayesian analysis of binary and polychotomous
#' response data. Journal of the American Statistical Association, 88(422), 669-679.
#'
#' @seealso
#' \code{\link[stats]{glm}} for classical probit regression,
#' \code{\link{B.Logit}} for Bayesian logistic regression,
#' \code{\link[RcppTN]{rtn}} for truncated normal random number generation.
#'
#' @export
#' @importFrom RcppTN rtn
#' @importFrom stats rnorm

B.Probit <- function(Y, X, mc = 5000, theta_prior = NULL) {

  # Input validation
  if (!is.numeric(Y)) {
    stop("Y must be a numeric vector")
  }

  if (!all(Y %in% c(0, 1))) {
    stop("Y must be binary (0 or 1)")
  }

  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data.frame")
  }

  # Convert to matrix and add intercept
  X <- as.matrix(X)
  Z <- cbind(1, X)

  nn <- nrow(X)
  pp <- ncol(X)

  if (length(Y) != nn) {
    stop("Number of observations in Y and X must match")
  }

  if (mc <= 0 || !is.numeric(mc)) {
    stop("mc must be a positive integer")
  }

  # Initialize parameters
  theta <- rep(0, pp+1)
  W     <- rep(0, nn)  # Latent variables

  # Set default prior if not specified
  if (is.null(theta_prior)) {
    theta_prior <- diag(pp+1) * 1/100
  } else if (length(theta_prior) == 1) {
    # If scalar, create diagonal precision matrix
    theta_prior <- diag(pp+1) * theta_prior
  } else if (!is.matrix(theta_prior)) {
    stop("theta_prior must be NULL, a scalar, or a square matrix")
  } else if (nrow(theta_prior) != (pp+1) || ncol(theta_prior) != (pp+1)) {
    stop("theta_prior matrix must have dimensions (p+1) x (p+1) where p is ncol(X)")
  }

  # Pre-compute truncation bounds for efficiency
  # For Y_i = 1: W_i ~ N(mu, 1) truncated to (0, Inf)
  # For Y_i = 0: W_i ~ N(mu, 1) truncated to (-Inf, 0]
  lower_bounds <- ifelse(Y == 1, 0, -Inf)
  upper_bounds <- ifelse(Y == 1, Inf, 0)

  # Storage for posterior samples
  post.theta <- matrix(NA, mc, pp+1)
  colnames(post.theta) <- c("(Intercept)", paste0("X", 1:pp))

  # MCMC sampling
  for (iter in 1:mc) {
    # Step 1: Sample all latent variables W from truncated normal distributions
    mu_W <- drop(tcrossprod(theta, Z))

    # Vectorized sampling from truncated normal distributions
    W <- RcppTN::rtn(
      .mean = mu_W,
      .sd   = rep(1,nn),
      .low  = lower_bounds,
      .high = upper_bounds
    )

    # Step 2: Sample regression coefficients given latent variables
    # This is standard Bayesian linear regression: W ~ N(Z*theta, I)
    Q_theta <- crossprod(Z) + theta_prior  # Z'Z + prior precision
    b_theta <- crossprod(Z, W)             # Z'W (where W is the latent variable vector)

    # Sample from multivariate normal using Cholesky decomposition
    U_theta <- chol(Q_theta)
    theta <- drop(base::backsolve(U_theta,
                                  base::backsolve(U_theta, b_theta, transpose=TRUE) + stats::rnorm(pp+1)))

    # Store samples
    post.theta[iter, ] <- theta
  }

  return(post.theta)
}

#' Bayesian Logistic Regression using Polya-Gamma Data Augmentation
#'
#' @title Bayesian Logistic Regression with Polya-Gamma Latent Variables
#'
#' @description
#' Performs Bayesian logistic regression using the Polya-Gamma data augmentation
#' scheme proposed by Polson, Scott, and Windle (2013). This method provides
#' an exact and efficient MCMC algorithm for binary logistic regression by
#' introducing Polya-Gamma auxiliary variables that render the likelihood
#' conditionally Gaussian.
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
#'   Each row represents one MCMC draw of the regression coefficients.
#'   Columns are named with "(Intercept)" for the intercept term and
#'   "X1", "X2", etc. for the predictor coefficients.
#'
#' @details
#' The function implements the Polya-Gamma data augmentation strategy which
#' represents the logistic likelihood as a scale mixture of Gaussians. The
#' algorithm alternates between:
#' \enumerate{
#'   \item Sampling Polya-Gamma auxiliary variables w_i ~ PG(1, X_i^T * theta)
#'   \item Sampling regression coefficients theta from a multivariate normal distribution
#' }
#'
#' The method is exact (no approximation) and often more efficient than
#' alternative data augmentation schemes, especially for hierarchical models.
#'
#' The prior for regression coefficients is theta ~ N(0, V), where V^(-1) is
#' specified by \code{theta_prior}. The default prior is relatively uninformative
#' with precision 0.01 for each coefficient.
#'
#' @note
#' The function requires the \code{pgdraw} package for sampling Polya-Gamma
#' random variables.
#'
#' @examples
#' \dontrun{
#' # Simulate binary logistic regression data
#' set.seed(123)
#' n <- 200
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' X <- cbind(X1, X2)
#'
#' # True coefficients: intercept = 0.5, X1 = 1.2, X2 = -0.8
#' eta <- 0.5 + 1.2 * X1 - 0.8 * X2
#' prob <- plogis(eta)  # logistic transformation
#' Y <- rbinom(n, 1, prob)
#'
#' # Fit Bayesian logistic regression
#' posterior_samples <- B.Logit(Y, X, mc = 1000)
#'
#' # Summarize posterior
#' colMeans(posterior_samples)  # posterior means
#' apply(posterior_samples, 2, sd)  # posterior standard deviations
#' }
#'
#' @references
#' Polson, N. G., Scott, J. G., & Windle, J. (2013). Bayesian inference for
#' logistic models using Polya-Gamma latent variables. Journal of the American
#' Statistical Association, 108(504), 1339-1349.
#'
#' @seealso
#' \code{\link[stats]{glm}} for classical logistic regression,
#' \code{\link[pgdraw]{pgdraw}} for Polya-Gamma random number generation.
#'
#' @export
#' @importFrom pgdraw pgdraw
#' @importFrom stats rnorm

B.Logit <- function(Y, X, mc = 5000, theta_prior = NULL) {

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
  theta   <- rep(0, pp+1)

  # Transform response for Pólya-Gamma augmentation
  iota <- Y - 0.5

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

  # Storage for posterior samples
  post.theta <- matrix(NA, mc, pp+1)
  colnames(post.theta) <- c("(Intercept)", paste0("X", 1:pp))

  # MCMC sampling
  for (iter in 1:mc) {
    # Sample Pólya-Gamma auxiliary variables
    b_w <- drop( tcrossprod(theta, Z) )
    w   <- pgdraw::pgdraw(1, b_w)

    # Sample regression coefficients
    b_theta <- drop( crossprod(Z, iota) )
    Q_theta <- crossprod(Z, w*Z) + theta_prior
    U_theta <- chol(Q_theta)
    theta   <- drop( base::backsolve(U_theta,
                                     base::backsolve(U_theta, b_theta, transpose=TRUE) + stats::rnorm(pp+1)) )

    # Store samples
    post.theta[iter, ] <- theta
  }

  return(post.theta)
}

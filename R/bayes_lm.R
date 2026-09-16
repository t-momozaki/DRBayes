#' Bayesian Linear Regression using MCMC
#'
#' @description
#' Performs Bayesian linear regression using Gibbs sampling with a normal prior
#' on the regression coefficients and an inverse gamma prior on the error
#' variance. The model is Y ~ N(beta_0 1_n + X beta, sigma^2 I) with
#' theta = (beta_0, beta) distributed as N(0, V) and sigma^2 ~ InvGamma(a, b),
#' where V^(-1) is the precision matrix specified by the user. The two priors
#' are independent of each other, so the full conditional of sigma^2 sees the
#' coefficients only through the residual sum of squares. The default precision
#' is read off the data, which is what keeps it weakly informative whatever the
#' units of Y and of the columns of X; see the section "Prior Specification".
#'
#' @param Y A numeric vector of continuous outcomes (response variable).
#' @param X A matrix or data.frame of covariates. Each row corresponds to an
#'   observation and each column to a predictor variable. An intercept term
#'   will be automatically added.
#' @param mc A positive integer specifying the number of MCMC iterations
#'   (default: 5000).
#' @param chains A positive integer, the number of Markov chains to run
#'   (default: 4). R-hat compares chains, so at least four are recommended.
#' @param init Starting values for the coefficient vector, either NULL
#'   (default) or a list of length \code{chains} whose elements are numeric
#'   vectors of length \code{ncol(X) + 1}, the intercept followed by the
#'   coefficients of the columns of X. When NULL the starting values are
#'   generated from the data and deliberately overdispersed, so that the chains
#'   can reveal sensitivity to where they were started: the ordinary least
#'   squares estimate is computed with \code{\link[stats]{lm.fit}}, chain one
#'   starts there, and each further chain starts at a normal draw centred on it
#'   with a standard deviation of four times that estimate's own standard
#'   error. If the least squares fit is unusable, for instance because the
#'   design is rank deficient, the starting values are drawn from the prior
#'   instead. The scan below draws the coefficients first, from a full
#'   conditional that does not depend on their previous value, so a chain is
#'   started by drawing the error variance from its own full conditional given
#'   these coefficients: a starting value far from the data begins its chain at
#'   a large residual sum of squares, and so at a large error variance. Chains
#'   started in different places therefore stay different for as long as it
#'   takes them to converge, which is what R-hat needs in order to say anything
#'   at all (Vehtari et al., 2021).
#' @param theta_prior Either a scalar or a square matrix specifying the prior
#'   precision of the regression coefficients, that is V^(-1) in
#'   theta ~ N(0, V). A scalar asks for a diagonal precision matrix with every
#'   element equal to it. The precision is on the absolute scale of the data,
#'   the reciprocal of a variance in the units of Y squared, as it is in
#'   \code{\link{bayes_logit}} and \code{\link{bayes_probit}}. If NULL
#'   (default), a precision scaled to the data is used; see "Prior
#'   Specification".
#' @param sigma_prior A numeric vector of length 2 specifying the shape (a) and
#'   scale (b) parameters for the inverse-gamma prior on the error variance
#'   sigma^2. The prior is sigma^2 ~ InvGamma(a, b) with mean b/(a-1) for
#'   a > 1 (default: c(1, 1)).
#'
#' @return An array of posterior samples for the regression coefficients with
#'   dimensions \code{mc} x \code{chains} x (\code{ncol(X)} + 1), indexed by
#'   iteration, then chain, then parameter. This is the layout the posterior
#'   package calls a draws_array. The third dimension is named with
#'   "(Intercept)" for the intercept term and with the column names of X, or
#'   "X1", "X2", etc. when X has none. The first two dimensions are unnamed.
#'
#' @details
#' The function implements a Gibbs sampler for Bayesian linear regression with
#' the following hierarchical model:
#' \deqn{Y_i = \beta_0 + X_i^T \beta + \varepsilon_i, \quad
#'   \varepsilon_i \sim N(0, \sigma^2)}
#' \deqn{\theta = (\beta_0, \beta)^T \sim N(0, V)}
#' \deqn{\sigma^2 \sim InvGamma(a, b)}
#'
#' The algorithm alternates between:
#' \enumerate{
#'   \item Sampling regression coefficients theta from their full conditional
#'         (multivariate normal with precision Z'Z / sigma^2 + V^(-1), where
#'         Z = cbind(1, X))
#'   \item Sampling error variance sigma^2 from its full conditional
#'         (inverse-gamma with shape a + n/2 and scale b + RSS/2, RSS being the
#'         residual sum of squares at the current coefficients)
#' }
#'
#' Both full conditionals are available in closed form, so each step is a
#' single draw and no tuning is needed.
#'
#' @section Prior Specification:
#'
#' \strong{Regression Coefficients:} The prior theta ~ N(0, V) is specified via
#' \code{theta_prior}, which sets V^(-1) (the precision matrix):
#' \itemize{
#'   \item Scalar: creates the diagonal precision matrix V^(-1) = theta_prior I
#'   \item Matrix: uses the provided matrix as the precision V^(-1)
#'   \item NULL: uses the default described next
#' }
#'
#' A single fixed number cannot be weakly informative for data on every scale,
#' because a prior precision is compared with Z'Z / sigma^2 and both sides move
#' with the units of the data. On the misspecified outcome model of Orihara et
#' al. (2025), where the residual variance is about 550 and the treatment
#' coefficient about 110, a fixed precision of 1/100 pulls the posterior mean
#' of that coefficient about 2.5 towards zero. The default is therefore scaled
#' to the data, in the manner of the autoscaled priors of rstanarm and of the
#' scaled default of Gelman et al. (2008):
#' \itemize{
#'   \item each slope is given prior standard deviation 2.5 sd(Y) / s_j, with
#'     s_j the standard deviation of column j of X, so that the default depends
#'     neither on the units of Y nor on the units of any column of X;
#'   \item the intercept is given prior standard deviation
#'     100 (|mean(Y)| + sd(Y)), which is effectively flat.
#' }
#'
#' The intercept needs the separate treatment because scaling alone would not
#' protect it. Shrinking a coefficient towards zero is harmless only where zero
#' is a plausible value, and the intercept carries the location of Y, which in
#' the simulation of Orihara et al. (2025) is about 100 while the residual
#' standard deviation is 1. Centring the columns of X would deal with that too,
#' but it would change what the reported intercept means, so instead the
#' intercept prior is made wide enough to ignore: the precision it contributes,
#' 10^(-4) / (|mean(Y)| + sd(Y))^2, is smaller than the intercept entry of
#' Z'Z / sigma^2 by a factor of at least 10^4 n whenever sigma^2 is no larger
#' than the variance of Y, which holds for any fit no worse than the sample
#' mean.
#'
#' A default read off the data has to survive data that carry no scale, so
#' three fallbacks are built in. A constant Y, or fewer than two observations,
#' leaves sd(Y) zero or undefined, and the unit scale is used in its place; a
#' constant column of X has no spread, and its root mean square is used
#' instead, that being the scale its coefficient is still measured in; a column
#' that is entirely zero falls back on the scale of Y alone. The default is
#' therefore proper, and the coefficient precision matrix stays positive
#' definite, on any design whose columns and outcome share a working range of
#' scales, including a collinear one or one with fewer rows than columns.
#' A precision is the square of a reciprocal scale, so it can leave the range of
#' double precision arithmetic in either direction once Y and a column of X are
#' more than about 1e150 apart. Each one is therefore checked once it has been
#' computed, and a fit whose default cannot be represented stops with a message
#' naming the coefficient responsible, rather than running on a prior that has
#' silently gone improper at one end or degenerate at the other. Name
#' \code{theta_prior} to fix the precision yourself if your data reach that far
#' apart.
#'
#' \strong{Error Variance:} The prior sigma^2 ~ InvGamma(a, b) is specified via
#' \code{sigma_prior = c(a, b)}, and has density proportional to
#' (sigma^2)^(-a-1) exp(-b / sigma^2). Common choices:
#' \itemize{
#'   \item c(1, 1): the default. Proper and heavy tailed, with mode 1/2, median
#'     about 1.44 and no finite mean, so it is weakly informative for an error
#'     variance of order one and is quickly dominated by the data.
#'   \item c(0.01, 0.01): an approximation to the improper Jeffreys prior
#'     1/sigma^2, which is the limit as a and b go to zero. It is not that
#'     prior, and differs from it noticeably when n is small.
#'   \item c(2, 1): proper with mean 1 and infinite variance.
#' }
#'
#' @note
#' This function requires the \code{extraDistr} package for sampling from
#' the inverse-gamma distribution. The current implementation only returns
#' posterior samples for regression coefficients.
#'
#' @examples
#' # Simulate linear regression data
#' set.seed(123)
#' n <- 100
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' X <- cbind(X1, X2)
#'
#' # True coefficients: intercept = 2, X1 = 1.5, X2 = -0.8, sigma^2 = 1
#' Y <- 2 + 1.5 * X1 - 0.8 * X2 + rnorm(n, 0, 1)
#'
#' # Fit Bayesian linear regression with default priors
#' draws <- bayes_lm(Y, X, mc = 1000, chains = 4)
#' dim(draws)  # iterations, chains, parameters
#'
#' # Summarize posterior, pooling the chains
#' apply(draws, 3, mean)  # posterior means
#' apply(draws, 3, sd)  # posterior standard deviations
#' apply(draws, 3, quantile, c(0.025, 0.975))  # 95% credible intervals
#'
#' # Pool the chains into an (iterations * chains) by parameters matrix
#' posterior_samples <- apply(draws, 3, as.vector)
#'
#' @references
#' Gelman, A., Jakulin, A., Pittau, M. G., & Su, Y.-S. (2008). A weakly
#' informative default prior distribution for logistic and other regression
#' models. The Annals of Applied Statistics, 2(4), 1360-1383.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
#' Rubin, D. B. (2013). Bayesian data analysis. CRC press.
#'
#' Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
#' Causal Inference via Posterior Coupling. arXiv preprint arXiv:2506.04868.
#'
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Burkner, P.-C.
#' (2021). Rank-normalization, folding, and localization: an improved R-hat for
#' assessing convergence of MCMC. Bayesian Analysis, 16(2), 667-718.
#'
#' @seealso
#' \code{\link[stats]{lm}} for classical linear regression,
#' \code{\link[extraDistr]{rinvgamma}} for inverse-gamma distribution.
#'
#' @export
#' @importFrom extraDistr rinvgamma
#' @importFrom stats rnorm

bayes_lm <- function(Y, X, mc = 5000, chains = 4L, init = NULL,
                     theta_prior = NULL, sigma_prior = c(1, 1)) {

  X      <- validate_sampler_inputs(Y, X, mc)
  chains <- validate_chains(chains)
  Z      <- cbind(1, X)
  nn     <- nrow(X)
  pp     <- ncol(X)

  prior <- validate_prior_precision(theta_prior, pp + 1L,
                                    default = default_coef_precision(Y, X))
  if (length(sigma_prior) != 2L || !is.numeric(sigma_prior) ||
      any(sigma_prior <= 0)) {
    stop("sigma_prior must be two positive numbers, the shape and scale of ",
         "the inverse gamma prior on the error variance")
  }
  a <- sigma_prior[1]
  b <- sigma_prior[2]

  inits <- resolve_inits(init, chains, Y, Z, family = "gaussian",
                         prior_sd = prior_sd_from_precision(prior))

  post.theta <- draw_storage(mc, chains, X)

  # Z and Y do not change, so their cross products are computed once.
  ZtZ <- crossprod(Z)
  ZtY <- drop(crossprod(Z, Y))
  check_cross_products(ZtZ, ZtY, dimnames(post.theta)[[3]])
  ref <- rss_reference(Y, Z, ZtZ, ZtY, prior)

  # The coefficient prior does not involve sigma^2, so the error variance is
  # informed by the residuals alone and its shape parameter never changes.
  a_sigma <- a + nn / 2

  draw_sigma2 <- function(theta) {
    # Cancellation can leave the residual sum of squares a rounding error
    # below zero when the fit is exact, and a scale has to stay positive.
    b_sigma <- max(b + rss_at(theta, ref) / 2, .Machine$double.eps)
    extraDistr::rinvgamma(1, a_sigma, b_sigma)
  }

  # The chains are independent and could run in parallel, but they run in a
  # plain loop on purpose: the coupling step of drbayes_pc accounts for some
  # 98% of a call, so worker startup and data copying would cost more than the
  # sampling they replace.
  for (chain in seq_len(chains)) {
    # The scan draws theta before reading it, so a chain is started by drawing
    # the error variance from its full conditional given the starting
    # coefficients. That is how init reaches this sampler.
    theta  <- inits[[chain]]
    sigma2 <- draw_sigma2(theta)

    for (iter in seq_len(mc)) {
      theta  <- draw_normal_precision(ZtZ / sigma2 + prior, ZtY / sigma2, iter)
      sigma2 <- draw_sigma2(theta)

      post.theta[iter, chain, ] <- theta
    }
  }

  post.theta
}


#' Reference fit for a numerically stable residual sum of squares
#'
#' The residual sum of squares is needed at every iteration, and recomputing
#' `sum((Y - Z theta)^2)` costs O(n p) each time. The usual way out,
#' `Y'Y - 2 theta'Z'Y + theta'Z'Z theta`, costs O(p^2) and is the same number
#' in exact arithmetic, but it subtracts quantities of order `Y'Y` from one
#' another to leave one of order the residual sum of squares. When Y sits far
#' from the origin there is nothing left of it: for `Y = 1e8 + N(0, 1)` and
#' n = 200 the three terms are of order 1e18 and the residual sum of squares of
#' order 200, so double precision keeps none of its digits and the error
#' variance collapses to zero.
#'
#' Accumulating the same quadratic around a reference fit instead of around the
#' origin keeps the O(p^2) cost and loses nothing, because all three terms are
#' then of the order of the answer. The reference is the posterior mode of the
#' coefficients under a ridge, which exists for any design because the ridge is
#' positive definite.
#'
#' @param Y Response passed to the sampler.
#' @param Z Design matrix including its intercept column.
#' @param ZtZ Cross product `crossprod(Z)`.
#' @param ZtY Cross product `crossprod(Z, Y)`, as a vector.
#' @param ridge Positive definite matrix added to `ZtZ` to place the reference.
#'
#' @return A list holding the reference coefficients, the residual sum of
#'   squares there, `Z'` times the residuals there, and `ZtZ`.
#'
#' @keywords internal
#' @noRd
rss_reference <- function(Y, Z, ZtZ, ZtY, ridge) {
  # The Cholesky factor rather than solve(): the ridged cross product is
  # positive definite by construction, while solve() refuses any matrix whose
  # reciprocal condition number is below the machine epsilon, which a collinear
  # design carrying a prior far weaker than its own cross products reaches. A
  # reference is only a point to expand around, so neither route is worth
  # failing over.
  theta <- tryCatch({
    U <- chol(ZtZ + ridge)
    drop(backsolve(U, backsolve(U, ZtY, transpose = TRUE)))
  }, error = function(e) NULL)

  if (is.null(theta) || !all(is.finite(theta))) {
    # Least squares instead, through the pivoted QR of lm.fit(), which answers
    # a rank deficient design by aliasing the columns it drops rather than by
    # failing. It is the better reference of the two where it exists, since Z'r
    # then vanishes and the expansion in rss_at() has no linear term at all.
    theta <- unname(stats::lm.fit(Z, Y)$coefficients)
    theta[is.na(theta)] <- 0
  }

  resid <- Y - drop(Z %*% theta)
  list(theta = theta, rss = sum(resid^2),
       Zt_resid = drop(crossprod(Z, resid)), ZtZ = ZtZ)
}


#' Refuse a design whose cross products have overflowed
#'
#' `crossprod()` returns Inf for a column whose squares are larger than a
#' double can hold. Every step after it then sees NaN, and the first to
#' complain is the Cholesky factorisation of the coefficient precision, whose
#' message blames collinearity or an underflowed shrinkage parameter. Neither
#' is the cause: the design is unusable because a column is on a scale double
#' precision arithmetic cannot square, and dividing that column by its own
#' spread fixes it. Say so here, where the column can still be named.
#'
#' @param ZtZ Cross product `crossprod(Z)` of the design, intercept included.
#' @param ZtY Cross product `crossprod(Z, Y)`, as a vector.
#' @param column_names Names of the columns of `Z`, the intercept first.
#'
#' @return `NULL`, invisibly. Called for the error it raises.
#'
#' @keywords internal
#' @noRd
check_cross_products <- function(ZtZ, ZtY, column_names) {
  bad <- !is.finite(ZtY) | apply(!is.finite(ZtZ), 1L, any)
  if (!any(bad)) {
    return(invisible(NULL))
  }
  stop("The cross products of the design are not finite, so there is no ",
       "posterior to sample from. A column whose squares leave the range of ",
       "double precision arithmetic does this, and here the affected rows of ",
       "X'X are ", toString(column_names[bad]), ". Divide those columns of X, ",
       "or Y itself, by their standard deviations and fit again: under the ",
       "default prior that changes nothing but the units the coefficients are ",
       "reported in.", call. = FALSE)
}


#' Residual sum of squares at a coefficient vector
#'
#' Expands `||Y - Z theta||^2` around the reference fit of [rss_reference()]:
#' with `d = theta - theta_ref` it is
#' `rss_ref - 2 d'Z'r_ref + d'Z'Z d`.
#'
#' @param theta Coefficient vector.
#' @param ref Value returned by [rss_reference()].
#'
#' @return The residual sum of squares, a single number.
#'
#' @keywords internal
#' @noRd
rss_at <- function(theta, ref) {
  delta <- theta - ref$theta
  ref$rss - 2 * sum(delta * ref$Zt_resid) +
    drop(crossprod(delta, ref$ZtZ %*% delta))
}


#' Weakly informative coefficient precision on the scale of the data
#'
#' The default prior of [bayes_lm()] and [bayes_lm_hs()] has to be weakly
#' informative whatever the units of Y and of the columns of X, so it is read
#' off the data, in the manner of the autoscaled priors of rstanarm. Each slope
#' is given prior standard deviation `spread * sd(Y) / s_j`, with `s_j` the
#' spread of column j of X; the intercept is given `flat * (|mean(Y)| + sd(Y))`,
#' wide enough to leave the location of Y entirely to the likelihood, which a
#' scaled prior on its own would not do.
#'
#' The fallbacks are what makes the default safe on degenerate data. A constant
#' Y, or a single observation, has no spread and is given the unit scale; a
#' constant column of X has no spread either, but its coefficient is still
#' measured in units of one over that column, which is what its root mean
#' square supplies; a column of zeros has neither and falls back on the scale
#' of Y. Every value returned is finite and strictly positive, and so is its
#' reciprocal, so the prior is proper, the coefficient precision matrix is
#' positive definite however rank deficient the design, and the prior standard
#' deviation that places the starting values can be read back off it.
#'
#' @param Y Response passed to the sampler.
#' @param X Covariate matrix, without its intercept column.
#' @param spread Prior standard deviation of a slope, in standard deviations of
#'   Y per standard deviation of the column.
#' @param flat Multiple of the magnitude of Y used for the intercept.
#'
#' @return A numeric vector of length `ncol(X) + 1`: the prior precision of the
#'   intercept, followed by one per column of X.
#'
#' @keywords internal
#' @noRd
default_coef_precision <- function(Y, X, spread = 2.5, flat = 100) {
  scale_y <- stats::sd(Y)
  if (!isTRUE(is.finite(scale_y)) || scale_y <= 0) {
    scale_y <- 1
  }

  column_names <- colnames(X)

  precision <- numeric(ncol(X) + 1L)
  precision[1L] <- usable_precision(flat * (abs(mean(Y)) + scale_y),
                                    "the intercept")
  for (j in seq_len(ncol(X))) {
    column <- X[, j]
    scale_x <- stats::sd(column)
    if (!isTRUE(is.finite(scale_x)) || scale_x <= 0) {
      scale_x <- sqrt(mean(column^2))
    }
    if (!isTRUE(is.finite(scale_x)) || scale_x <= 0) {
      scale_x <- 1
    }
    # The label is only ever used in an error, and R does not evaluate an
    # argument it does not reach, so building it costs nothing on a fit that
    # works.
    precision[j + 1L] <- usable_precision(
      spread * scale_y / scale_x,
      if (is.null(column_names) || !nzchar(column_names[j])) {
        paste0("column ", j, " of X")
      } else {
        paste0("column ", column_names[j], " of X")
      })
  }
  precision
}


#' Turn a prior standard deviation read off the data into a usable precision
#'
#' A precision is one over the square of a scale, and a scale read off the data
#' can be anything, so the square leaves the range of double precision
#' arithmetic in one direction or the other once the outcome and the columns of
#' the design are more than about 1e150 apart. Both ends look like a working
#' prior and behave like a broken one. A precision that has underflowed to zero
#' is an improper prior whose precision matrix has no Cholesky factor, and the
#' factorisation fails without saying which coefficient it was; a precision of
#' infinity raises nothing at all, and simply returns that coefficient as
#' exactly zero in every draw. Neither is worth guessing a value for, so the
#' precision is checked here, while the coefficient it belongs to can still be
#' named.
#'
#' @param prior_sd Prior standard deviation, on the scale of the data.
#' @param coefficient Which coefficient it belongs to, for the error message.
#'
#' @return One over `prior_sd` squared: a finite, strictly positive number
#'   whose reciprocal is also finite.
#'
#' @keywords internal
#' @noRd
usable_precision <- function(prior_sd, coefficient) {
  # (1 / sd)^2 rather than 1 / sd^2, because squaring the standard deviation
  # first overflows on data the reciprocal still handles.
  precision <- (1 / prior_sd)^2
  # The variance has to be representable as well as the precision: the prior
  # standard deviation is read back off the precision matrix to place the
  # starting values of the chains.
  if (isTRUE(is.finite(precision)) && precision > 0 &&
      isTRUE(is.finite(1 / precision))) {
    return(precision)
  }
  stop("The default prior for ", coefficient, " works out to a standard ",
       "deviation of ", format(prior_sd, digits = 3), ", whose precision is ",
       "too small, or too large, to hold alongside the variance it stands ",
       "for. The default is read off the scale of the data, so this means Y ",
       "sits enormously far from the origin, or the spread of a column of X ",
       "is enormously far from the spread of Y. Divide them by their standard ",
       "deviations and fit again; bayes_lm() will also take the precision ",
       "straight from you, through theta_prior.", call. = FALSE)
}

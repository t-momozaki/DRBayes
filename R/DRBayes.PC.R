#' Bayesian Doubly Robust Causal Inference via Posterior Coupling
#'
#' This function implements Bayesian doubly robust estimation for average treatment effects (ATE)
#' using posterior coupling. The method incorporates propensity score information via moment
#' conditions directly into the posterior distribution, avoiding the feedback problem and
#' enabling a fully Bayesian interpretation of DR estimation without requiring two-step estimation.
#'
#' The function supports both continuous outcomes (linear regression) and binary outcomes
#' (logistic and probit regression) through the family argument, using an intuitive formula
#' interface that integrates naturally with the R ecosystem.
#'
#' The function supports two modes of operation: internal posterior sampling using built-in
#' model functions, or external posterior samples from other Bayesian software packages
#' (e.g., Stan, JAGS, brms).
#'
#' @param outcome.formula Formula specifying the outcome model.
#'   For example: Y ~ A + X1 + X2 + A:X1. The response variable should be on the left side
#'   and predictors including treatment variable on the right side. Interactions and
#'   transformations are fully supported.
#' @param ps.formula Formula specifying the propensity score model.
#'   For example: A ~ X1 + X2 + X3. The treatment variable should be on the left side
#'   and confounders on the right side.
#' @param data Data.frame containing all variables specified in the formulas.
#' @param family Character string specifying the outcome model family. If NULL (default),
#'   automatically detects based on Y values (binary → "binomial", otherwise → "gaussian").
#'   Options are:
#'   \itemize{
#'     \item NULL (default): Automatically detect based on Y values
#'     \item "gaussian": Linear regression for continuous outcomes
#'     \item "binomial": Logistic or probit regression for binary outcomes
#'   }
#' @param link Character string specifying the link function. If NULL (default),
#'   automatically determined from family: "identity" for gaussian, "logit" for binomial.
#'   Options are: "identity", "logit", "probit".
#' @param outcome.model Function for fitting the outcome model. Should take arguments Y, X, and mc.
#'   Available functions include \code{\link{B.LM}} for continuous outcomes, \code{\link{B.Logit}}
#'   and \code{\link{B.Probit}} for binary outcomes, and \code{\link{HS.LM}}, \code{\link{HS.Logit}},
#'   and \code{\link{HS.Probit}} for horseshoe priors. Set to NULL when using external samples. Default is NULL.
#' @param ps.model Function for fitting the propensity score model. Should take arguments Y, X, and mc.
#'   Available functions include \code{\link{B.Logit}}, \code{\link{B.Probit}}, \code{\link{HS.Logit}},
#'   and \code{\link{HS.Probit}} for binary treatments. Set to NULL when using external samples. Default is NULL.
#' @param outcome.samples Matrix of posterior samples for outcome model coefficients. Each row
#'   represents one MCMC iteration and each column represents a coefficient. The first column
#'   should be the intercept, followed by coefficients for covariates in the outcome model.
#'   Set to NULL when using internal model functions. Default is NULL.
#' @param ps.samples Matrix of posterior samples for propensity score model coefficients. Each row
#'   represents one MCMC iteration and each column represents a coefficient. The first column
#'   should be the intercept, followed by coefficients for covariates in the propensity score model.
#'   Set to NULL when using internal model functions. Default is NULL.
#' @param mc Integer. Number of MCMC iterations for internal sampling. Ignored when using
#'   external samples. Default is 5000.
#' @param bn Integer. Number of burn-in iterations to discard. Applied to both internal
#'   and external samples. Default is 1000.
#' @param thin Integer. Thinning interval for MCMC samples. Applied to both internal
#'   and external samples. Default is 2.
#' @param outcome.priors List of prior distribution parameters to pass to the
#'   outcome model function. For example, when using \code{B.LM}, you can specify
#'   \code{list(theta_prior = 1/50, sigma_prior = c(2, 1))} to customize
#'   the prior distributions for regression coefficients and error variance.
#'   Ignored when using external samples. Default is an empty list (uses function defaults).
#' @param ps.priors List of prior distribution parameters to pass to the
#'   propensity score model function. Ignored when using external samples.
#'   Default is an empty list (uses function defaults).
#' @param na.action Character string specifying how to handle missing values.
#'   Options are "na.omit" (default), "na.fail", or "na.exclude".
#'
#' @return A list containing:
#' \describe{
#'   \item{g.comp}{Numeric vector of posterior samples for ATE using g-computation.}
#'   \item{pc}{Numeric vector of posterior samples for ATE using posterior coupling.}
#'   \item{family}{Character string indicating the outcome model family used.}
#'   \item{link}{Character string indicating the link function used.}
#'   \item{call}{The matched call.}
#'   \item{data_info}{Summary information about the processed data.}
#' }
#'
#' @details
#' **Formula Interface:**
#'
#' DRBayes.PC uses R's standard formula notation, making it intuitive and consistent with
#' other statistical functions like lm() and glm(). The formula interface automatically
#' handles missing values, creates appropriate design matrices, and manages factor variables.
#' \preformatted{
#' # Basic usage
#' result <- DRBayes.PC(
#'   outcome.formula = Y ~ A + X1 + X2,
#'   ps.formula = A ~ X1 + X2 + X3,
#'   data = my_data,
#'   outcome.model = B.LM,
#'   ps.model = B.Logit
#' )
#'
#' # With interactions
#' result <- DRBayes.PC(
#'   outcome.formula = Y ~ A + X1 + X2 + A:X1,
#'   ps.formula = A ~ X1 + X2 + X3 + I(X1^2),
#'   data = my_data,
#'   outcome.model = B.LM,
#'   ps.model = B.Logit
#' )
#' }
#'
#' **Model Families and Link Functions:**
#'
#' The function supports different outcome model types:
#' \itemize{
#'   \item \strong{Automatic detection} (family=NULL): Determines family based on Y values.
#'   If all Y values are 0 or 1, uses "binomial"; otherwise uses "gaussian".
#'   \item \strong{Gaussian family} (family="gaussian"): For continuous outcomes using linear regression
#'   with identity link. ATE is estimated as the mean difference E[Y|A=1,X] - E[Y|A=0,X].
#'   \item \strong{Binomial family} (family="binomial"): For binary outcomes using logistic regression
#'   (link="logit") or probit regression (link="probit"). ATE is estimated as the probability
#'   difference P(Y=1|A=1,X) - P(Y=1|A=0,X).
#' }
#'
#' **Sampling Methods:**
#'
#' The function automatically detects the sampling mode based on the provided arguments:
#' \itemize{
#'   \item \strong{Internal sampling}: Provide \code{outcome.model} and \code{ps.model}
#'   \item \strong{External samples}: Provide \code{outcome.samples} and \code{ps.samples}
#' }
#'
#' If both internal and external options are provided, external samples take precedence
#' with a warning message.
#'
#' **Internal Sampling Mode:**
#' The function generates posterior samples from separate outcome and propensity score models
#' using the provided model functions and performs the following steps:
#' \enumerate{
#'   \item Generates posterior samples from separate outcome and propensity score models
#'   \item Applies sequential Monte Carlo (SMC) to enforce the moment condition constraint
#'   \item Computes posterior distributions for the average treatment effect
#' }
#'
#' **External Samples Mode:**
#' When using external posterior samples (e.g., from Stan, JAGS, brms), ensure that:
#' \itemize{
#'   \item Samples are in matrix format with rows = iterations, columns = coefficients
#'   \item The first column contains intercept terms
#'   \item Coefficient ordering matches the column order in the design matrices created from formulas
#'   \item Both outcome.samples and ps.samples have the same number of iterations
#' }
#'
#' The moment condition enforced is:
#' \deqn{E[(A - \pi(X)) \cdot (Y - \mu(X)) / (\pi(X)(1-\pi(X)))] = 0}
#'
#' where \eqn{\pi(X)} is the propensity score and \eqn{\mu(X)} is the outcome model.
#'
#' @section Link Functions:
#' \itemize{
#'   \item \strong{Identity} (family="gaussian"): \eqn{\mu = X\beta}
#'   \item \strong{Logit} (family="binomial", link="logit"): \eqn{\mu = \text{logit}^{-1}(X\beta) = \frac{e^{X\beta}}{1 + e^{X\beta}}}
#'   \item \strong{Probit} (family="binomial", link="probit"): \eqn{\mu = \Phi(X\beta)} where \eqn{\Phi} is the standard normal CDF
#' }
#'
#' @section External Software Integration:
#' This function is designed to work with posterior samples from various Bayesian software:
#'
#' \strong{Stan/cmdstanr:}
#' \preformatted{
#' # Extract samples
#' outcome.samples <- fit$draws("beta", format = "matrix")
#' ps.samples <- fit$draws("gamma", format = "matrix")
#'
#' # Use in DRBayes.PC
#' result <- DRBayes.PC(
#'   outcome.formula = Y ~ A + X1 + X2 + A:X1,
#'   ps.formula = A ~ X1 + X2 + X3,
#'   data = data,
#'   outcome.samples = outcome.samples,
#'   ps.samples = ps.samples
#' )
#' }
#'
#' \strong{brms:}
#' \preformatted{
#' # Extract fixed effects
#' outcome.samples <- as.matrix(brms_fit)[, 1:ncol(model.matrix(outcome.formula, data))]
#' ps.samples <- as.matrix(brms_fit)[, 1:ncol(model.matrix(ps.formula, data))]
#'
#' # Use in DRBayes.PC
#' result <- DRBayes.PC(
#'   outcome.formula = Y ~ A + X1 + X2,
#'   ps.formula = A ~ X1 + X2 + X3,
#'   data = data,
#'   outcome.samples = outcome.samples,
#'   ps.samples = ps.samples
#' )
#' }
#'
#' \strong{JAGS/R2jags:}
#' \preformatted{
#' # Extract parameter samples
#' outcome.samples <- jags.fit$BUGSoutput$sims.matrix[, grep("beta", colnames(...))]
#' ps.samples <- jags.fit$BUGSoutput$sims.matrix[, grep("gamma", colnames(...))]
#'
#' # Use in DRBayes.PC
#' result <- DRBayes.PC(
#'   outcome.formula = Y ~ A + X1 + X2,
#'   ps.formula = A ~ X1 + X2 + X3,
#'   data = data,
#'   outcome.samples = outcome.samples,
#'   ps.samples = ps.samples
#' )
#' }
#'
#' @references
#' Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust Causal Inference via
#' Posterior Coupling. arXiv preprint arXiv:2506.04868.
#'
#' @examples
#' \dontrun{
#' # ========================================================================
#' # Example 1: Continuous outcome with internal sampling
#' # ========================================================================
#' set.seed(123)
#' n <- 200
#' data <- data.frame(
#'   X1 = rnorm(n),
#'   X2 = rnorm(n),
#'   X3 = rnorm(n)
#' )
#' # Generate treatment with confounding
#' data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$X1 - 0.3*data$X2))
#' # Generate continuous outcome with treatment effect and interactions
#' data$Y <- 1 + 1.5*data$A + 0.8*data$X1 + 0.6*data$X2 + 0.4*data$A*data$X1 + rnorm(n)
#'
#' # Fit using internal sampling
#' result_continuous <- DRBayes.PC(
#'   outcome.formula = Y ~ A + X1 + X2 + A:X1,
#'   ps.formula = A ~ X1 + X2 + X3,
#'   data = data,
#'   family = "gaussian",
#'   outcome.model = B.LM,
#'   ps.model = B.Logit,
#'   mc = 3000, bn = 1000, thin = 2
#' )
#'
#' # View results
#' cat("ATE estimate:", round(mean(result_continuous$pc), 3), "\n")
#' cat("Posterior SD:", round(sd(result_continuous$pc), 3), "\n")
#' cat("95% CI:", round(quantile(result_continuous$pc, c(0.025, 0.975)), 3), "\n")
#'
#' # ========================================================================
#' # Example 2: Binary outcome with logistic regression
#' # ========================================================================
#' # Generate binary outcome
#' data$Y_binary <- rbinom(n, 1, plogis(0.5 + 0.8*data$A + 0.4*data$X1 + 0.3*data$X2))
#'
#' result_binary <- DRBayes.PC(
#'   outcome.formula = Y_binary ~ A + X1 + X2,
#'   ps.formula = A ~ X1 + X2 + X3,
#'   data = data,
#'   family = "binomial",
#'   link = "logit",
#'   outcome.model = B.Logit,
#'   ps.model = B.Logit,
#'   mc = 3000, bn = 1000, thin = 2
#' )
#'
#' # Results represent probability differences for binary outcomes
#' cat("ATE (probability difference):", round(mean(result_binary$pc), 3), "\n")
#'
#' # ========================================================================
#' # Example 3: External sampling integration
#' # ========================================================================
#' # Simulate external posterior samples (e.g., from Stan, JAGS, brms)
#' n_samples <- 1000
#' n_outcome_params <- ncol(model.matrix(Y ~ A + X1 + X2 + A:X1, data))
#' n_ps_params <- ncol(model.matrix(A ~ X1 + X2 + X3, data))
#'
#' # Simulated external samples (replace with actual Stan/JAGS/brms output)
#' outcome.samples <- matrix(rnorm(n_samples * n_outcome_params),
#'                          nrow = n_samples, ncol = n_outcome_params)
#' ps.samples <- matrix(rnorm(n_samples * n_ps_params),
#'                     nrow = n_samples, ncol = n_ps_params)
#'
#' # Use external samples with DRBayes
#' result_external <- DRBayes.PC(
#'   outcome.formula = Y ~ A + X1 + X2 + A:X1,
#'   ps.formula = A ~ X1 + X2 + X3,
#'   data = data,
#'   outcome.samples = outcome.samples,
#'   ps.samples = ps.samples,
#'   bn = 100, thin = 2
#' )
#'
#' # ========================================================================
#' # Example 4: Complex formulas with transformations
#' # ========================================================================
#' result_complex <- DRBayes.PC(
#'   outcome.formula = Y ~ A + log(abs(X1) + 1) + I(X2^2) + poly(X3, 2) + A:X1,
#'   ps.formula = A ~ X1 + X2 + X3 + I(X1^2) + X2:X3,
#'   data = data,
#'   outcome.model = B.LM,
#'   ps.model = B.Logit,
#'   mc = 2000, bn = 500
#' )
#'
#' # ========================================================================
#' # Example 5: Missing data handling
#' # ========================================================================
#' # Create data with missing values
#' data_missing <- data
#' data_missing$X1[sample(n, 20)] <- NA
#' data_missing$X2[sample(n, 15)] <- NA
#'
#' result_missing <- DRBayes.PC(
#'   outcome.formula = Y ~ A + X1 + X2 + X3,
#'   ps.formula = A ~ X1 + X2 + X3,
#'   data = data_missing,
#'   na.action = "na.omit",
#'   outcome.model = B.LM,
#'   ps.model = B.Logit,
#'   mc = 2000, bn = 500
#' )
#'
#' cat("Original sample size:", n, "\n")
#' cat("Effective sample size:", result_missing$data_info$n_observations, "\n")
#'
#' # ========================================================================
#' # Compare results across examples
#' # ========================================================================
#' compare_results <- function(result, name) {
#'   cat(name, ":\n")
#'   cat("  ATE estimate:", round(mean(result$pc), 4), "\n")
#'   cat("  Posterior SD:", round(sd(result$pc), 4), "\n")
#'   cat("  95% CI: [", round(quantile(result$pc, c(0.025, 0.975)), 4), "]\n\n")
#' }
#'
#' compare_results(result_continuous, "Continuous outcome")
#' compare_results(result_binary, "Binary outcome")
#' compare_results(result_external, "External sampling")
#' compare_results(result_complex, "Complex formulas")
#' }
#'
#' @seealso
#' \code{\link{B.LM}} for Bayesian linear regression,
#' \code{\link{B.Logit}} for Bayesian logistic regression,
#' \code{\link{B.Probit}} for Bayesian probit regression,
#' \code{\link{HS.LM}} for horseshoe linear regression,
#' \code{\link{HS.Logit}} for horseshoe logistic regression,
#' \code{\link{HS.Probit}} for horseshoe probit regression,
#' \code{\link{DRBayes.BB}} for Bayesian bootstrap approach.
#'
#' @importFrom stats cov plogis pnorm model.matrix model.frame terms na.omit na.fail na.exclude
#' @importFrom mvtnorm rmvnorm
#' @importFrom MCMCpack rdirichlet
#'
#' @export
DRBayes.PC <- function(outcome.formula, ps.formula, data,
                       family = NULL, link = NULL,
                       outcome.model   = NULL, ps.model   = NULL,
                       outcome.samples = NULL, ps.samples = NULL,
                       mc = 5000, bn = 1000, thin = 2,
                       outcome.priors = list(), ps.priors = list(),
                       na.action = "na.omit") {

  # Store the call for reference
  call_info <- match.call()

  # ============================================================================
  # 1. Input validation and formula processing
  # ============================================================================

  # Validate formula arguments
    if (!inherits(outcome.formula, "formula")) {
      stop("outcome.formula must be a formula object (e.g., Y ~ A + X1 + X2)")
    }

    if (!inherits(ps.formula, "formula")) {
      stop("ps.formula must be a formula object (e.g., A ~ X1 + X2)")
    }

    if (!is.data.frame(data)) {
      stop("data must be a data.frame")
    }

    # Validate na.action
    na.action     <- match.arg(na.action, choices = c("na.omit", "na.fail", "na.exclude"))
    na.action.fun <- switch(na.action,
                           "na.omit" = na.omit,
                           "na.fail" = na.fail,
                           "na.exclude" = na.exclude)

    # Extract variable names from formulas
    outcome_terms <- terms(outcome.formula)
    ps_terms      <- terms(ps.formula)

    # Get response variable names
    outcome_response <- all.vars(outcome.formula)[1]
    ps_response      <- all.vars(ps.formula)[1]

    # Check if variables exist in data
    outcome_vars <- all.vars(outcome.formula)
    ps_vars      <- all.vars(ps.formula)
    all_vars     <- unique(c(outcome_vars, ps_vars))

    missing_vars <- setdiff(all_vars, names(data))
    if (length(missing_vars) > 0) {
      stop("The following variables are not found in data: ", paste(missing_vars, collapse = ", "))
    }

    # Create model frames with missing value handling
    outcome_frame <- model.frame(outcome.formula, data = data, na.action = na.action.fun)
    ps_frame      <- model.frame(ps.formula, data = data, na.action = na.action.fun)

    # Check if both frames have the same number of observations after NA handling
    if (nrow(outcome_frame) != nrow(ps_frame)) {
      warning("Different number of observations in outcome and propensity score models after handling missing values. ",
              "This may indicate different missing value patterns. Consider using complete cases.")
    }

    # Find common rows (intersection of non-missing cases)
    outcome_rownames <- as.numeric(rownames(outcome_frame))
    ps_rownames      <- as.numeric(rownames(ps_frame))
    common_rows      <- intersect(outcome_rownames, ps_rownames)

    if (length(common_rows) == 0) {
      stop("No common observations found between outcome and propensity score models after handling missing values")
    }

    # Subset to common rows
    common_data <- data[common_rows, ]

    # Re-create model frames with common data
    outcome_frame <- model.frame(outcome.formula, data = common_data, na.action = na.fail)
    ps_frame      <- model.frame(ps.formula,      data = common_data, na.action = na.fail)

    # Extract response variables
    Y <- outcome_frame[, 1]
    A <- ps_frame[, 1]

    # Validate treatment variable
    if (!all(A %in% c(0, 1))) {
      stop("Treatment variable (", ps_response, ") must be binary (0 or 1). Found values: ",
           toString(unique(A)))
    }

    # Create design matrices (without intercept since it's added later)
    X.lm <- model.matrix(outcome_terms, data = outcome_frame)[, -1, drop = FALSE]
    X.ps <- model.matrix(ps_terms,      data = ps_frame)[, -1, drop = FALSE]

    # Store data information
    data_info <- list(
      n_observations       = length(Y),
      n_treated            = sum(A),
      n_control            = sum(1 - A),
      outcome_variables    = ncol(X.lm),
      ps_variables         = ncol(X.ps),
      outcome_formula      = outcome.formula,
      ps_formula           = ps.formula,
      missing_observations = nrow(data) - length(Y)
    )

    # Store treatment variable name for potential outcomes matrix creation
    treatment_var_name <- all.vars(ps.formula)[1]

    cat("Formula processing complete:\n")
    cat("  Observations used:", length(Y), "out of", nrow(data), "\n")
    cat("  Treatment group:", sum(A), "| Control group:", sum(1-A), "\n")
    cat("  Outcome model variables:", ncol(X.lm), "\n")
    cat("  Propensity score variables:", ncol(X.ps), "\n")

  # ============================================================================
  # 2. Family and link function determination
  # ============================================================================

  # Determine family
  if (is.null(family)) {
    if (all(Y %in% c(0, 1))) {
      family <- "binomial"
      cat("Auto-detected binary outcome: using family='binomial'\n")
    } else {
      family <- "gaussian"
      cat("Auto-detected continuous outcome: using family='gaussian'\n")
    }
  } else {
    # Validate family argument
    family <- match.arg(family, choices = c("gaussian", "binomial"))
  }

  # Determine link function
  if (is.null(link)) {
    link <- switch(family,
                   "gaussian" = "identity",
                   "binomial" = "logit")
  } else {
    # Validate link argument
    valid_links <- switch(family,
                          "gaussian" = "identity",
                          "binomial" = c("logit", "probit"))
    if (!link %in% valid_links) {
      stop(sprintf("Link '%s' not supported for family '%s'. Valid links: %s",
                   link, family, paste(valid_links, collapse = ", ")))
    }
  }

  # Validate outcome values based on family
  if (family == "binomial") {
    if (!all(Y %in% c(0, 1))) {
      stop("For family='binomial', Y must be binary (0 or 1)")
    }
  }

  # Define inverse link function
  inverse_link <- function(eta) {
    switch(link,
           "identity" = eta,
           "logit"    = plogis(eta),
           "probit"   = pnorm(eta))
  }

  # ============================================================================
  # 3. Sampling method determination
  # ============================================================================

  # Determine sampling method
  use_external <- !is.null(outcome.samples) && !is.null(ps.samples)
  use_internal <- !is.null(outcome.model) && !is.null(ps.model)

  # Check argument consistency
  if (use_external && use_internal) {
    warning("Both external samples and models provided. Using external samples.")
    use_internal <- FALSE
  } else if (!use_external && !use_internal) {
    stop("Either provide (outcome.model, ps.model) or (outcome.samples, ps.samples)")
  } else if (use_external && (is.null(outcome.samples) || is.null(ps.samples))) {
    stop("Both outcome.samples and ps.samples must be provided for external sampling")
  } else if (use_internal && (is.null(outcome.model) || is.null(ps.model))) {
    stop("Both outcome.model and ps.model must be provided for internal sampling")
  }

  # ============================================================================
  # 4. Design matrices preparation
  # ============================================================================

  # Convert to matrices and add intercepts
  X.lm <- as.matrix(X.lm)
  Z.lm <- cbind(1, X.lm)

  X.ps <- as.matrix(X.ps)
  Z.ps <- cbind(1, X.ps)

  # Create design matrices for potential outcomes using model.matrix()
  # This correctly handles interactions and transformations
  data_A1 <- common_data
  data_A0 <- common_data
  data_A1[[treatment_var_name]] <- 1
  data_A0[[treatment_var_name]] <- 0

  # Generate correct design matrices using the same terms object
  Z.lm1 <- model.matrix(outcome_terms, data = data_A1)
  Z.lm0 <- model.matrix(outcome_terms, data = data_A0)

  nn    <- nrow(X.lm)
  pp.lm <- ncol(X.lm)
  pp.ps <- ncol(X.ps)

  # ============================================================================
  # 5. Obtain posterior samples
  # ============================================================================

  if (use_external) {
    # Use externally provided samples
    cat("Using externally provided posterior samples...\n")

    # Input validation
    if (!is.matrix(outcome.samples) && !is.data.frame(outcome.samples)) {
      stop("outcome.samples must be a matrix or data.frame")
    }
    if (!is.matrix(ps.samples) && !is.data.frame(ps.samples)) {
      stop("ps.samples must be a matrix or data.frame")
    }

    outcome.samples <- as.matrix(outcome.samples)
    ps.samples      <- as.matrix(ps.samples)

    # Check if sample sizes match
    if (nrow(outcome.samples) != nrow(ps.samples)) {
      stop("outcome.samples and ps.samples must have the same number of rows (iterations)")
    }

    # Check dimensions
    expected_otc_cols <- ncol(Z.lm)   # intercept + covariates
    expected_ps_cols  <- ncol(Z.ps)   # intercept + covariates

    if (ncol(outcome.samples) != expected_otc_cols) {
      stop(sprintf("outcome.samples must have %d columns (intercept + %d covariates), but has %d",
                   expected_otc_cols, pp.lm, ncol(outcome.samples)))
    }

    if (ncol(ps.samples) != expected_ps_cols) {
      stop(sprintf("ps.samples must have %d columns (intercept + %d covariates), but has %d",
                   expected_ps_cols, pp.ps, ncol(ps.samples)))
    }

    # Apply burn-in and thinning to external samples
    n_external_samples <- nrow(outcome.samples)

    if (n_external_samples <= (bn + 100)) {
      warning("External samples may be too few for burn-in. Using all samples.")
      betas.otc <- outcome.samples
      betas.ps  <- ps.samples
    } else {
      # Apply burn-in and thinning to external samples
      indices <- seq(bn+1, n_external_samples, by=thin)
      if (length(indices) == 0) {
        indices <- seq(bn+1, n_external_samples)
      }
      betas.otc <- outcome.samples[indices, , drop=FALSE]
      betas.ps  <- ps.samples[indices, , drop=FALSE]
    }

  } else {
    # Generate samples using internal model functions
    cat("Generating posterior samples using provided models...\n")

    # Prepare arguments for outcome model
    outcome.args <- list(Y = Y, X = X.lm, mc = mc)
    outcome.args <- c(outcome.args, outcome.priors)

    # Prepare arguments for propensity score model
    ps.args <- list(Y = A, X = X.ps, mc = mc)
    ps.args <- c(ps.args, ps.priors)

    # Call models with additional arguments
    post.otc <- do.call(outcome.model, outcome.args)
    post.ps  <- do.call(ps.model, ps.args)

    # Apply thinning
    betas.otc <- post.otc[seq(bn+2, mc, by=thin),]
    betas.ps  <- post.ps[seq(bn+2, mc, by=thin),]
  }

  num_iterations <- nrow(betas.otc)

  cat(sprintf("Using %d posterior samples for analysis...\n", num_iterations))

  # ============================================================================
  # 6. Posterior coupling with sequential Monte Carlo
  # ============================================================================

  # Calculate predicted values using inverse link function
  ps.mat <- plogis( tcrossprod(Z.ps, betas.ps) )

  # For outcome model, apply appropriate inverse link function
  eta.mat <- tcrossprod(Z.lm, betas.otc)
  mu.mat  <- inverse_link(eta.mat)

  # Calculate constraint violation (moment condition)
  BB <- colMeans( (A - ps.mat) * (Y - mu.mat)/(ps.mat * (1-ps.mat)) )

  # Initialize SMC parameters
  Inv_mat <- 0
  yy      <- 0
  for(iter in seq_len(num_iterations)){
    BB_vec <- BB[iter]

    exp_val <- exp(BB_vec*0)
    Inv_mat <- Inv_mat + c(exp_val) * (BB_vec %*% t(BB_vec))
    yy      <- yy + c(exp_val) * BB_vec
  }
  lam.new <- 0 - solve(Inv_mat + diag((1e-5), nrow(Inv_mat)), yy)

  # Determine lambda sequence direction
  if(lam.new<0){  lam.s <- (-1)*seq(0, 10, by=0.01) }
  if(lam.new>0){  lam.s <- seq(0, 10, by=0.01) }

  L <- length(lam.s)
  Post.pars.new <- cbind(betas.ps, betas.otc)
  aa <- 0.99  # Smoothing parameter
  BB_mean <- mean(BB)

  # Sequential Monte Carlo loop
  lam.diff <- diff(lam.s)
  for(ll in 1:(L-1)){
    lam.sa <- lam.diff[ll]

    # Calculate importance weights
    val <- lam.sa * BB
    val <- val - max(val)  # Numerical stability
    ww  <- exp(val)/sum(exp(val))

    # Resample particles
    resample_idx <- sample.int(num_iterations, size = num_iterations, replace = TRUE, prob = ww)

    # Apply smoothing kernel
    ep <- mvtnorm::rmvnorm(num_iterations, sigma = (1-aa^2)*stats::cov(Post.pars.new))
    Post.pars.new <- sweep(aa*Post.pars.new[resample_idx,]+ep, 2, (1-aa)*colMeans(Post.pars.new), "+")

    # Extract updated parameters
    betas.otc.new <- Post.pars.new[,-seq_len(pp.ps+1)]
    betas.ps.new  <- Post.pars.new[,seq_len(pp.ps+1)]

    # Recalculate moment condition with updated parameters
    ps.mat.new  <- plogis( tcrossprod(Z.ps, betas.ps.new) )
    eta.mat.new <- tcrossprod(Z.lm, betas.otc.new)
    mu.mat.new  <- inverse_link(eta.mat.new)
    BB          <- colMeans( (A - ps.mat.new) * (Y - mu.mat.new)/(ps.mat.new * (1-ps.mat.new)) )
    BB_mean_new <- mean(BB)

    # Check for convergence (sign change in moment condition)
    sgn <- (BB_mean_new*BB_mean<0)
    if(sgn){ break }
  }

  # ============================================================================
  # 7. Compute posterior of ATE
  # ============================================================================

  # ATE posterior for g-computation
  ww.g <- MCMCpack::rdirichlet(num_iterations, rep(1, nn))

  # Calculate potential outcomes using inverse link function
  eta.ym1.g  <- tcrossprod(betas.otc, Z.lm1)
  eta.ym0.g  <- tcrossprod(betas.otc, Z.lm0)
  post.ym1.g <- inverse_link(eta.ym1.g)
  post.ym0.g <- inverse_link(eta.ym0.g)

  post.mu1.g <- rowSums(ww.g * post.ym1.g)
  post.mu0.g <- rowSums(ww.g * post.ym0.g)
  post.ate.g <- post.mu1.g - post.mu0.g

  # ATE posterior for posterior coupling
  eta.ym1.pc  <- tcrossprod(betas.otc.new, Z.lm1)
  eta.ym0.pc  <- tcrossprod(betas.otc.new, Z.lm0)
  post.ym1.pc <- inverse_link(eta.ym1.pc)
  post.ym0.pc <- inverse_link(eta.ym0.pc)

  post.mu1.pc <- rowMeans(post.ym1.pc)
  post.mu0.pc <- rowMeans(post.ym0.pc)
  post.ate.pc <- post.mu1.pc - post.mu0.pc

  # Return results with family and link information
  result <- list(
    g.comp = post.ate.g,
    pc     = post.ate.pc,
    family = family,
    link   = link,
    call   = call_info,
    data_info = data_info
  )

  class(result) <- c("DRBayes", "list")
  return(result)
}

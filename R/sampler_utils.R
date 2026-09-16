#' Validate a normal prior precision supplied as a scalar or a matrix
#'
#' Returns the precision as a `d` by `d` matrix. The checks are ordered so that
#' a one by one matrix is treated as a matrix rather than as a scalar, and so
#' that an indefinite precision is rejected here rather than surfacing thousands
#' of iterations later as a Cholesky failure inside the sampler.
#'
#' @param prior A positive scalar, a `d` by `d` symmetric positive definite
#'   matrix, or NULL for the weakly informative default.
#' @param d Dimension of the coefficient vector, including the intercept.
#' @param default Precision used when `prior` is NULL.
#' @param arg Argument name, used in error messages.
#'
#' @return A `d` by `d` precision matrix.
#'
#' @keywords internal
#' @noRd
validate_prior_precision <- function(prior, d, default = 1 / 100,
                                     arg = "theta_prior") {
  if (is.null(prior)) {
    return(diag(default, d))
  }

  if (is.matrix(prior)) {
    if (nrow(prior) != d || ncol(prior) != d) {
      stop(arg, " must be a ", d, " by ", d, " matrix, but is ",
           nrow(prior), " by ", ncol(prior))
    }
    if (!isSymmetric(unname(prior), tol = 1e-8)) {
      stop(arg, " must be symmetric")
    }
    if (any(!is.finite(prior))) {
      stop(arg, " must be finite")
    }
    eigenvalues <- eigen(prior, symmetric = TRUE, only.values = TRUE)$values
    if (min(eigenvalues) <= 0) {
      stop(arg, " must be positive definite, but its smallest eigenvalue is ",
           format(min(eigenvalues), digits = 3), ". A precision matrix that is ",
           "not positive definite does not define a proper prior.")
    }
    return(prior)
  }

  if (length(prior) != 1L || !is.numeric(prior) || !is.finite(prior)) {
    stop(arg, " must be NULL, a single positive number, or a ", d, " by ", d,
         " positive definite matrix")
  }
  if (prior <= 0) {
    stop(arg, " is a prior PRECISION, so it must be positive; ", prior,
         " would give an improper prior. Pass 1/variance.")
  }
  diag(prior, d)
}


#' Columns of a design matrix that involve a given variable
#'
#' Finds every column derived from `variable`, including the ones contributed by
#' interactions such as `A:X1`. Used to tell a shrinkage prior which coefficients
#' must not be shrunk: in a causal model that is the treatment main effect
#' together with every effect modification term, since shrinking those attenuates
#' the treatment effect.
#'
#' @param terms_obj Terms object of the fitted model frame.
#' @param Z Design matrix built from `terms_obj`, carrying its `assign` attribute.
#' @param variable Name of the variable to look for.
#'
#' @return Integer column indices into `Z`.
#'
#' @keywords internal
#' @noRd
columns_involving <- function(terms_obj, Z, variable) {
  factors <- attr(terms_obj, "factors")
  if (is.null(factors) || !(variable %in% rownames(factors))) {
    return(integer(0))
  }
  involving <- which(factors[variable, ] != 0)
  assign    <- attr(Z, "assign")
  if (is.null(assign)) {
    return(integer(0))
  }
  which(assign %in% involving)
}


#' Storage array for posterior draws
#'
#' Allocates as double rather than letting `array(NA, ...)` create a logical
#' array that is copied and coerced on the first assignment, and keeps the
#' column names the caller supplied instead of replacing them with X1, X2, ...
#' The iteration by chain by parameter layout is what the posterior package
#' calls a draws_array, so the result can be passed to convergence diagnostics
#' without being rearranged first.
#'
#' @param mc Number of iterations.
#' @param chains Number of chains.
#' @param X The covariate matrix the sampler was given.
#'
#' @keywords internal
#' @noRd
draw_storage <- function(mc, chains, X) {
  pp <- ncol(X)
  nm <- colnames(X)
  if (is.null(nm) || !all(nzchar(nm))) {
    # paste0("X", integer(0)) is "X", not character(0), so an intercept-only
    # design would otherwise be handed one name too many.
    nm <- if (pp > 0L) paste0("X", seq_len(pp)) else character(0)
  }
  array(NA_real_, dim = c(mc, chains, pp + 1L),
        dimnames = list(NULL, NULL, c("(Intercept)", nm)))
}


#' Validate the requested number of chains
#'
#' @param chains Number of Markov chains to run.
#'
#' @return The number of chains, as an integer.
#'
#' @keywords internal
#' @noRd
validate_chains <- function(chains) {
  # as.integer() would turn anything past the integer bound into NA with only
  # a warning, so the bound is part of the contract rather than left to it.
  if (!is.numeric(chains) || length(chains) != 1L || !is.finite(chains) ||
      chains != round(chains) || chains < 1L ||
      chains > .Machine$integer.max) {
    stop("chains must be a single positive integer no larger than ",
         .Machine$integer.max, ", the number of Markov chains to run. Four or ",
         "more are recommended: R-hat compares chains against each other, so ",
         "one chain cannot support it.")
  }
  as.integer(chains)
}


#' Marginal prior standard deviations implied by a precision matrix
#'
#' @param prior A positive definite precision matrix.
#'
#' @keywords internal
#' @noRd
prior_sd_from_precision <- function(prior) {
  sqrt(diag(chol2inv(chol(prior))))
}


#' Quick frequentist estimate used to place the starting values
#'
#' Returns the coefficients and their standard errors, or NULL when the fit is
#' unusable. A rank deficient design, a logistic fit that has not converged and
#' complete separation all end up here, and the caller answers them by falling
#' back on the prior rather than by failing: an initialisation problem must not
#' stop a sampler that would have run perfectly well from a poorer start.
#'
#' @param Y Response passed to the sampler.
#' @param Z Design matrix including its intercept column.
#' @param family Either "gaussian" or "binomial".
#' @param link Link for the binomial fit, "logit" or "probit".
#'
#' @return A list with `coefficients` and `se`, or NULL.
#'
#' @keywords internal
#' @noRd
quick_coef_fit <- function(Y, Z, family, link) {
  fit <- tryCatch(
    # Warnings such as "fitted probabilities numerically 0 or 1" come from a
    # fit the user never asked for, so they are not passed on. A fit that is
    # actually unusable is caught by the checks below instead.
    suppressWarnings(
      if (family == "gaussian") {
        stats::lm.fit(Z, Y)
      } else {
        stats::glm.fit(Z, Y, family = stats::binomial(link = link))
      }),
    error = function(e) NULL)

  d <- ncol(Z)
  if (is.null(fit) || !isTRUE(fit$rank == d)) {
    return(NULL)
  }
  if (family == "binomial" && !isTRUE(fit$converged)) {
    return(NULL)
  }
  coefficients <- unname(fit$coefficients)
  if (!all(is.finite(coefficients))) {
    return(NULL)
  }

  unscaled <- tryCatch(chol2inv(qr.R(fit$qr)), error = function(e) NULL)
  if (is.null(unscaled)) {
    return(NULL)
  }
  se <- numeric(d)
  # qr.R() is in pivoted order, so the variances belong to the pivoted columns.
  se[fit$qr$pivot] <- sqrt(pmax(diag(unscaled), 0))
  if (family == "gaussian") {
    se <- se * sqrt(sum(fit$residuals^2) / max(length(Y) - d, 1L))
  }
  # Under separation the coefficients stay finite while the standard errors
  # blow up, and a starting value drawn on that scale is no better than noise.
  if (!all(is.finite(se)) || any(se <= 0) || max(se) > 1e6) {
    return(NULL)
  }
  list(coefficients = coefficients, se = se)
}


#' Starting values for each chain
#'
#' Chains that all start at the same point cannot reveal multimodality or
#' sensitivity to the starting value, so a convergence diagnostic computed from
#' them understates the problem (Vehtari et al. 2021, section 2). The automatic
#' starting values are therefore overdispersed relative to the posterior: chain
#' one sits at the quick frequentist estimate and every other chain is drawn
#' around it with a standard deviation of `spread` times that estimate's own
#' standard error, which is the posterior scale up to the usual asymptotic
#' approximation.
#'
#' Random numbers are taken from the caller's stream, so the starting values
#' are reproducible under a single `set.seed()` and no seed is set here.
#'
#' @param init NULL, or a list of `chains` numeric vectors of length `ncol(Z)`.
#' @param chains Number of chains.
#' @param Y Response passed to the sampler.
#' @param Z Design matrix including its intercept column.
#' @param family Either "gaussian" or "binomial".
#' @param link Link for the binomial fit, "logit" or "probit".
#' @param prior_sd Prior standard deviation of each coefficient, used when the
#'   frequentist fit is unusable. A scalar or a vector of length `ncol(Z)`.
#' @param spread Multiple of the standard error used to disperse the chains.
#'
#' @return A list of `chains` numeric vectors of length `ncol(Z)`.
#'
#' @keywords internal
#' @noRd
resolve_inits <- function(init, chains, Y, Z, family = "gaussian",
                          link = "logit", prior_sd = 10, spread = 4) {
  d <- ncol(Z)

  if (!is.null(init)) {
    if (!is.list(init)) {
      stop("init must be NULL or a list holding one starting value per chain, ",
           "but is of class ", class(init)[1L], ". A single starting value ",
           "still has to be wrapped in a list, as in list(rep(0, ", d, ")).")
    }
    if (length(init) != chains) {
      stop("init has ", length(init), " element(s) but chains is ", chains,
           ", so ", chains, " starting values are needed, one per chain")
    }
    for (i in seq_len(chains)) {
      value <- init[[i]]
      if (!is.numeric(value) || length(value) != d) {
        stop("init[[", i, "]] must be a numeric vector of length ", d,
             ", the intercept followed by the ", d - 1L, " column(s) of X, ",
             "but is of class ", class(value)[1L], " and length ",
             length(value))
      }
      if (!all(is.finite(value))) {
        stop("init[[", i, "]] must be finite; found ", sum(!is.finite(value)),
             " missing or infinite value(s)")
      }
    }
    return(lapply(init, function(value) as.numeric(unname(value))))
  }

  fit    <- quick_coef_fit(Y, Z, family, link)
  centre <- if (is.null(fit)) numeric(d) else fit$coefficients
  scale  <- if (is.null(fit)) rep_len(prior_sd, d) else spread * fit$se

  inits <- vector("list", chains)
  # Keeping chain one at the estimate itself guarantees that one chain starts
  # somewhere the data support, however wide the dispersion of the others.
  inits[[1L]] <- centre
  for (chain in seq_len(chains)[-1L]) {
    inits[[chain]] <- centre + scale * stats::rnorm(d)
  }
  inits
}


#' Shared input checks for the Gibbs samplers
#'
#' @return The covariate matrix, coerced and validated.
#'
#' @keywords internal
#' @noRd
validate_sampler_inputs <- function(Y, X, mc, binary = FALSE) {
  if (!is.numeric(Y)) {
    stop("Y must be a numeric vector")
  }
  if (binary && !all(Y %in% c(0, 1))) {
    stop("Y must be binary (0 or 1)")
  }
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data.frame")
  }
  X <- as.matrix(X)
  if (!is.numeric(X)) {
    stop("X must be numeric. Expand factors into dummy columns first, for ",
         "example with model.matrix().")
  }
  # Not just NA: a single NaN or Inf reaching pgdraw::pgdraw() sends its
  # rejection sampler into a loop that never returns, so the logistic samplers
  # hang rather than fail. Refuse non-finite input here, before that can happen.
  if (!all(is.finite(X))) {
    stop("X must be finite; found ", sum(!is.finite(X)),
         " missing or infinite value(s)")
  }
  if (!all(is.finite(Y))) {
    stop("Y must be finite; found ", sum(!is.finite(Y)),
         " missing or infinite value(s)")
  }
  # An intercept-only model, ncol(X) == 0, is legitimate: it is what
  # ps.formula = A ~ 1 asks for. draw_storage handles it via seq_len().
  if (nrow(X) != length(Y)) {
    stop("X has ", nrow(X), " rows but Y has ", length(Y), " elements")
  }
  if (!is.numeric(mc) || length(mc) != 1L || mc != round(mc) || mc < 1L) {
    stop("mc must be a single positive integer")
  }
  X
}


#' Validate the set of coefficients exempted from shrinkage
#'
#' @param unshrunk Integer column indices into `X`, or a character vector of
#'   column names.
#' @param pp Number of columns of `X`.
#'
#' @return Integer column indices, sorted and deduplicated.
#'
#' @keywords internal
#' @noRd
validate_unshrunk <- function(unshrunk, pp) {
  if (is.null(unshrunk) || length(unshrunk) == 0L) {
    return(integer(0))
  }
  if (is.character(unshrunk)) {
    stop("unshrunk must be given as column indices, not names")
  }
  if (!is.numeric(unshrunk) || any(unshrunk != round(unshrunk))) {
    stop("unshrunk must be whole numbers indexing the columns of X")
  }
  unshrunk <- sort(unique(as.integer(unshrunk)))
  if (any(unshrunk < 1L) || any(unshrunk > pp)) {
    stop("unshrunk indexes columns outside X, which has ", pp, " columns: ",
         toString(unshrunk[unshrunk < 1L | unshrunk > pp]))
  }
  unshrunk
}


#' Prior scale for the global shrinkage parameter
#'
#' Piironen and Vehtari (2017) recommend centring the prior on tau at
#' \eqn{\tau_0 = \frac{p_0}{p - p_0} \frac{\sigma}{\sqrt{n}}}, where \eqn{p_0}
#' is a guess at the number of non-zero coefficients and \eqn{\sigma} is the
#' error scale of the observation model. When `sigma` is NULL the caller's model
#' supplies its own scale elsewhere and only the dimensional part is used.
#'
#' @keywords internal
#' @noRd
resolve_tau_prior <- function(tau_prior, p, n, p0, sigma) {
  if (!is.null(tau_prior)) {
    if (!is.numeric(tau_prior) || length(tau_prior) != 1L || tau_prior <= 0 ||
        !is.finite(tau_prior)) {
      stop("tau_prior must be a single positive number")
    }
    return(tau_prior)
  }
  if (!is.numeric(p0) || length(p0) != 1L || p0 < 1 || p0 != round(p0)) {
    stop("p0 must be a single positive integer, a guess at how many ",
         "coefficients are non-zero")
  }
  if (p <= p0) {
    # Fewer coefficients than the guessed number of signals leaves nothing to
    # shrink towards, so fall back to the unit-scale half-Cauchy.
    return(1)
  }
  scale <- if (is.null(sigma)) 1 else sigma
  p0 / (p - p0) * scale / sqrt(n)
}


#' Draw from a multivariate normal given its precision matrix
#'
#' Uses the Cholesky factor of the precision, so no explicit inverse is formed.
#' A failure here means the precision is not positive definite, which is
#' actionable information the raw chol() error does not convey.
#'
#' @param Q Precision matrix.
#' @param b Precision times mean.
#' @param iter Iteration number, used in the error message.
#'
#' @keywords internal
#' @noRd
draw_normal_precision <- function(Q, b, iter) {
  U <- tryCatch(
    chol(Q),
    error = function(e) {
      stop("The coefficient precision matrix stopped being positive definite ",
           "at iteration ", iter, ". This usually means a shrinkage parameter ",
           "underflowed or the design matrix is collinear. Original message: ",
           conditionMessage(e), call. = FALSE)
    })
  d <- ncol(Q)
  drop(backsolve(U, backsolve(U, b, transpose = TRUE) + stats::rnorm(d)))
}


#' Evaluate an expression under a seed without disturbing the caller
#'
#' Sets `seed`, evaluates `expr`, then puts the random number generator back
#' exactly as it was, including its kind and the case where the caller had never
#' drawn a random number at all. A function that simply calls [set.seed()] and
#' returns would replace the caller's stream, so a script that seeds itself at
#' the top would lose control of everything after the call.
#'
#' @param seed A single number.
#' @param expr Expression to evaluate.
#'
#' @keywords internal
#' @noRd
with_preserved_rng <- function(seed, expr) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
            add = TRUE)
  }
  old_kind <- RNGkind()
  # Restore the kind LAST, because switching kind re-seeds the generator.
  #
  # How many components there are is not fixed. R has had three since 3.6.0 and
  # R-devel added a fourth in 2026, for the binomial method, so naming them
  # positionally would quietly stop restoring the last one on a newer R.
  # Whatever this build reports is what goes back.
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE, after = FALSE)

  set.seed(seed)
  expr
}


#' Default for a NULL value
#'
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

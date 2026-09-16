#' Bayesian bootstrap doubly robust estimator
#'
#' The estimator of Saarela et al. (2016), as displayed in section 4.3 of
#' Orihara, Momozaki and Sugasawa (2025). For each Dirichlet draw \eqn{\xi} it
#' refits both models with those weights and evaluates
#' \deqn{\sum_i \xi_i \left\{ m_1(X_i) - m_0(X_i) +
#'   \frac{A_i - e(X_i)}{e(X_i)(1 - e(X_i))} (Y_i - m_{A_i}(X_i)) \right\}.}
#'
#' Both regressions are solved on design matrices built once, outside the loop.
#' Refitting through the formula interface would rebuild the terms object and
#' the model matrix on every draw from data that never changes, and would leave
#' the fitted coefficients to be matched against a hand-built design matrix by
#' position, which breaks silently on factors and on rank-deficient fits.
#'
#' @param md Output of [prepare_model_data()].
#' @param num_iterations Number of Dirichlet draws.
#' @param family,link Outcome model family and link.
#' @param ps.link Propensity score link.
#' @param trim Propensity scores are confined to `[trim, 1 - trim]`. Zero leaves
#'   them untouched.
#' @param verbose Whether to report progress.
#'
#' @return Numeric vector of `num_iterations` posterior draws of the estimand.
#'
#' @keywords internal
#' @noRd
bayes_bootstrap_dr <- function(md, num_iterations, family, link, ps.link,
                               trim = 0, verbose = TRUE) {

  Y  <- md$Y
  A  <- md$A
  nn <- length(Y)

  inv_link    <- link_function(link)
  ps_inv_link <- link_function(ps.link)

  ps_family      <- stats::quasibinomial(link = ps.link)
  outcome_family <- if (family == "binomial") {
    stats::quasibinomial(link = link)
  } else {
    NULL
  }

  post.ate    <- numeric(num_iterations)
  n_truncated <- 0L

  for (iter in seq_len(num_iterations)) {
    # Bayesian bootstrap weights. Scaling by nn does not move either fit, since
    # both are weighted least squares problems, but it keeps the weights on the
    # scale of counts that glm.fit's convergence test expects.
    xi <- dirichlet_weights(nn)
    w  <- nn * xi

    alpha <- fit_weighted(md$Z.ps, A, w, ps_family, "propensity score")
    beta  <- fit_weighted(md$Z.lm, Y, w, outcome_family, "outcome")

    ps <- ps_inv_link(drop(md$Z.ps %*% alpha))
    if (trim > 0) {
      clipped     <- ps < trim | ps > 1 - trim
      n_truncated <- n_truncated + sum(clipped)
      ps          <- pmin(pmax(ps, trim), 1 - trim)
    }
    if (any(!is.finite(ps)) || any(ps <= 0 | ps >= 1)) {
      stop("The propensity score reached 0 or 1 at Dirichlet draw ", iter,
           ", so the augmentation term is not finite. This is a positivity ",
           "violation: check the overlap between the treatment groups, or set ",
           "trim to confine the propensity score away from 0 and 1.")
    }

    mu <- inv_link(drop(md$Z.lm  %*% beta))
    y1 <- inv_link(drop(md$Z.lm1 %*% beta))
    y0 <- inv_link(drop(md$Z.lm0 %*% beta))

    # The augmentation multiplier of equation (2.3). For binary A this is the
    # same quantity as A/e - (1 - A)/(1 - e).
    cc <- (A - ps) / (ps * (1 - ps))

    post.ate[iter] <- sum(xi * (y1 - y0)) + sum(xi * (Y - mu) * cc)
  }

  if (n_truncated > 0 && verbose) {
    message(n_truncated, " of ", nn * num_iterations,
            " fitted propensity scores were truncated to [", trim, ", ",
            1 - trim, "].")
  }

  post.ate
}


#' Draw one vector of Bayesian bootstrap weights
#'
#' A Dirichlet(1, ..., 1) vector on the simplex of dimension `nn`. Independent
#' unit exponentials divided by their sum have exactly that distribution, since
#' Exp(1) is Gamma(1, 1) and a Dirichlet vector is a vector of independent
#' gammas with the corresponding shapes divided by its own sum. Drawing it this
#' way keeps the generator in stats and avoids the matrix allocation and the
#' matrix product that a general Dirichlet routine spends on a single draw.
#'
#' @param nn Length of the weight vector.
#'
#' @return Numeric vector of `nn` positive weights summing to one.
#'
#' @keywords internal
#' @noRd
dirichlet_weights <- function(nn) {
  g <- stats::rexp(nn)
  g / sum(g)
}


#' Solve one weighted regression on a prebuilt design matrix
#'
#' @param Z Design matrix including its intercept column.
#' @param y Response.
#' @param w Weights.
#' @param family A family object, or NULL for a Gaussian fit.
#' @param what Name used in error messages.
#'
#' @keywords internal
#' @noRd
fit_weighted <- function(Z, y, w, family, what) {
  if (is.null(family)) {
    # The Gaussian score is linear in the coefficients, so this is exact. The
    # pivoted QR also reports rank deficiency instead of returning silent NAs.
    s   <- sqrt(w)
    fit <- stats::.lm.fit(s * Z, s * y)
    if (fit$rank < ncol(Z)) {
      stop(rank_deficiency_message(Z, fit$rank, what))
    }
    return(fit$coefficients)
  }

  fit <- stats::glm.fit(Z, y, weights = w, family = family)
  if (fit$rank < ncol(Z) || anyNA(fit$coefficients)) {
    stop(rank_deficiency_message(Z, fit$rank, what))
  }
  if (!fit$converged) {
    stop("The ", what, " model did not converge on a Bayesian bootstrap draw. ",
         "This usually means the treatment groups are nearly separated by the ",
         "covariates in that model.")
  }
  unname(fit$coefficients)
}


#' @keywords internal
#' @noRd
rank_deficiency_message <- function(Z, rank, what) {
  paste0("The ", what, " design matrix is rank deficient: ", ncol(Z),
         " columns but rank ", rank, ". Some coefficients are not estimable, ",
         "so the fitted values would be arbitrary. Drop the collinear columns ",
         "from the formula.")
}


#' @keywords internal
#' @noRd
link_function <- function(link) {
  switch(link,
         "identity" = function(eta) eta,
         "logit"    = stats::plogis,
         "probit"   = stats::pnorm,
         stop("Unsupported link: ", link))
}


#' Resolve the outcome family and link, and the propensity score link
#'
#' @keywords internal
#' @noRd
resolve_family_link <- function(Y, family, link, ps.link, verbose = TRUE) {
  if (is.null(family)) {
    family <- if (all(Y %in% c(0, 1))) "binomial" else "gaussian"
    if (verbose) {
      message("Auto-detected ",
              if (family == "binomial") "binary" else "continuous",
              " outcome: using family = \"", family, "\".")
    }
  } else {
    family <- match.arg(family, c("gaussian", "binomial"))
  }

  valid_links <- switch(family,
                        "gaussian" = "identity",
                        "binomial" = c("logit", "probit"))
  if (is.null(link)) {
    link <- valid_links[1L]
  } else if (!link %in% valid_links) {
    stop("Link \"", link, "\" is not supported for family \"", family,
         "\". Valid links: ", paste(valid_links, collapse = ", "))
  }

  if (family == "binomial" && !all(Y %in% c(0, 1))) {
    stop("For family = \"binomial\", the outcome must be binary (0 or 1)")
  }

  list(family  = family,
       link    = link,
       ps.link = match.arg(ps.link, c("logit", "probit")))
}

#' Build the model data shared by the causal estimators
#'
#' Turns a pair of formulas and a data frame into the response vectors, design
#' matrices and counterfactual design matrices that both [drbayes_pc()] and the
#' Bayesian bootstrap estimators need. Keeping this in one place is what stops
#' the two estimators from disagreeing about how a formula is interpreted.
#'
#' The counterfactual design matrices are built through the *fitted* terms
#' object, carrying its `predvars`, `xlev` and `contrasts`. Rebuilding them from
#' the raw data instead would re-derive data-dependent bases such as
#' [stats::poly()] and re-map factor levels, so the counterfactual columns would
#' no longer be on the same scale as the coefficients they multiply.
#'
#' @param outcome.formula Formula for the outcome model, e.g. `Y ~ A + X1`.
#' @param ps.formula Formula for the propensity score model, e.g. `A ~ X1`.
#' @param data Data frame holding every variable used by either formula.
#' @param na.action One of `"na.omit"`, `"na.fail"` or `"na.exclude"`.
#'
#' @return A list with the response vectors `Y` and `A`; the design matrices
#'   `Z.lm` and `Z.ps` including their intercept column and `X.lm` and `X.ps`
#'   without it; the counterfactual outcome design matrices `Z.lm1` and `Z.lm0`;
#'   the treatment variable name `treatment`; and `data_info`.
#'
#' @keywords internal
#' @noRd
prepare_model_data <- function(outcome.formula, ps.formula, data,
                               na.action = "na.omit") {

  if (!inherits(outcome.formula, "formula")) {
    stop("outcome.formula must be a formula object (e.g., Y ~ A + X1 + X2)")
  }
  if (!inherits(ps.formula, "formula")) {
    stop("ps.formula must be a formula object (e.g., A ~ X1 + X2)")
  }
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }

  na.action     <- match.arg(na.action, c("na.omit", "na.fail", "na.exclude"))
  na.action.fun <- switch(na.action,
                          "na.omit"    = stats::na.omit,
                          "na.fail"    = stats::na.fail,
                          "na.exclude" = stats::na.exclude)

  # Expand any "." against the data before looking at the variables, otherwise
  # a formula such as Y ~ A + . reports "." itself as a missing variable.
  outcome.formula <- stats::formula(stats::terms(outcome.formula, data = data))
  ps.formula      <- stats::formula(stats::terms(ps.formula, data = data))

  all_vars     <- unique(c(all.vars(outcome.formula), all.vars(ps.formula)))
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("The following variables are not found in data: ",
         paste(missing_vars, collapse = ", "))
  }

  # Both models must be fitted on the same rows. Missing values are resolved by
  # intersecting row NAMES; using the names as positions would silently select
  # different observations whenever the row names are not 1:nrow(data), which is
  # the normal state of a data frame that has been through subset() or merge().
  outcome_frame <- stats::model.frame(outcome.formula, data = data,
                                      na.action = na.action.fun)
  ps_frame      <- stats::model.frame(ps.formula, data = data,
                                      na.action = na.action.fun)

  if (nrow(outcome_frame) != nrow(ps_frame)) {
    warning("The outcome and propensity score models have different missing ",
            "value patterns; only the ", length(intersect(rownames(outcome_frame),
                                                          rownames(ps_frame))),
            " observations complete in both are used.")
  }

  common_rows <- intersect(rownames(outcome_frame), rownames(ps_frame))
  if (length(common_rows) == 0) {
    stop("No observation is complete in both the outcome and the propensity ",
         "score model")
  }

  common_data   <- data[common_rows, , drop = FALSE]
  outcome_frame <- stats::model.frame(outcome.formula, data = common_data,
                                      na.action = stats::na.fail)
  ps_frame      <- stats::model.frame(ps.formula, data = common_data,
                                      na.action = stats::na.fail)

  # Use the terms attached to the model frames: unlike terms(formula) these
  # carry predvars, which is what pins poly(), scale() and ns() to the basis
  # that was fitted rather than one re-derived from whatever data is handed in.
  outcome_terms <- stats::terms(outcome_frame)
  ps_terms      <- stats::terms(ps_frame)

  for (nm in c("outcome", "ps")) {
    tt <- if (nm == "outcome") outcome_terms else ps_terms
    if (!identical(attr(tt, "intercept"), 1L)) {
      stop(nm, ".formula must include an intercept; posterior coupling ",
           "parameterises both models with one. Remove the '- 1' or '+ 0' term.")
    }
  }

  Y         <- stats::model.response(outcome_frame)
  treatment <- all.vars(ps.formula)[1]
  A         <- encode_treatment(stats::model.response(ps_frame), treatment)

  Z.lm <- stats::model.matrix(outcome_terms, outcome_frame)
  Z.ps <- stats::model.matrix(ps_terms, ps_frame)

  # Counterfactual design matrices: set the treatment to 1 and to 0 in the data
  # and rebuild through the fitted terms, so that every column derived from the
  # treatment, including interactions such as A:X1, is recomputed.
  Z.lm1 <- counterfactual_matrix(outcome_terms, outcome_frame, common_data,
                                 treatment, Z.lm, value = 1)
  Z.lm0 <- counterfactual_matrix(outcome_terms, outcome_frame, common_data,
                                 treatment, Z.lm, value = 0)

  list(
    Y         = Y,
    A         = A,
    Z.lm      = Z.lm,
    Z.ps      = Z.ps,
    X.lm      = Z.lm[, -1, drop = FALSE],
    X.ps      = Z.ps[, -1, drop = FALSE],
    Z.lm1     = Z.lm1,
    Z.lm0     = Z.lm0,
    treatment = treatment,
    # Columns of X.lm derived from the treatment: its main effect and every
    # interaction with it. A shrinkage prior must leave these alone, since
    # shrinking them attenuates the treatment effect.
    treatment_columns = columns_involving(outcome_terms, Z.lm, treatment) - 1L,
    data_info = list(
      n_observations       = length(Y),
      n_treated            = sum(A),
      n_control            = sum(1 - A),
      n_outcome_terms      = ncol(Z.lm) - 1L,
      n_ps_terms           = ncol(Z.ps) - 1L,
      outcome_formula      = outcome.formula,
      ps_formula           = ps.formula,
      missing_observations = nrow(data) - length(Y)
    )
  )
}


#' Coerce a treatment variable to a 0/1 integer vector
#'
#' Accepts numeric 0/1, logical, or a two-level factor or character vector. For
#' the latter the second level counts as treated, matching how
#' [stats::model.matrix()] dummy-codes it.
#'
#' @keywords internal
#' @noRd
encode_treatment <- function(a, name) {
  if (is.factor(a) || is.character(a)) {
    a      <- as.factor(a)
    levs   <- levels(a)
    if (length(levs) != 2L) {
      stop("Treatment variable (", name, ") must have exactly two levels, but ",
           "has ", length(levs), ": ", toString(levs))
    }
    message("Treating '", levs[2L], "' as treated and '", levs[1L],
            "' as control for ", name, ".")
    return(as.integer(a == levs[2L]))
  }
  if (is.logical(a)) {
    return(as.integer(a))
  }
  if (!is.numeric(a) || !all(a %in% c(0, 1))) {
    stop("Treatment variable (", name, ") must be binary. Supply 0/1, a ",
         "logical, or a two-level factor. Found: ",
         toString(utils::head(unique(a), 5)))
  }
  as.integer(a)
}


#' Rebuild an outcome design matrix with the treatment fixed at one value
#'
#' @keywords internal
#' @noRd
counterfactual_matrix <- function(fitted_terms, fitted_frame, data, treatment,
                                  Z, value) {
  cf_data <- data
  original <- data[[treatment]]
  cf_data[[treatment]] <- if (is.factor(original)) {
    factor(levels(original)[value + 1L], levels = levels(original))
  } else if (is.character(original)) {
    sort(unique(original))[value + 1L]
  } else if (is.logical(original)) {
    as.logical(value)
  } else {
    value
  }

  rhs <- stats::delete.response(fitted_terms)
  mf  <- stats::model.frame(rhs, cf_data,
                            xlev      = stats::.getXlevels(fitted_terms, fitted_frame),
                            na.action = stats::na.fail)
  Z_cf <- stats::model.matrix(rhs, mf, contrasts.arg = attr(Z, "contrasts"))

  if (!identical(colnames(Z_cf), colnames(Z))) {
    stop("The counterfactual design matrix does not match the fitted one. ",
         "Setting ", treatment, " to ", value, " produced columns ",
         toString(setdiff(colnames(Z_cf), colnames(Z))), " and dropped ",
         toString(setdiff(colnames(Z), colnames(Z_cf))), ".")
  }
  Z_cf
}

#' Functions renamed in DRBayes 0.1.0
#'
#' @description
#' The exported functions were renamed before the first CRAN release, to one
#' naming scheme and away from dots. A dot in a function name is how R marks an
#' S3 method, so `print.DRBayes` and a function called `DRBayes.PC` would have
#' sat in the same namespace with no way to tell them apart.
#'
#' Each old name still works and forwards to its replacement, with a warning.
#' They will be removed in the next release.
#'
#' \tabular{ll}{
#'   \strong{Old} \tab \strong{New} \cr
#'   `DRBayes.PC()`  \tab [drbayes_pc()] \cr
#'   `DRBayes.BB()`  \tab [drbayes_bb()] \cr
#'   `DRBayes.BB.nleqslv()` \tab [drbayes_bb()] \cr
#'   `B.LM()`        \tab [bayes_lm()] \cr
#'   `B.Logit()`     \tab [bayes_logit()] \cr
#'   `B.Probit()`    \tab [bayes_probit()] \cr
#'   `HS.LM()`       \tab [bayes_lm_hs()] \cr
#'   `HS.Logit()`    \tab [bayes_logit_hs()] \cr
#'   `HS.Probit()`   \tab [bayes_probit_hs()]
#' }
#'
#' `DRBayes.BB.nleqslv()` solved the same two weighted regressions with a
#' general nonlinear root finder. It computed the same estimator as
#' `DRBayes.BB()` to about 1e-08, the difference being only the root finder's
#' stopping tolerance, so it forwards there rather than being reimplemented.
#'
#' @param ... Passed to the replacement function.
#'
#' @return Whatever the replacement function returns.
#'
#' @name DRBayes-deprecated
NULL

#' @rdname DRBayes-deprecated
#' @export
DRBayes.PC <- function(...) {
  .Deprecated("drbayes_pc")
  drbayes_pc(...)
}

#' @rdname DRBayes-deprecated
#' @export
DRBayes.BB <- function(...) {
  .Deprecated("drbayes_bb")
  drbayes_bb(...)
}

#' @rdname DRBayes-deprecated
#' @export
DRBayes.BB.nleqslv <- function(...) {
  .Deprecated("drbayes_bb")
  drbayes_bb(...)
}

#' @rdname DRBayes-deprecated
#' @export
B.LM <- function(...) {
  .Deprecated("bayes_lm")
  bayes_lm(...)
}

#' @rdname DRBayes-deprecated
#' @export
B.Logit <- function(...) {
  .Deprecated("bayes_logit")
  bayes_logit(...)
}

#' @rdname DRBayes-deprecated
#' @export
B.Probit <- function(...) {
  .Deprecated("bayes_probit")
  bayes_probit(...)
}

#' @rdname DRBayes-deprecated
#' @export
HS.LM <- function(...) {
  .Deprecated("bayes_lm_hs")
  bayes_lm_hs(...)
}

#' @rdname DRBayes-deprecated
#' @export
HS.Logit <- function(...) {
  .Deprecated("bayes_logit_hs")
  bayes_logit_hs(...)
}

#' @rdname DRBayes-deprecated
#' @export
HS.Probit <- function(...) {
  .Deprecated("bayes_probit_hs")
  bayes_probit_hs(...)
}

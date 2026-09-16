#' Choose the sampler for a model from its family, link and prior
#'
#' The model to fit is already fully described by the formula, the family and
#' the link, and given those the sampling algorithm follows: a Gaussian
#' likelihood with an identity link has a conjugate Gibbs sampler, a logistic
#' one is handled by Polya-Gamma augmentation, and a probit one by the Albert
#' and Chib latent variable scheme. The only choice left to the user is the
#' prior, which is why that, and not a function, is what they pass.
#'
#' @param link One of `"identity"`, `"logit"` or `"probit"`.
#' @param prior One of `"normal"` or `"horseshoe"`.
#'
#' @return The sampler function.
#'
#' @keywords internal
#' @noRd
select_sampler <- function(link, prior) {
  prior <- match.arg(prior, c("normal", "horseshoe"))
  samplers <- list(
    normal    = list(identity = bayes_lm,    logit = bayes_logit,
                     probit   = bayes_probit),
    horseshoe = list(identity = bayes_lm_hs, logit = bayes_logit_hs,
                     probit   = bayes_probit_hs))
  sampler <- samplers[[prior]][[link]]
  if (is.null(sampler)) {
    stop("No sampler for a ", link, " link with a ", prior, " prior")
  }
  sampler
}


#' Name of the sampler chosen, for reporting back to the user
#'
#' @keywords internal
#' @noRd
sampler_name <- function(link, prior) {
  base <- switch(link, identity = "bayes_lm", logit = "bayes_logit",
                 probit = "bayes_probit")
  if (prior == "horseshoe") paste0(base, "_hs") else base
}


#' The link a built-in sampler draws its coefficients on
#'
#' Used to refuse a combination that would evaluate draws with the wrong inverse
#' link, which is silent and produces a plausible but wrong answer.
#'
#' @keywords internal
#' @noRd
sampler_link <- function(model) {
  if (is.null(model)) return(NULL)
  known <- list(
    identity = list(bayes_lm, bayes_lm_hs),
    logit    = list(bayes_logit, bayes_logit_hs),
    probit   = list(bayes_probit, bayes_probit_hs))
  for (lk in names(known)) {
    if (any(vapply(known[[lk]], identical, logical(1), model))) {
      return(lk)
    }
  }
  NULL
}

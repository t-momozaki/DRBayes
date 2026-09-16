#' Fit One of the Package Models with Stan Instead of the Gibbs Sampler
#'
#' @description
#' Fits the same six models the Gibbs samplers fit, using the No-U-Turn sampler
#' through CmdStan, and returns the draws in the layout the Gibbs samplers
#' return. The result can be passed straight to the \code{outcome.samples} or
#' \code{ps.samples} argument of \code{\link{drbayes_pc}}, which does not care
#' which backend produced it.
#'
#' This is an option, never the default. See the section "Why Stan is not the
#' default" below.
#'
#' @param model Character, which model to fit. One of "lm", "logit", "probit",
#'   "lm_hs", "logit_hs" and "probit_hs", naming the counterparts of
#'   \code{\link{bayes_lm}}, \code{\link{bayes_logit}},
#'   \code{\link{bayes_probit}}, \code{\link{bayes_lm_hs}},
#'   \code{\link{bayes_logit_hs}} and \code{\link{bayes_probit_hs}}.
#' @param Y A numeric vector of outcomes. Continuous for "lm" and "lm_hs",
#'   binary for the other four.
#' @param X A matrix or data.frame of covariates, one row per observation. An
#'   intercept term is added automatically, as in the Gibbs samplers.
#' @param mc A positive integer, the number of draws kept per chain after
#'   warm-up (default: 5000).
#' @param chains A positive integer, the number of Markov chains (default: 4).
#' @param warmup A positive integer, the number of warm-up iterations per chain
#'   (default: 1000). Warm-up draws are discarded by the sampler and never
#'   appear in the returned array, so there is no burn-in left to remove. When
#'   the result is handed to \code{\link{drbayes_pc}}, set that function's
#'   \code{bn} to 0 unless you want to throw good draws away.
#' @param theta_prior Prior precision of the coefficients under the normal
#'   prior, a positive scalar, a square matrix, or NULL for the weakly
#'   informative default of 1/100 times the identity. Used by "lm", "logit" and
#'   "probit".
#' @param sigma_prior A numeric vector of length 2, the shape and scale of the
#'   inverse gamma prior on the error variance (default: c(1, 1)). Used by "lm"
#'   and "lm_hs".
#' @param unshrunk Integer column indices into \code{X} naming the coefficients
#'   that keep a normal prior instead of the horseshoe. In a causal model those
#'   are the treatment main effect and every effect modification term, because
#'   shrinking them attenuates the treatment effect. Used by the three
#'   horseshoe models.
#' @param beta0_prior Prior precision of the intercept under the horseshoe
#'   models (default: 1/100).
#' @param unshrunk_prior Prior precision of the coefficients listed in
#'   \code{unshrunk} (default: 1/100).
#' @param tau_prior Scale of the half-Cauchy prior on the global shrinkage
#'   parameter, or NULL (default) to use the value implied by \code{p0}
#'   following Piironen and Vehtari (2017), exactly as the Gibbs samplers do.
#' @param p0 A positive integer, a guess at how many coefficients are non-zero,
#'   used only when \code{tau_prior} is NULL (default: 5).
#' @param seed A single number passed to CmdStan, or NULL (default) to draw one
#'   from the caller's random number stream so that the fit is reproducible
#'   under a single \code{\link{set.seed}}.
#' @param adapt_delta Target average acceptance probability, or NULL (default)
#'   for 0.8 under the normal priors and 0.99 under the horseshoe, where the
#'   posterior geometry is harder and a larger value avoids most divergences.
#' @param max_treedepth Maximum depth of the NUTS trajectory tree (default: 10).
#' @param parallel_chains Number of chains to run at the same time (default: 1).
#'   Set it to \code{chains} on a machine with the cores to spare.
#' @param cache_dir Directory in which to keep the compiled model when the
#'   package was installed without compiled Stan executables. Defaults to
#'   \code{getOption("DRBayes.stan_cache")} and, failing that, to a
#'   subdirectory of the session's temporary directory, so nothing is written
#'   outside the session unless you ask for it. Point it at a permanent
#'   directory to compile once instead of once per session.
#' @param quiet Logical, whether to suppress the compiler and sampler progress
#'   messages (default: TRUE).
#' @param ... Further arguments passed to the \code{$sample()} method of the
#'   CmdStan model, for instance \code{init} or \code{refresh}.
#'
#' @return An array of posterior draws of the regression coefficients with
#'   dimensions \code{mc} by \code{chains} by \code{ncol(X) + 1}, indexed by
#'   iteration, then chain, then parameter, with "(Intercept)" and the column
#'   names of X on the third margin. This is the same layout, and the same
#'   dimnames, that the Gibbs samplers return.
#'
#'   Two attributes carry the information a Gibbs sampler cannot report:
#'   \describe{
#'     \item{nuts_diagnostics}{A data frame with one row per chain and columns
#'       \code{chain}, \code{divergent} (post-warm-up divergent transitions),
#'       \code{max_treedepth} (iterations that hit \code{max_treedepth}),
#'       \code{ebfmi} (the energy Bayesian fraction of missing information),
#'       \code{step_size} and \code{accept_stat}.}
#'     \item{stan_info}{A list recording the model name, the number of warm-up
#'       iterations, the seed, the CmdStan version and the total sampling time
#'       in seconds.}
#'   }
#'   They are attributes rather than list elements so that the return value is
#'   still an array, and so passing it to \code{\link{drbayes_pc}} needs no
#'   unwrapping. Attributes on an object that gets subset are dropped, so read
#'   them off the value this function returns rather than off a slice of it.
#'
#' @section Why Stan is not the default:
#' The CRAN build machines do not have CmdStan, so the macOS and Windows
#' binaries CRAN distributes cannot carry compiled Stan models, and users of
#' those binaries would have to install from source to get them. Most R users
#' install binaries, so making Stan the default would leave the majority with a
#' default path that does not work.
#'
#' Speed is not the reason to switch either. About 98 percent of a
#' \code{\link{drbayes_pc}} call is the coupling sweep rather than the model
#' fit, and on the data generating mechanism of the paper one Stan fit takes
#' roughly 1.9 seconds against 1.0 for the corresponding Gibbs sampler.
#'
#' The reason to offer Stan is capability. NUTS reports divergent transitions,
#' E-BFMI and treedepth saturation, so a fit that has gone wrong says so
#' instead of quietly returning biased draws. Its non-centred parameterisation
#' of the horseshoe also handles the funnel shaped posterior that a Gibbs
#' sampler explores only slowly.
#'
#' @section Divergences under the horseshoe:
#' The horseshoe models often report a handful of divergent transitions even at
#' \code{adapt_delta = 0.99}. The half-Cauchy tails of the local scales leave a
#' narrow neck in the posterior that the non-centred parameterisation reduces
#' but does not remove. The corresponding Gibbs samplers have the same trouble
#' with the same region and simply cannot report it. Raise \code{adapt_delta}
#' further if the warning persists, and compare the posterior means against
#' \code{\link{bayes_lm_hs}} before trusting either.
#'
#' @section Requirements:
#' CmdStan and the \code{cmdstanr} package must be installed. If they are not,
#' this function stops and names the two steps. It never falls back to a
#' different sampler: a user who asked for Stan must not be handed something
#' else without being told.
#'
#' If the installed package does not carry compiled Stan executables, which is
#' the case for a binary installation, the model is compiled on first use. That
#' takes about half a minute per model and is then cached, for the session by
#' default or permanently if you set \code{cache_dir}.
#'
#' @examples
#' # CmdStan is not installed with the package, so the example runs only where
#' # it is present. Everything else in DRBayes works without it. It is wrapped
#' # in donttest because loading a compiled Stan model costs several seconds
#' # before any sampling starts, which is longer than a routine example should
#' # take; the code itself is valid and does run.
#' \donttest{
#' if (requireNamespace("instantiate", quietly = TRUE) &&
#'     instantiate::stan_cmdstan_exists()) {
#'   set.seed(1)
#'   n <- 120
#'   X <- cbind(A = rbinom(n, 1, 0.5), X1 = rnorm(n))
#'   Y <- 1 + 2 * X[, "A"] - 0.5 * X[, "X1"] + rnorm(n)
#'
#'   draws <- bayes_stan("lm", Y, X, mc = 250, warmup = 250, chains = 2)
#'   print(dim(draws))
#'   print(attr(draws, "nuts_diagnostics"))
#'
#'   # The same array a Gibbs sampler would have produced, so drbayes_pc takes
#'   # it as it stands.
#'   print(apply(draws, 3, mean))
#' }
#' }
#'
#' @references
#' Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B.,
#' Betancourt, M., Brubaker, M., Guo, J., Li, P., and Riddell, A. (2017).
#' Stan: A probabilistic programming language. Journal of Statistical Software,
#' 76(1).
#'
#' Piironen, J. and Vehtari, A. (2017). Sparsity information and regularization
#' in the horseshoe and other shrinkage priors. Electronic Journal of
#' Statistics, 11(2), 5018-5051.
#'
#' @seealso \code{\link{bayes_lm}} and the other Gibbs samplers, which fit the
#'   same models without needing CmdStan, and \code{\link{drbayes_pc}}, which
#'   consumes the draws.
#'
#' @export
bayes_stan <- function(model, Y, X, mc = 5000, chains = 4L, warmup = 1000,
                       theta_prior = NULL, sigma_prior = c(1, 1),
                       unshrunk = integer(0), beta0_prior = NULL,
                       unshrunk_prior = NULL, tau_prior = NULL, p0 = 5,
                       seed = NULL, adapt_delta = NULL, max_treedepth = 10L,
                       parallel_chains = 1L, cache_dir = NULL, quiet = TRUE,
                       ...) {

  model  <- validate_stan_model_name(model)
  binary <- model %in% c("logit", "probit", "logit_hs", "probit_hs")
  X      <- validate_sampler_inputs(Y, X, mc, binary = binary)
  chains <- validate_chains(chains)
  warmup <- validate_stan_count(warmup, "warmup")
  parallel_chains <- validate_stan_count(parallel_chains, "parallel_chains")

  stan_data <- build_stan_data(model, Y, X, theta_prior = theta_prior,
                               sigma_prior = sigma_prior, unshrunk = unshrunk,
                               beta0_prior = beta0_prior,
                               unshrunk_prior = unshrunk_prior,
                               tau_prior = tau_prior, p0 = p0)

  if (is.null(adapt_delta)) {
    # The horseshoe posterior is funnel shaped even under the non-centred
    # parameterisation, and 0.8 leaves hundreds of divergences in its neck.
    # Measured on a 250 by 11 design, 0.8 gave 211 divergences and 0.99 gave
    # 33, for about 60 percent more sampling time.
    adapt_delta <- if (endsWith(model, "_hs")) 0.99 else 0.8
  }
  if (!is.numeric(adapt_delta) || length(adapt_delta) != 1L ||
      !is.finite(adapt_delta) || adapt_delta <= 0 || adapt_delta >= 1) {
    stop("adapt_delta must be a single number strictly between 0 and 1, the ",
         "target acceptance probability, or NULL for the default")
  }
  if (is.null(seed)) {
    # Taken from the caller's stream so the whole fit is reproducible under one
    # set.seed(), the same contract the Gibbs samplers' starting values follow.
    seed <- sample.int(.Machine$integer.max, 1L)
  }

  compiled <- stan_model_object(model, cache_dir = cache_dir, quiet = quiet)

  # Each fit gets its own directory for the CmdStan output. Left to itself
  # cmdstanr names those files from the model, the clock to the minute and a
  # suffix drawn from R's random number stream, so two fits of the same model
  # started from the same set.seed() in the same minute collide and the second
  # one fails to read its own draws. Everything needed is read below, so the
  # directory goes away with the call rather than filling up the session.
  output_dir <- tempfile("DRBayes-stan-")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  fit <- compiled$sample(data = stan_data, seed = seed, chains = chains,
                         parallel_chains = parallel_chains,
                         output_dir = output_dir,
                         iter_warmup = warmup, iter_sampling = mc,
                         adapt_delta = adapt_delta,
                         max_treedepth = max_treedepth,
                         refresh = if (quiet) 0 else NULL,
                         show_messages = !quiet, show_exceptions = !quiet,
                         # Under quiet the diagnostics are still computed and
                         # still warned about, just once, by this function
                         # rather than twice.
                         diagnostics = if (quiet) NULL else
                           c("divergences", "treedepth", "ebfmi"),
                         ...)

  draws <- stan_theta_array(fit, mc = mc, chains = chains, X = X)

  diagnostics <- stan_nuts_diagnostics(fit, chains)
  warn_about_nuts_diagnostics(diagnostics, adapt_delta, max_treedepth)

  attr(draws, "nuts_diagnostics") <- diagnostics
  attr(draws, "stan_info") <- list(
    model = model, warmup = warmup, seed = seed,
    cmdstan_version = as.character(instantiate::stan_cmdstan_version()),
    time_seconds = unname(fit$time()$total))
  draws
}


#' Names of the Stan models shipped with the package
#'
#' Kept in one place because the file names in inst/stan, the model argument of
#' [bayes_stan()] and the data blocks all have to agree.
#'
#' @keywords internal
#' @noRd
stan_model_names <- function() {
  c("lm", "logit", "probit", "lm_hs", "logit_hs", "probit_hs")
}


#' Validate the requested model and turn it into a Stan file name
#'
#' @return The model name, without the "bayes_" prefix the files carry.
#'
#' @keywords internal
#' @noRd
validate_stan_model_name <- function(model) {
  if (!is.character(model) || length(model) != 1L || is.na(model)) {
    stop("model must be a single character string naming one of the models: ",
         toString(stan_model_names()))
  }
  # Accepting the function names too, since "bayes_lm" is what a user reading
  # the rest of the package has in mind.
  model <- sub("^bayes_", "", model)
  if (!(model %in% stan_model_names())) {
    stop("model must be one of ", toString(stan_model_names()), ", but is \"",
         model, "\"")
  }
  model
}


#' Validate a positive integer argument of the Stan backend
#'
#' @keywords internal
#' @noRd
validate_stan_count <- function(value, arg) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value != round(value) || value < 1L) {
    stop(arg, " must be a single positive integer")
  }
  as.integer(value)
}


#' Data list for one of the Stan models
#'
#' Mirrors the prior handling of the corresponding Gibbs sampler, reusing the
#' same validators, so that a user who moves from one backend to the other gets
#' the same posterior and the same error messages.
#'
#' @keywords internal
#' @noRd
build_stan_data <- function(model, Y, X, theta_prior, sigma_prior, unshrunk,
                            beta0_prior, unshrunk_prior, tau_prior, p0) {
  Z  <- cbind(1, X)
  nn <- nrow(X)
  pp <- ncol(X)

  stan_data <- list(N = nn, P = pp + 1L, Z = unname(Z))
  stan_data$Y <- if (model %in% c("lm", "lm_hs")) {
    as.numeric(Y)
  } else {
    as.integer(Y)
  }

  if (endsWith(model, "_hs")) {
    unshrunk <- validate_unshrunk(unshrunk, pp)
    shrunk   <- setdiff(seq_len(pp), unshrunk)
    if (length(shrunk) == 0L) {
      stop("Every column of X is listed in unshrunk, so there is nothing for ",
           "the horseshoe prior to shrink. Use model = \"",
           sub("_hs$", "", model), "\" instead.")
    }
    stan_data$K <- length(unshrunk)
    stan_data$M <- length(shrunk)
    # as.array() keeps a single index from being written to JSON as a scalar,
    # which Stan would then refuse as the wrong shape.
    stan_data$unshrunk_idx <- as.array(as.integer(unshrunk) + 1L)
    stan_data$shrunk_idx   <- as.array(as.integer(shrunk) + 1L)

    # Only the Gaussian outcome has a free scale, so only there do the fixed
    # variances need a default read off the data. The logit and probit latent
    # scales are fixed, which is what keeps 1/100 weakly informative for them.
    scaled <- if (identical(model, "lm_hs")) {
      default_coef_precision(Y, X)
    } else {
      rep(1 / 100, pp + 1L)
    }
    stan_data$beta0_prior <- validate_stan_precision(
      beta0_prior %||% scaled[1L], "beta0_prior")
    stan_data$unshrunk_prior <- as.array(validate_hs_precision(
      unshrunk_prior, scaled[unshrunk + 1L], "unshrunk_prior"))
    # The error scale each Gibbs sampler assumes when it places tau, so the two
    # backends start the horseshoe from the same global scale.
    sigma <- switch(model, lm_hs = NULL, logit_hs = 2, probit_hs = 1)
    stan_data$tau_prior <- resolve_tau_prior(tau_prior, p = length(shrunk),
                                             n = nn, p0 = p0, sigma = sigma)
  } else {
    # Matches the Gibbs samplers: the Gaussian outcome takes its default
    # precision from the scale of the data, the binary ones keep 1/100.
    default <- if (identical(model, "lm")) {
      default_coef_precision(Y, X)
    } else {
      1 / 100
    }
    stan_data$prior_precision <- unname(
      validate_prior_precision(theta_prior, pp + 1L, default = default))
  }

  if (model %in% c("lm", "lm_hs")) {
    if (length(sigma_prior) != 2L || !is.numeric(sigma_prior) ||
        any(sigma_prior <= 0)) {
      stop("sigma_prior must be two positive numbers, the shape and scale of ",
           "the inverse gamma prior on the error variance")
    }
    stan_data$sigma_a <- sigma_prior[1]
    stan_data$sigma_b <- sigma_prior[2]
  }

  stan_data
}


#' Validate a scalar prior precision
#'
#' @keywords internal
#' @noRd
validate_stan_precision <- function(value, arg) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value <= 0) {
    stop(arg, " is a prior PRECISION, so it must be a single positive number. ",
         "Pass 1/variance.")
  }
  as.numeric(value)
}


#' Compiled models already built in this session
#'
#' Compilation is checked, not repeated, when a model is asked for twice, and
#' this environment saves even the check.
#'
#' @keywords internal
#' @noRd
stan_model_env <- new.env(parent = emptyenv())


#' Whether a working CmdStan installation can be reached
#'
#' Wrapped in its own function so that the failure path can be tested on a
#' machine that does have CmdStan.
#'
#' @keywords internal
#' @noRd
cmdstan_is_available <- function() {
  requireNamespace("instantiate", quietly = TRUE) &&
    isTRUE(instantiate::stan_cmdstan_exists())
}


#' Stop unless CmdStan and cmdstanr are both available
#'
#' @keywords internal
#' @noRd
assert_cmdstan <- function() {
  if (cmdstan_is_available()) {
    return(invisible(TRUE))
  }
  if (!requireNamespace("instantiate", quietly = TRUE)) {
    stop("The Stan backend needs the instantiate package. Install it with ",
         "install.packages(\"instantiate\").", call. = FALSE)
  }
  stop("The Stan backend needs CmdStan, which was not found. Two steps ",
       "install it:\n",
       "  1. install.packages(\"cmdstanr\", repos = ",
       "c(\"https://stan-dev.r-universe.dev\", getOption(\"repos\")))\n",
       "  2. cmdstanr::install_cmdstan()\n",
       "Nothing is substituted for Stan here: use bayes_lm(), bayes_logit(), ",
       "bayes_probit() or their _hs versions if you want the Gibbs samplers ",
       "instead.", call. = FALSE)
}


#' Directory holding the compiled model when the installation has none
#'
#' Defaults inside the session's temporary directory, because a package must
#' not write outside it without being asked. Setting the option or the argument
#' is that request, and buys compilation once instead of once per session.
#'
#' @keywords internal
#' @noRd
resolve_stan_cache_dir <- function(cache_dir) {
  if (is.null(cache_dir)) {
    cache_dir <- getOption("DRBayes.stan_cache",
                           file.path(tempdir(), "DRBayes-stan"))
  }
  if (!is.character(cache_dir) || length(cache_dir) != 1L ||
      is.na(cache_dir) || !nzchar(cache_dir)) {
    stop("cache_dir must be a single directory path, or NULL for the default")
  }
  if (!dir.exists(cache_dir) &&
      !dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create the Stan cache directory \"", cache_dir,
         "\". Pass cache_dir, or set options(DRBayes.stan_cache = ...), to a ",
         "directory that can be written to.")
  }
  cache_dir
}


#' The cmdstanr model object for one of the shipped Stan files
#'
#' Two situations have to be told apart. When the package was configured with
#' instantiate and installed from source on a machine that had CmdStan, the
#' executables are already in bin/stan and are used as they are. Otherwise, and
#' that is what a binary installation looks like, the .stan file shipped in
#' inst/stan is copied out and compiled on demand, because compiling into the
#' library would write to a directory that may well be read only.
#'
#' @keywords internal
#' @noRd
stan_model_object <- function(model, cache_dir = NULL, quiet = TRUE) {
  assert_cmdstan()
  name <- paste0("bayes_", model)

  prebuilt <- system.file("bin", "stan", paste0(name, ".stan"),
                          package = "DRBayes")
  exe_name <- if (identical(.Platform$OS.type, "windows")) {
    paste0(name, ".exe")
  } else {
    name
  }
  if (nzchar(prebuilt) &&
      file.exists(file.path(dirname(prebuilt), exe_name))) {
    return(instantiate::stan_package_model(name = name, package = "DRBayes",
                                           compile = FALSE))
  }

  cache_dir <- resolve_stan_cache_dir(cache_dir)
  key <- paste(name, cache_dir, sep = "\r")
  if (!is.null(stan_model_env[[key]])) {
    return(stan_model_env[[key]])
  }

  # bin/stan carries the sources as well as the executables, so it is still the
  # place to read from when the install ran on a machine without CmdStan and
  # the compilation step was skipped.
  shipped <- if (nzchar(prebuilt)) {
    prebuilt
  } else {
    system.file("stan", paste0(name, ".stan"), package = "DRBayes")
  }
  if (!nzchar(shipped)) {
    stop("The Stan model file for \"", model, "\" is missing from the ",
         "installed DRBayes package. Reinstall the package.", call. = FALSE)
  }
  target <- file.path(cache_dir, paste0(name, ".stan"))
  # copy.date keeps the executable newer than its source, which is how cmdstanr
  # decides that a cached model needs no recompiling.
  if (!file.exists(target) ||
      !identical(readLines(target, warn = FALSE),
                 readLines(shipped, warn = FALSE))) {
    file.copy(shipped, target, overwrite = TRUE, copy.date = TRUE)
  }

  # instantiate knows where CmdStan is even when it was pinned at install time,
  # which cmdstanr on its own does not.
  path_old <- if (nzchar(cmdstanr::cmdstan_path())) {
    cmdstanr::cmdstan_path()
  } else {
    NULL
  }
  path_new <- instantiate::stan_cmdstan_path()
  if (nzchar(path_new) && !identical(path_new, path_old)) {
    suppressMessages(cmdstanr::set_cmdstan_path(path_new))
    if (!is.null(path_old)) {
      on.exit(suppressMessages(cmdstanr::set_cmdstan_path(path_old)),
              add = TRUE)
    }
  }

  compiled <- cmdstanr::cmdstan_model(stan_file = target, quiet = quiet)
  stan_model_env[[key]] <- compiled
  compiled
}


#' Coefficient draws in the layout the Gibbs samplers return
#'
#' @keywords internal
#' @noRd
stan_theta_array <- function(fit, mc, chains, X) {
  drawn <- fit$draws(variables = "theta")
  out   <- draw_storage(mc, chains, X)

  if (!identical(dim(drawn), dim(out))) {
    stop("Stan returned draws of dimension ", toString(dim(drawn)),
         " where ", toString(dim(out)), " was expected. This is a bug in the ",
         "DRBayes Stan backend; please report it.", call. = FALSE)
  }
  expected <- paste0("theta[", seq_len(dim(out)[3L]), "]")
  if (!identical(unname(dimnames(drawn)[[3L]]), expected)) {
    stop("Stan returned the coefficients in an unexpected order (",
         toString(dimnames(drawn)[[3L]]), "). This is a bug in the DRBayes ",
         "Stan backend; please report it.", call. = FALSE)
  }

  # Assigning through [] keeps the storage array's own dimnames, which are the
  # ones the Gibbs samplers use, and drops the posterior classes cmdstanr adds.
  out[] <- as.numeric(drawn)
  out
}


#' NUTS diagnostics, one row per chain
#'
#' These are what a Gibbs sampler cannot report, and the reason the Stan
#' backend exists. Divergent transitions and low E-BFMI both say the sampler
#' failed to explore part of the posterior, so the draws are biased in a way
#' R-hat and effective sample size will not necessarily reveal.
#'
#' @keywords internal
#' @noRd
stan_nuts_diagnostics <- function(fit, chains) {
  summary <- fit$diagnostic_summary(quiet = TRUE)
  sampler <- fit$sampler_diagnostics()

  chain_mean <- function(variable) {
    if (!(variable %in% dimnames(sampler)[[3L]])) {
      return(rep(NA_real_, chains))
    }
    apply(sampler[, , variable, drop = FALSE], 2L, mean)
  }

  data.frame(
    chain = seq_len(chains),
    divergent = as.integer(summary$num_divergent),
    max_treedepth = as.integer(summary$num_max_treedepth),
    ebfmi = as.numeric(summary$ebfmi),
    step_size = chain_mean("stepsize__"),
    accept_stat = chain_mean("accept_stat__"),
    row.names = NULL)
}


#' Warn when the NUTS diagnostics say the draws cannot be trusted
#'
#' A warning rather than an error: the draws are returned either way, because
#' the user may be diagnosing the problem and needs to look at them. What must
#' not happen is that the problem goes unmentioned.
#'
#' @keywords internal
#' @noRd
warn_about_nuts_diagnostics <- function(diagnostics, adapt_delta,
                                        max_treedepth) {
  problems <- character(0)

  divergent <- sum(diagnostics$divergent, na.rm = TRUE)
  if (divergent > 0L) {
    problems <- c(problems, paste0(
      divergent, " divergent transition(s) after warm-up. The sampler ",
      "missed part of the posterior, so the draws are biased. Raise ",
      "adapt_delta above ", adapt_delta, "."))
  }

  saturated <- sum(diagnostics$max_treedepth, na.rm = TRUE)
  if (saturated > 0L) {
    problems <- c(problems, paste0(
      saturated, " iteration(s) hit the maximum treedepth of ", max_treedepth,
      ". The draws are valid but inefficient. Raise max_treedepth."))
  }

  # Betancourt (2018) treats E-BFMI below 0.3 as a sign that the momentum
  # resampling cannot carry the chain across the energy distribution.
  low_ebfmi <- which(diagnostics$ebfmi < 0.3)
  if (length(low_ebfmi) > 0L) {
    problems <- c(problems, paste0(
      "E-BFMI below 0.3 in chain(s) ", toString(diagnostics$chain[low_ebfmi]),
      ". The chain is exploring the tails of the posterior slowly; ",
      "reparameterising the model usually helps more than running longer."))
  }

  if (length(problems) > 0L) {
    warning("Stan reported problems with this fit:\n  ",
            paste(problems, collapse = "\n  "),
            "\nSee attr(draws, \"nuts_diagnostics\") for the per chain counts.",
            call. = FALSE)
  }
  invisible(NULL)
}

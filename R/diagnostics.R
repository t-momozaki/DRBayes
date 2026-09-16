#' Rank normalise draws with Blom's transform
#'
#' Replaces each draw by its rank among the pooled draws of all chains and maps
#' that rank to a normal score,
#' \eqn{z = \Phi^{-1}((r - 3/8)/(S + 1/4))}, equation (4.1) of Vehtari et al.
#' (2021).
#'
#' Split-R-hat and the effective sample size are built from means and variances,
#' so they are only defined when the target has a finite mean and variance. The
#' posteriors coupled by this package can be heavy tailed, and a Cauchy-like
#' marginal makes the raw statistics numerically meaningless while still
#' returning a value close to one. Ranks exist whatever the tails do, and for a
#' near normal marginal the normal scores reproduce the untransformed
#' diagnostics, so nothing is lost by always transforming.
#'
#' Ties take the average rank. That keeps the number of distinct values of a
#' discrete quantity intact, which matters because the tail effective sample
#' size is computed from a zero-one indicator.
#'
#' @param x Draws as a matrix of iterations by chains, or a vector.
#'
#' @return An object of the same shape as `x` holding the normal scores.
#'
#' @keywords internal
#' @noRd
rank_normalise <- function(x) {
  r <- rank(as.vector(x), ties.method = "average", na.last = "keep")
  z <- stats::qnorm((r - 3 / 8) / (length(r) + 1 / 4))
  if (!is.null(dim(x))) {
    z <- array(z, dim = dim(x), dimnames = dimnames(x))
  }
  z
}


#' Split every chain in half
#'
#' Doubles the number of chains by cutting each one at its midpoint, as in
#' Section 3.1 of Vehtari et al. (2021). Without this a chain that drifts
#' steadily in one direction is invisible to R-hat: its running mean matches the
#' other chains well enough that the between-chain variance stays small, even
#' though no part of it is stationary. Comparing a chain against its own second
#' half exposes that drift.
#'
#' @param x Draws as a matrix of iterations by chains, or a vector.
#'
#' @return A matrix with half as many rows and twice as many columns. With an
#'   odd number of iterations the middle draw of each chain is dropped so that
#'   the two halves are the same length.
#'
#' @keywords internal
#' @noRd
split_chains <- function(x) {
  x <- as.matrix(x)
  n_iter <- nrow(x)
  if (n_iter < 2L) {
    return(x)
  }
  half <- n_iter %/% 2L
  cbind(x[seq_len(half), , drop = FALSE],
        x[(n_iter - half + 1L):n_iter, , drop = FALSE])
}


#' Are all draws numerically identical
#'
#' A parameter that never moves has zero within-chain variance, which would make
#' R-hat and the effective sample size divide by zero. Callers report NA for
#' such a parameter instead.
#'
#' @keywords internal
#' @noRd
draws_are_constant <- function(x) {
  abs(max(x) - min(x)) < .Machine$double.eps
}


#' Potential scale reduction factor of already prepared draws
#'
#' Equations (3.1) to (3.4) of Vehtari et al. (2021), applied as given. The
#' draws are used exactly as supplied, so the caller is responsible for having
#' split and transformed them; [rhat()] does both.
#'
#' @param x Draws as a matrix of iterations by chains.
#'
#' @return The scalar R-hat, or NA if the draws are constant, not all finite, or
#'   too short or too few to give the statistic any meaning.
#'
#' @keywords internal
#' @noRd
rhat_basic <- function(x) {
  x <- as.matrix(x)
  if (!all(is.finite(x)) || draws_are_constant(x)) {
    return(NA_real_)
  }
  n_iter  <- nrow(x)
  n_chain <- ncol(x)
  if (n_iter < 2L || n_chain < 2L) {
    return(NA_real_)
  }

  chain_mean <- colMeans(x)
  chain_var  <- apply(x, 2L, stats::var)
  b <- n_iter * stats::var(chain_mean)
  w <- mean(chain_var)
  var_plus <- (n_iter - 1) / n_iter * w + b / n_iter
  sqrt(var_plus / w)
}


#' Rank normalised, folded split-R-hat
#'
#' The convergence statistic to report: the larger of split-R-hat on the rank
#' normalised draws and split-R-hat on the rank normalised folded draws
#' \eqn{\zeta = |\theta - \mathrm{median}(\theta)|}, equation (4.2) of Vehtari
#' et al. (2021).
#'
#' Folding is what catches chains that agree on location but disagree on scale,
#' for instance when one chain is stuck near the middle of the distribution and
#' never visits the tails. Plain R-hat compares chain means, so it sees nothing
#' wrong there; the folded draws turn that difference in spread into a
#' difference in location, which R-hat can see.
#'
#' @param x Draws as a matrix of iterations by chains, or a vector.
#'
#' @return The scalar R-hat, or NA for degenerate draws.
#'
#' @keywords internal
#' @noRd
rhat <- function(x) {
  x <- as.matrix(x)
  # Rank after splitting rather than before: with an odd number of iterations
  # splitting discards the middle draw of each chain, and the normal scores must
  # be those of the draws that are actually used.
  rhat_bulk <- rhat_basic(rank_normalise(split_chains(x)))
  folded    <- abs(x - stats::median(x))
  rhat_fold <- rhat_basic(rank_normalise(split_chains(folded)))
  max(rhat_bulk, rhat_fold)
}


#' Autocovariance at every lag by fast Fourier transform
#'
#' Returns the biased estimate, that is the one with divisor N rather than
#' N - t, which Geyer (1992) recommends and Section 3.2 of Vehtari et al. (2021)
#' adopts: the unbiased version has so much variance at long lags that the
#' truncation rule below becomes unstable.
#'
#' The series is padded to at least twice its length before transforming,
#' because a discrete Fourier transform treats the input as periodic and without
#' the padding the tail of the series would wrap around and contaminate the
#' short lags. Normalising by the lag zero term rather than by the transform
#' length avoids relying on the scaling convention of `fft(inverse = TRUE)`.
#'
#' @param x A single chain, as a numeric vector.
#'
#' @return A numeric vector of length `length(x)` holding lags 0 to N - 1.
#'
#' @keywords internal
#' @noRd
autocovariance_fft <- function(x) {
  n <- length(x)
  x_var <- stats::var(x)
  if (x_var == 0) {
    return(rep.int(0, n))
  }
  padded_length <- 2L * stats::nextn(n)
  centred <- c(x - mean(x), rep.int(0, padded_length - n))
  acov <- Re(stats::fft(Mod(stats::fft(centred))^2, inverse = TRUE)[seq_len(n)])
  acov / acov[1L] * x_var * (n - 1) / n
}


#' Effective sample size of a set of draws
#'
#' Equations (3.10) to (3.13) of Vehtari et al. (2021). Chains are split before
#' anything else, so that a drifting chain is penalised here as well as in
#' R-hat.
#'
#' Autocorrelations are pooled across chains with equation (3.10),
#' \eqn{\hat\rho_t = 1 - (W - \frac{1}{M}\sum_m s_m^2 \hat\rho_{t,m}) /
#' \widehat{\mathrm{var}}^+}. Using the multi-chain variance in the denominator
#' is what makes the estimate conservative when the chains have not mixed: if
#' they sit in different places the pooled autocorrelation stays near
#' \eqn{1 - \hat{R}^{-2}} at every lag and the effective sample size collapses
#' towards the number of modes found, rather than reporting each badly mixed
#' chain as if it carried independent information.
#'
#' The sum over lags is truncated by Geyer's initial positive sequence rule,
#' applied to the pair sums \eqn{P_t = \hat\rho_{2t} + \hat\rho_{2t+1}}, which
#' are positive, monotone and convex for a reversible Markov chain up to
#' estimator noise. The monotone step removes what noise remains.
#'
#' @param x Draws as a matrix of iterations by chains, or a vector.
#' @param split Split the chains first. Set to FALSE only when the caller has
#'   already split them, as [ess_bulk()] has.
#'
#' @return The scalar effective sample size, or NA for degenerate draws.
#'
#' @keywords internal
#' @noRd
ess_basic <- function(x, split = TRUE) {
  x <- if (split) split_chains(x) else as.matrix(x)
  n_iter  <- nrow(x)
  n_chain <- ncol(x)
  if (n_iter < 3L || !all(is.finite(x)) || draws_are_constant(x)) {
    return(NA_real_)
  }

  acov <- vapply(seq_len(n_chain), function(m) autocovariance_fft(x[, m]),
                 numeric(n_iter))

  # The lag zero autocovariance has divisor N, so scaling by N / (N - 1) turns
  # every lag into s_m^2 * rho_{t,m}, the quantity equation (3.10) averages.
  scaled_acov <- rowMeans(acov) * n_iter / (n_iter - 1)
  w <- scaled_acov[1L]
  var_plus <- w * (n_iter - 1) / n_iter
  if (n_chain > 1L) {
    var_plus <- var_plus + stats::var(colMeans(x))
  }
  rho <- 1 - (w - scaled_acov) / var_plus

  # Geyer's rule needs the pair one past the truncation point, and the longest
  # lags are the noisiest, so the last few are never considered.
  last_pair <- max(0L, (n_iter - 4L) %/% 2L)
  pair_index <- 0:last_pair
  pairs <- rho[2L * pair_index + 1L] + rho[2L * pair_index + 2L]

  # Initial positive sequence: keep pairs up to the last one before the first
  # non-positive pair, always leaving one pair in hand for the even lag term.
  non_positive <- which(pairs[-1L] <= 0)
  k <- if (length(non_positive) > 0L) non_positive[1L] - 1L else last_pair - 1L
  k <- max(k, 0L)

  # Initial monotone sequence: a pair sum that rises above its predecessor can
  # only be noise, so pull it back down.
  kept <- cummin(pairs[seq_len(k + 1L)])

  tau <- -1 + 2 * sum(kept)
  # Section 3.2: average the sum truncated at the usual odd lag with the sum
  # truncated at the next even lag. The second of those is tau plus twice the
  # next even lag correlation, so their average is tau plus that correlation
  # once. For a strongly antithetic chain the individual odd lag correlations
  # stay large even where the pairs vanish, and averaging the two truncations
  # is the more stable of the estimates. A negative value here belongs to the
  # pair that failed the positivity test and is therefore noise; letting it
  # shrink tau would inflate the reported sample size.
  next_even <- 2L * (k + 1L) + 1L
  if (next_even <= n_iter && rho[next_even] > 0) {
    tau <- tau + rho[next_even]
  }

  s_total <- n_iter * n_chain
  # A very short antithetic chain can produce an absurdly small tau by chance.
  # Capping the reported size at S log10(S) keeps such an estimate finite.
  tau <- max(tau, 1 / log10(s_total))
  s_total / tau
}


#' Bulk effective sample size
#'
#' The effective sample size of the rank normalised draws, Section 4.1 of
#' Vehtari et al. (2021). This is the one to read alongside R-hat, since it is
#' well defined even when the marginal posterior has no finite mean or variance.
#' It is not the right quantity for the Monte Carlo error of the posterior mean,
#' which is what [mcse_mean()] uses the untransformed draws for.
#'
#' @param x Draws as a matrix of iterations by chains, or a vector.
#'
#' @return The scalar bulk effective sample size, or NA for degenerate draws.
#'
#' @keywords internal
#' @noRd
ess_bulk <- function(x) {
  ess_basic(rank_normalise(split_chains(x)), split = FALSE)
}


#' Tail effective sample size
#'
#' The smaller of the effective sample sizes of the 5\% and 95\% quantiles,
#' Section 4.3 of Vehtari et al. (2021).
#'
#' A quantile is not an expectation, but the cumulative probability at that
#' quantile is: the indicator \eqn{I(\theta \le \theta_\alpha)} turns the draws
#' into zeros and ones whose mean estimates it, and equation (4.4) shows the
#' quantile inherits the effective sample size of that mean. Reporting the tails
#' separately matters because sampling efficiency is rarely uniform over the
#' parameter space, and a posterior interval can be far less reliable than the
#' posterior mean computed from the same draws.
#'
#' @param x Draws as a matrix of iterations by chains, or a vector.
#'
#' @return The scalar tail effective sample size, or NA for degenerate draws.
#'
#' @keywords internal
#' @noRd
ess_tail <- function(x) {
  x <- as.matrix(x)
  if (!all(is.finite(x))) {
    return(NA_real_)
  }
  # The indicator is formed from the pooled draws before splitting, so both
  # halves of every chain are compared against the same quantile.
  ess_lower <- ess_basic(x <= stats::quantile(x, 0.05))
  ess_upper <- ess_basic(x <= stats::quantile(x, 0.95))
  min(ess_lower, ess_upper)
}


#' Monte Carlo standard error of the posterior mean
#'
#' The posterior standard deviation divided by the square root of the effective
#' sample size, equation (3.5) of Vehtari et al. (2021). Computed on the
#' untransformed draws, because it must be on the scale of the parameter to be
#' comparable with the posterior standard deviation reported next to it.
#'
#' @param x Draws as a matrix of iterations by chains, or a vector.
#'
#' @return The scalar Monte Carlo standard error, or NA for degenerate draws.
#'
#' @keywords internal
#' @noRd
mcse_mean <- function(x) {
  stats::sd(as.vector(x)) / sqrt(ess_basic(x))
}


#' Coerce posterior draws to an iterations by chains by parameters array
#'
#' Accepts the three-dimensional array the diagnostics are defined on, and also
#' the iterations by parameters matrix the package's own Gibbs samplers return,
#' which is the single chain case.
#'
#' @keywords internal
#' @noRd
as_draws_array <- function(draws, arg = "draws") {
  if (is.data.frame(draws)) {
    draws <- as.matrix(draws)
  }
  if (is.matrix(draws)) {
    parameters <- colnames(draws)
    draws <- array(draws, dim = c(nrow(draws), 1L, ncol(draws)))
    dimnames(draws) <- list(NULL, NULL, parameters)
  }
  if (!is.array(draws) || length(dim(draws)) != 3L) {
    stop(arg, " must be a three-dimensional array of iterations by chains by ",
         "parameters, or a matrix of iterations by parameters for a single ",
         "chain. Combine several chains with, for example, ",
         "array(c(chain1, chain2), dim = c(n_iter, 2, n_par)).")
  }
  if (!is.numeric(draws)) {
    stop(arg, " must be numeric")
  }
  if (dim(draws)[1L] < 4L) {
    stop(arg, " has only ", dim(draws)[1L], " iteration(s) per chain. The ",
         "split diagnostics need at least four, and are only informative with ",
         "far more than that.")
  }
  if (dim(draws)[3L] < 1L) {
    stop(arg, " contains no parameters")
  }
  if (is.null(dimnames(draws)[[3L]])) {
    dimnames(draws)[[3L]] <- paste0("par", seq_len(dim(draws)[3L]))
  }
  draws
}


#' Convergence diagnostics for posterior draws
#'
#' Summarises a set of MCMC draws with the rank normalised split-R-hat, the bulk
#' and tail effective sample sizes and the Monte Carlo standard error of the
#' posterior mean, following Vehtari, Gelman, Simpson, Carpenter and Burkner
#' (2021).
#'
#' @details
#' Posterior coupling tilts particles drawn from two separate posteriors, one
#' for the outcome model and one for the propensity score model. The doubly
#' robust estimate is only as good as those two inputs: if either sampler has
#' not converged, the tilted posterior is wrong and the estimate built on it is
#' wrong with it, in a way that the coupling itself cannot reveal. Run these
#' diagnostics on the draws from each model before coupling them.
#'
#' The statistics reported are
#' \describe{
#'   \item{`rhat`}{The larger of split-R-hat on the rank normalised draws and
#'     split-R-hat on the rank normalised folded draws. Rank normalisation keeps
#'     the statistic meaningful for a heavy tailed posterior, where R-hat
#'     computed from second moments returns a value near one whatever the chains
#'     are doing. Folding catches chains that share a location but differ in
#'     scale, which plain R-hat cannot see.}
#'   \item{`ess_bulk`}{The effective sample size of the rank normalised draws,
#'     which describes how well the centre of the distribution has been
#'     explored.}
#'   \item{`ess_tail`}{The smaller of the effective sample sizes of the 5\% and
#'     95\% quantiles, which describes how reliable a posterior interval is. It
#'     is often much smaller than the bulk value.}
#'   \item{`mcse_mean`}{The standard deviation of the draws divided by the
#'     square root of their effective sample size, on the scale of the
#'     parameter.}
#' }
#'
#' Vehtari et al. (2021), Section 2, recommend running at least four chains and
#' using the draws only if `rhat` is below 1.01 and both effective sample sizes
#' exceed 400. The R-hat threshold is much tighter than the 1.1 of Gelman and
#' Rubin (1992) because R-hat is not in fact a potential scale reduction factor
#' and can dip below 1.1 well before the chains have converged. The threshold of
#' 400 is the point at which the variances and autocorrelations that R-hat and
#' the effective sample size are themselves built from become stable, which is
#' why an effective sample size below 400 makes a small R-hat uninformative
#' rather than reassuring.
#'
#' A parameter that never moves, or whose draws are not all finite, is reported
#' with NA diagnostics rather than a value obtained by dividing by zero.
#'
#' @param draws Posterior draws, either a three-dimensional numeric array of
#'   iterations by chains by parameters with the parameter names as the third
#'   dimension name, or a matrix of iterations by parameters for a single chain.
#'
#' @return A data frame with one row per parameter and the columns `parameter`,
#'   `mean`, `sd`, `rhat`, `ess_bulk`, `ess_tail` and `mcse_mean`.
#'
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C.
#' (2021). Rank-normalization, folding, and localization: an improved R-hat for
#' assessing convergence of MCMC. \emph{Bayesian Analysis} 16(2), 667-718.
#' \doi{10.1214/20-BA1221}
#'
#' @examples
#' set.seed(1)
#' draws <- array(stats::rnorm(500 * 4 * 2), dim = c(500, 4, 2),
#'                dimnames = list(NULL, NULL, c("alpha", "beta")))
#' convergence_diagnostics(draws)
#'
#' @export
convergence_diagnostics <- function(draws) {
  draws <- as_draws_array(draws)
  parameters <- dimnames(draws)[[3L]]
  n_iter  <- dim(draws)[1L]
  n_chain <- dim(draws)[2L]

  summarise_one <- function(j) {
    x <- matrix(draws[, , j], nrow = n_iter, ncol = n_chain)
    c(mean      = mean(x),
      sd        = stats::sd(as.vector(x)),
      rhat      = rhat(x),
      ess_bulk  = ess_bulk(x),
      ess_tail  = ess_tail(x),
      mcse_mean = mcse_mean(x))
  }

  stats_matrix <- vapply(seq_along(parameters), summarise_one, numeric(6L))

  out <- data.frame(parameter = parameters,
                    t(stats_matrix),
                    row.names = NULL,
                    stringsAsFactors = FALSE)
  out
}


#' Parameters that fail the recommended convergence thresholds
#'
#' Reads a summary from [convergence_diagnostics()] and reports which parameters
#' fall short of the thresholds of Vehtari et al. (2021), so that a caller can
#' assemble a warning from the result.
#'
#' @details
#' Section 3.2 of Vehtari et al. (2021) warns about a trap that a naive check
#' walks straight into. R-hat is a ratio of variance estimates, and those
#' estimates are only stable once there are enough effectively independent
#' draws to compute them from; the authors put that at an average of 50 per
#' split chain, or 400 in total for the recommended four chains. Below that,
#' R-hat is noise, and a value under 1.01 is not evidence of anything. The
#' `rhat_unreliable` flag marks this case: when it is TRUE, the R-hat column
#' must not be read as reassurance, whatever it says.
#'
#' A parameter with NA diagnostics counts as a failure. Constant or non-finite
#' draws are exactly the situation in which convergence cannot be assessed, and
#' silently passing them would defeat the purpose of the check.
#'
#' @param diagnostics A data frame as returned by [convergence_diagnostics()].
#' @param rhat_max Largest acceptable R-hat. Defaults to 1.01.
#' @param ess_min Smallest acceptable bulk and tail effective sample size.
#'   Defaults to 400.
#'
#' @return A list with elements `ok`, TRUE when nothing failed; `parameters`,
#'   the names of all failing parameters; `rhat`, those failing the R-hat
#'   threshold; `ess`, those failing either effective sample size threshold; and
#'   `rhat_unreliable`, TRUE when some parameter has an effective sample size
#'   below `ess_min`, so that its R-hat cannot be trusted in either direction.
#'
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C.
#' (2021). Rank-normalization, folding, and localization: an improved R-hat for
#' assessing convergence of MCMC. \emph{Bayesian Analysis} 16(2), 667-718.
#' \doi{10.1214/20-BA1221}
#'
#' @keywords internal
#' @noRd
flag_convergence_failures <- function(diagnostics, rhat_max = 1.01,
                                      ess_min = 400) {
  required <- c("parameter", "rhat", "ess_bulk", "ess_tail")
  missing_columns <- setdiff(required, names(diagnostics))
  if (!is.data.frame(diagnostics) || length(missing_columns) > 0L) {
    stop("diagnostics must be a data frame as returned by ",
         "convergence_diagnostics(); it is missing the column(s) ",
         toString(missing_columns))
  }
  if (!is.numeric(rhat_max) || length(rhat_max) != 1L || !is.finite(rhat_max) ||
      rhat_max <= 0) {
    stop("rhat_max must be a single positive number")
  }
  if (!is.numeric(ess_min) || length(ess_min) != 1L || !is.finite(ess_min) ||
      ess_min <= 0) {
    stop("ess_min must be a single positive number")
  }

  parameters <- as.character(diagnostics$parameter)
  # NA fails: it means the parameter is constant or holds a non-finite draw, and
  # neither can be certified as converged.
  fails_rhat <- is.na(diagnostics$rhat) | diagnostics$rhat > rhat_max
  low_ess    <- is.na(diagnostics$ess_bulk) | diagnostics$ess_bulk < ess_min |
    is.na(diagnostics$ess_tail) | diagnostics$ess_tail < ess_min

  list(
    ok              = !any(fails_rhat | low_ess),
    parameters      = parameters[fails_rhat | low_ess],
    rhat            = parameters[fails_rhat],
    ess             = parameters[low_ess],
    rhat_unreliable = any(low_ess)
  )
}

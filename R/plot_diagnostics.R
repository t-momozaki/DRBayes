#' Rank plots of posterior draws
#'
#' @description
#' Draws the rank histograms proposed by Vehtari et al. (2021, section 4.5).
#' The draws of one parameter are ranked over all chains pooled together, and
#' the ranks belonging to each chain are then shown as a separate histogram
#' panel on a shared rank axis. If every chain is targeting the same posterior,
#' the ranks within each chain are uniformly spread, so every panel should sit
#' flat along the reference line drawn at the uniform expectation. A chain that
#' has settled in a different location leans towards the low or the high ranks;
#' a chain with a different scale is heaped in the middle or at both ends.
#'
#' @details
#' Rank plots are proposed as a replacement for trace plots, not as a companion
#' to them: unlike trace plots, they "don't tend to squeeze to a fuzzy mess when
#' used with long chains" (Vehtari et al., 2021, section 4.5). Departure from
#' uniformity is judged against a fixed reference height, which does not become
#' harder to read as the chains get longer.
#'
#' Posterior coupling tilts draws taken from two separate posteriors, one for
#' the outcome model and one for the propensity score model. The doubly robust
#' estimate is only as trustworthy as the worse of the two, so run this on the
#' parameters of both models rather than on the treatment effect alone.
#'
#' Note that the object returned by the estimators in this package is a set of
#' posterior draws of a derived quantity, not a set of MCMC chains, and a plot
#' of it against the iteration index is not a trace plot. This function needs
#' genuine per-chain draws: run the sampler more than once, from dispersed
#' starting points, and stack the results into the array described under
#' `draws`.
#'
#' @param draws A numeric array of posterior draws with dimensions iterations by
#'   chains by parameters. A matrix of iterations by chains is accepted as the
#'   single parameter case. Chain and parameter names are taken from the second
#'   and third `dimnames` when present.
#' @param parameter The parameter to plot, given either as a name matched
#'   against the parameter `dimnames` or as a column index into the third
#'   dimension. Defaults to the first parameter.
#' @param bins Number of histogram bins, spread evenly over the rank axis.
#'   Defaults to 20.
#' @param plot Set to FALSE to compute the rank counts without opening or
#'   drawing on a graphics device. The return value is unchanged.
#' @param main Overall title placed above the panels. Defaults to the parameter
#'   name.
#' @param col,border Fill and border colours of the histogram bars.
#' @param ref_col Colour of the uniform reference line.
#'
#' @return Invisibly, a list with components
#'   \describe{
#'     \item{counts}{A `bins` by chains matrix of rank counts, with the chain
#'       names as column names.}
#'     \item{breaks}{The `bins + 1` bin boundaries on the rank axis.}
#'     \item{expected}{The count expected in every cell of `counts` when the
#'       chains agree, namely the number of iterations divided by `bins`.}
#'     \item{parameter}{The name of the parameter that was plotted.}
#'   }
#'   These are everything needed to redraw the plot, and to test the
#'   computation where no graphics device is available.
#'
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C.
#' (2021). Rank-normalization, folding, and localization: an improved R-hat for
#' assessing convergence of MCMC. Bayesian Analysis, 16(2), 667-718.
#' \doi{10.1214/20-BA1221}
#'
#' @examples
#' set.seed(1)
#' draws <- array(stats::rnorm(4000), dim = c(1000, 4, 1),
#'                dimnames = list(NULL, paste0("chain", 1:4), "ATE"))
#'
#' # Chains that agree: four flat panels.
#' rank_plot(draws, "ATE")
#'
#' # A chain stuck at a different location: its panel tilts, and the other
#' # three tilt the other way to compensate.
#' draws[, 4, 1] <- draws[, 4, 1] + 1
#' rank_plot(draws, "ATE")
#'
#' # The counts alone, with no drawing at all.
#' rank_plot(draws, "ATE", plot = FALSE)$counts
#'
#' @export
rank_plot <- function(draws, parameter = 1L, bins = 20L, plot = TRUE,
                      main = NULL, col = "grey70", border = "white",
                      ref_col = "red") {

  draws <- validate_rank_plot_draws(draws)
  index <- resolve_parameter_index(parameter, dimnames(draws)[[3L]],
                                   dim(draws)[3L])

  if (!is.numeric(bins) || length(bins) != 1L || !is.finite(bins) ||
      bins != round(bins) || bins < 2L) {
    stop("bins must be a single whole number of at least 2")
  }
  bins <- as.integer(bins)
  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("plot must be TRUE or FALSE")
  }

  theta <- draws[, , index, drop = FALSE]
  dim(theta) <- dim(draws)[1:2]
  n_iter  <- nrow(theta)
  n_chain <- ncol(theta)

  if (!all(is.finite(theta))) {
    stop("Parameter '", dimnames(draws)[[3L]][index], "' contains ",
         sum(!is.finite(theta)), " missing or infinite draw(s). Ranking them ",
         "would silently pile them at one end of the rank axis; drop or ",
         "investigate them first.")
  }
  if (bins > n_iter) {
    stop("bins is ", bins, " but each chain has only ", n_iter,
         " iterations, so most bins would be empty for every chain. Use at ",
         "most ", n_iter, " bins.")
  }

  # Ranking over the pooled draws is what makes the panels comparable: each bin
  # holds the same number of pooled draws by construction, so any departure
  # from a flat panel is a statement about one chain relative to the others.
  # Average ranks for ties keep discrete or rounded draws from being ordered by
  # their position in the array.
  ranks <- rank(as.vector(theta), ties.method = "average")

  n_total <- n_iter * n_chain
  breaks  <- seq(0, n_total, length.out = bins + 1L)
  bin_of  <- .bincode(ranks, breaks, right = TRUE, include.lowest = TRUE)
  dim(bin_of) <- c(n_iter, n_chain)

  chain_names <- dimnames(draws)[[2L]]
  if (is.null(chain_names)) {
    chain_names <- paste0("chain ", seq_len(n_chain))
  }

  counts <- vapply(seq_len(n_chain),
                   function(m) tabulate(bin_of[, m], nbins = bins),
                   integer(bins))
  dim(counts) <- c(bins, n_chain)
  dimnames(counts) <- list(NULL, chain_names)

  expected  <- n_iter / bins
  parameter_name <- dimnames(draws)[[3L]][index]

  result <- list(counts    = counts,
                 breaks    = breaks,
                 expected  = expected,
                 parameter = parameter_name)

  if (plot) {
    draw_rank_panels(counts, breaks, expected, chain_names,
                     main = if (is.null(main)) parameter_name else main,
                     col = col, border = border, ref_col = ref_col)
  }

  invisible(result)
}


#' Coerce and check a draws array for the rank plot
#'
#' A matrix is read as iterations by chains for a single unnamed parameter,
#' which is the shape people reach for when they have only one quantity of
#' interest. Everything downstream can then assume three dimensions.
#'
#' @return The draws as a three dimensional array with a parameter name in
#'   `dimnames[[3]]`.
#'
#' @keywords internal
#' @noRd
validate_rank_plot_draws <- function(draws) {
  if (!is.numeric(draws)) {
    stop("draws must be a numeric array of iterations by chains by parameters")
  }
  if (is.matrix(draws)) {
    nms   <- dimnames(draws)
    draws <- array(draws, dim = c(dim(draws), 1L),
                   dimnames = list(nms[[1L]], nms[[2L]], "parameter"))
  }
  if (length(dim(draws)) != 3L) {
    stop("draws must have 3 dimensions, iterations by chains by parameters, ",
         "but has ", length(dim(draws)), ". A single chain still needs a ",
         "chain dimension: use array(x, dim = c(length(x), 1, 1)).")
  }
  if (dim(draws)[1L] < 2L) {
    stop("draws must have at least 2 iterations per chain, but has ",
         dim(draws)[1L])
  }
  if (dim(draws)[2L] < 1L || dim(draws)[3L] < 1L) {
    stop("draws must have at least one chain and one parameter, but has ",
         dim(draws)[2L], " chain(s) and ", dim(draws)[3L], " parameter(s)")
  }
  if (is.null(dimnames(draws)[[3L]])) {
    dimnames(draws)[[3L]] <- paste0("parameter ", seq_len(dim(draws)[3L]))
  }
  draws
}


#' Resolve a parameter given as a name or an index
#'
#' @param parameter A single name or index supplied by the user.
#' @param names Parameter names taken from the draws array.
#' @param n_param Number of parameters in the draws array.
#'
#' @return A single integer index.
#'
#' @keywords internal
#' @noRd
resolve_parameter_index <- function(parameter, names, n_param) {
  if (length(parameter) != 1L) {
    stop("parameter must be a single name or index, but has length ",
         length(parameter), ". Call the function once per parameter.")
  }
  if (is.character(parameter)) {
    index <- match(parameter, names)
    if (is.na(index)) {
      stop("There is no parameter called '", parameter, "' in draws. ",
           "Available: ", toString(names))
    }
    return(index)
  }
  if (!is.numeric(parameter) || !is.finite(parameter) ||
      parameter != round(parameter)) {
    stop("parameter must be a parameter name or a whole number index")
  }
  parameter <- as.integer(parameter)
  if (parameter < 1L || parameter > n_param) {
    stop("parameter index ", parameter, " is outside draws, which has ",
         n_param, " parameter(s)")
  }
  parameter
}


#' Draw one rank histogram panel per chain
#'
#' Kept apart from the computation so that `rank_plot(plot = FALSE)` touches no
#' graphics device at all. The bars are drawn with `rect()` rather than
#' `hist()` because the counts have already been formed from the pooled ranks
#' and must not be recomputed per panel: the bins have to stay identical across
#' panels for the panels to be comparable.
#'
#' @keywords internal
#' @noRd
draw_rank_panels <- function(counts, breaks, expected, chain_names, main,
                             col, border, ref_col) {
  n_chain <- ncol(counts)
  n_col   <- ceiling(sqrt(n_chain))
  n_row   <- ceiling(n_chain / n_col)

  # CRAN policy: leave the user's graphics state as we found it, whatever
  # happens in between.
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(n_row, n_col),
                mar   = c(3.1, 3.6, 2.1, 1.1),
                oma   = c(2.4, 0, 2.6, 0),
                mgp   = c(2.2, 0.7, 0))

  # A shared vertical axis, as well as a shared rank axis, so that a tall bar
  # in one panel cannot be mistaken for a tall bar in another.
  y_max <- max(max(counts), expected) * 1.1

  for (m in seq_len(n_chain)) {
    graphics::plot(NULL,
                   xlim = range(breaks), ylim = c(0, y_max),
                   xlab = "", ylab = "Count", main = chain_names[m],
                   xaxs = "i", yaxs = "i", bty = "n")
    graphics::rect(breaks[-length(breaks)], 0, breaks[-1L], counts[, m],
                   col = col, border = border)
    graphics::abline(h = expected, col = ref_col, lty = 2)
    graphics::box(bty = "l")
  }

  graphics::mtext("Rank of draw within the pooled chains", side = 1,
                  outer = TRUE, line = 1)
  graphics::mtext(main, side = 3, outer = TRUE, line = 0.8, font = 2)
  invisible(NULL)
}

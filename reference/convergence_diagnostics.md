# Convergence diagnostics for posterior draws

Summarises a set of MCMC draws with the rank normalised split-R-hat, the
bulk and tail effective sample sizes and the Monte Carlo standard error
of the posterior mean, following Vehtari, Gelman, Simpson, Carpenter and
Burkner (2021).

## Usage

``` r
convergence_diagnostics(draws)
```

## Arguments

- draws:

  Posterior draws, either a three-dimensional numeric array of
  iterations by chains by parameters with the parameter names as the
  third dimension name, or a matrix of iterations by parameters for a
  single chain.

## Value

A data frame with one row per parameter and the columns `parameter`,
`mean`, `sd`, `rhat`, `ess_bulk`, `ess_tail` and `mcse_mean`.

## Details

Posterior coupling tilts particles drawn from two separate posteriors,
one for the outcome model and one for the propensity score model. The
doubly robust estimate is only as good as those two inputs: if either
sampler has not converged, the tilted posterior is wrong and the
estimate built on it is wrong with it, in a way that the coupling itself
cannot reveal. Run these diagnostics on the draws from each model before
coupling them.

The statistics reported are

- `rhat`:

  The larger of split-R-hat on the rank normalised draws and split-R-hat
  on the rank normalised folded draws. Rank normalisation keeps the
  statistic meaningful for a heavy tailed posterior, where R-hat
  computed from second moments returns a value near one whatever the
  chains are doing. Folding catches chains that share a location but
  differ in scale, which plain R-hat cannot see.

- `ess_bulk`:

  The effective sample size of the rank normalised draws, which
  describes how well the centre of the distribution has been explored.

- `ess_tail`:

  The smaller of the effective sample sizes of the 5\\ 95\\ is often
  much smaller than the bulk value.

- `mcse_mean`:

  The standard deviation of the draws divided by the square root of
  their effective sample size, on the scale of the parameter.

Vehtari et al. (2021), Section 2, recommend running at least four chains
and using the draws only if `rhat` is below 1.01 and both effective
sample sizes exceed 400. The R-hat threshold is much tighter than the
1.1 of Gelman and Rubin (1992) because R-hat is not in fact a potential
scale reduction factor and can dip below 1.1 well before the chains have
converged. The threshold of 400 is the point at which the variances and
autocorrelations that R-hat and the effective sample size are themselves
built from become stable, which is why an effective sample size below
400 makes a small R-hat uninformative rather than reassuring.

A parameter that never moves, or whose draws are not all finite, is
reported with NA diagnostics rather than a value obtained by dividing by
zero.

## References

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C.
(2021). Rank-normalization, folding, and localization: an improved R-hat
for assessing convergence of MCMC. *Bayesian Analysis* 16(2), 667-718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

## Examples

``` r
set.seed(1)
draws <- array(stats::rnorm(500 * 4 * 2), dim = c(500, 4, 2),
               dimnames = list(NULL, NULL, c("alpha", "beta")))
convergence_diagnostics(draws)
#>   parameter        mean       sd      rhat ess_bulk ess_tail  mcse_mean
#> 1     alpha -0.01395503 1.037195 0.9996235 2081.234 2014.468 0.02272433
#> 2      beta  0.01601564 1.034629 1.0005338 1876.717 2045.842 0.02387183
```

# Sensitivity analysis for unmeasured confounding

Turns the posterior draws of a
[`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
fit into the posterior they would have been had the outcome model
carried an unmeasured confounding bias \\\xi\\. The original draws are
reweighted by the importance sampling of Algorithm 3 of Orihara,
Momozaki and Sugasawa (2025), or of Algorithm A.1 in their Appendix E,
and the reweighted draws are then coupled again by Algorithm 2, which is
what keeps the answer doubly robust. No model is refitted.

## Usage

``` r
drbayes_sensitivity(
  object,
  xi,
  M = 100,
  method = c("common", "per-observation"),
  sigma = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `"DRBayes"` from
  [`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md),
  fitted with `control = drbayes_control(keep_particles = TRUE)`. The
  untilted draws it keeps are what the weights are applied to, and the
  design matrices and responses it keeps with them are what the weights
  are computed from, so a fit without them cannot be used.

- xi:

  The sensitivity distribution \\g()\\, as either a function of `n`
  returning `n` draws, such as the closure
  [`xi_triangular()`](https://t-momozaki.github.io/DRBayes/reference/xi_triangular.md)
  builds, or a numeric vector of draws to resample. A single number is a
  point mass, and `xi = 0` leaves the weights uniform.

- M:

  Integer. Number of draws from `xi` used for the Monte Carlo integral
  in the weight. Default 100. The integrand concentrates far more
  sharply than \\g()\\ does, so far fewer than `M` of the draws carry
  the integral; how many is reported as `xi_ess` and warned about when
  it is small. Raising `M` raises that count in proportion, and the
  running time of the weights with it, while the share of the draws that
  counts stays where it is: on the gaussian example of this package's
  tests, `M` of 100, 1000 and 5000 left a median of 5.8, 65 and 347
  effective draws, 5.8 to 6.9 percent of `M` throughout. What can be
  bought is therefore the usual square root: ten times the work for a
  third of the Monte Carlo error. A \\g()\\ with less mass far from zero
  costs nothing instead.

  Algorithm 3 reads as `M` fresh draws for each of the `S` posterior
  draws. They are taken once and reused across the posterior draws here,
  which makes the analysis exact for the `M` atom distribution that was
  actually drawn, rather than an unbiased estimate of the analysis under
  \\g()\\ itself, and lets the result report which distribution that
  was. What it costs is that those atoms are themselves random: at
  `M = 200` on the gaussian example, repeating the call moves the
  g-formula shift by 3 percent of itself where fresh draws move it by
  0.3 percent. In `ate` it is not detectable, because the coupling gives
  most of the shift back and most of that scatter with it, leaving it
  well under the scatter of the sweep itself.

- method:

  `"common"` (default) for one \\\xi\\ shared by all observations,
  Algorithm A.1, or `"per-observation"` for one \\\xi_i\\ per
  observation, Algorithm 3. See Details.

- sigma:

  Residual standard deviation of a gaussian outcome model, which the
  likelihood ratio needs and the coefficient draws do not carry. Left
  NULL, one value is drawn for each coefficient draw from the
  conditional posterior of \\\sigma^2\\ given those coefficients, an
  inverse gamma with shape \\n/2\\ and rate half the residual sum of
  squares at them, so that the posterior uncertainty about the residual
  scale is carried into the weights rather than assumed away. That
  conditional is the outcome likelihood under a flat prior on \\\log
  \sigma^2\\, and it is used whatever sampler produced the draws: a fit
  records which prior its sampler put on \\\sigma^2\\ but not the prior
  itself, so that prior is not reproduced, which at the default of
  [`bayes_lm()`](https://t-momozaki.github.io/DRBayes/reference/bayes_lm.md)
  is worth 0.2 percent of a posterior standard deviation of \\\sigma\\.
  A single number fixes the scale instead, which understates the spread
  of the reweighted posterior by the amount quoted in Details; a vector
  of one value per posterior draw, from a sampler that returns its
  \\\sigma^2\\ chain, is used as it stands and is the exact thing the
  default approximates. Ignored for a binary outcome.

- ...:

  Not used, and refused rather than ignored so that a mistyped argument
  is not swallowed silently.

## Value

An object of class `"DRBayes_sensitivity"`, a list containing:

- ate:

  Numeric vector of draws of the average treatment effect under the
  assumed unmeasured confounding, reweighted, resampled and coupled
  again. They carry equal weight.

- ate_gcomp:

  The same draws before the coupling: the g-formula average treatment
  effect under the assumed unmeasured confounding, which is where the
  whole of the shift is. See Details.

- ate_original:

  The `pc` draws of the fit that was reweighted, for comparison.

- weights:

  Numeric vector of self-normalised sensitivity weights on the untilted
  draws, one per draw, summing to one. The resampling has already
  applied them; they are kept because they are what the diagnostics
  below are computed from.

- log_weights:

  The weights of Algorithm 3 or A.1 on the log scale, before
  normalising. They run to large negative numbers, which is why they are
  never exponentiated on their own.

- ess:

  Effective sample size of the weights, \\1/\sum w^2\\. The number of
  the original draws the answer really rests on.

- smc:

  The diagnostics of the coupling of the reweighted draws, as
  [`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
  reports them in its own `smc`: the tilting parameter `lambda`, the
  posterior mean `B_mean` of the moment condition (3.4) and the `tol` it
  is judged against, the number of steps taken, and the particle
  diagnostics. `converged` is whether the sensitivity analysis is still
  doubly robust, which takes the reweighting into account as well as the
  coupling: `converged_coupling` is the coupling's own verdict, and it
  is held to `ess_sensitivity` reaching `ess_floor`, which is `ess_frac`
  of the draws. See "Reading the result".

- xi:

  The draws from \\g()\\ that were used: a vector of length `M` for
  `method = "common"`, or a treated observations by `M` matrix for
  `method = "per-observation"`.

- xi_ess:

  The effective number of those draws that carries the Monte Carlo
  integral in the weight, as the `median` and the `min` over the
  integrals computed. Usually far smaller than `M`, and equal to it only
  when every draw of \\\xi\\ counts equally, as under a point mass.

- n_saturated:

  How many (treated observation, posterior draw) pairs the shift took
  outside the range of the outcome, `saturated` `of` the pairs there
  are, and so zero for a gaussian outcome. See "The shift is on the
  conditional mean".

- method, M, n_draws:

  The settings the weights were computed with.

- sigma:

  The residual standard deviations the gaussian likelihood ratio used,
  one per draw unless a single number was supplied, or NA for a binary
  outcome, whose likelihood has no scale parameter.

- call:

  The matched call.

## Details

**What xi is.**

Without strong ignorability the contrast the data identify is not the
average treatment effect but \$\$E\[E\[Y \| A = 1, X\] - E\[Y \| A = 0,
X\]\] = \tau + \xi,\$\$ so \\\xi\\ is the part of the observed contrast
that unmeasured confounding contributes, and \\\xi = 0\\ is strong
ignorability (Daniels, Linero and Roy, 2023, and Section 6.2 of the
paper). It enters through the shifted outcome model \\m_A(X; \beta, \xi)
= m_A(X; \beta) + A\xi\\: only treated observations are shifted, so only
they contribute to the weights.

**The shift is on the conditional mean.**

\\m_a(X)\\ is \\E\[Y \| A = a, X\]\\, so the shift is added to the
fitted mean and not to the linear predictor. For a gaussian outcome
those are the same thing. For a binary one they are not, and it is the
mean that is meant: \\\xi\\ is a shift of the probability, which is what
puts it in the units of \\\tau\\ in the display above and on the risk
difference scale that Section 6.3 reports. Added to the log odds instead
it would move the average treatment effect by several times less. On the
binary example of this package's tests, whose treated observations have
a mean fitted risk of 0.17, \\E\[\xi\] = 0.013\\ moves the g-formula
effect by -0.013 as a shift of the mean and by -0.002 as a shift of the
log odds, a factor of 8.

A probability cannot leave \\\[0, 1\]\\, so where \\m_A(X; \beta) +
\xi\\ would, the shifted mean is held at the boundary and that
observation carries a bias smaller than \\\xi\\. How many fitted values
that happens to is warned about and reported as `n_saturated`. It is a
statement about the outcome and not a numerical convenience: a treated
group whose fitted risk is 0.1 has no room for a risk difference of -0.5
of unmeasured confounding, so the \\Tri(-0.5, 0)\\ of Section 6.3 is a
bias a rare binary outcome cannot carry, and this reports how far short
of it the analysis fell.

A positive \\\xi\\ means part of what was read as a treatment effect was
confounding, so the reweighted posterior of the average treatment effect
is smaller than the original one, by about \\E\[\xi\]\\ while the
reweighting is still well resolved: on the gaussian example of this
package's tests \\E\[\xi\] = 0.033\\ moves it by -0.033, and on the
binary one \\E\[\xi\] = 0.013\\ moves it by -0.013. Beyond that the
shift falls short, because importance sampling cannot move a posterior
further than the draws it started from reach. And the coupling that
follows the reweighting gives most of even that back, for the reason set
out two sections below.

**Which algorithm.**

`method = "common"` is Algorithm A.1: one \\\xi\\ is shared by every
observation, the sum over observations happens inside the integral, and
the whole uncertainty in \\g()\\ is carried into the reweighted
posterior, which therefore comes out wider than the original by about
the variance of \\\xi\\. This is the default, and it is what the single
sensitivity parameter \\\xi = \Delta\\ of Section 6.2 means.

`method = "per-observation"` is Algorithm 3: each observation has its
own \\\xi_i\\ drawn from \\g()\\, and each integral is done separately
before the product over observations. Independent biases average out
over the sample, so this shifts the posterior by \\E\[\xi\]\\ while
adding almost none of the spread of \\g()\\; use it when the bias really
does vary from one observation to the next.

**Reweight first, then couple.**

The weights of both algorithms are ratios of outcome model likelihoods,
so they belong to the untilted posterior \\p_n(\alpha, \beta \| D)\\ of
equation (3.3) and that is what they are applied to here. The reweighted
draws are then resampled to equal weight and coupled again, which is the
order Section 6.2 of the paper prescribes: the posterior samples under
unmeasured confounding are obtained first, and the doubly robust
posterior of the average treatment effect follows from them through
Algorithm 1 or 2.

Reweighting the coupled draws instead would leave the moment condition
(3.4) violated, because importance weights move the particle cloud and
nothing moves it back. On the gaussian example of this package's tests
that costs a factor of 14 to 55 on \\\|mean B_n\|\\ relative to the
tolerance the fit itself was held to, over sensitivity distributions
from \\Tri(0, 0.05)\\ to \\Tri(0, 0.5)\\; coupling the reweighted draws
instead brings every one of them back inside that tolerance. The
recoupling is what the `$smc` component of the result reports, and its
`converged` entry is the verdict on the analysis as a whole rather than
on the coupling alone; see "Reading the result". It is printed with the
result.

**What the coupling does to the shift.**

The moment condition (3.4) is written with the unshifted mean function
\\m_A(X; \beta)\\ and the propensity score of the same fit, and both of
those describe the data as observed, unmeasured confounding included.
Imposing it on the reweighted draws therefore pulls the treated fitted
values back to what the propensity score weighted data say, and most of
the shift with them. Measured on the gaussian example, average treatment
effect 2.052:

|  |  |  |  |  |
|----|----|----|----|----|
| **g()** | **\\E\[\xi\]\\** | **after reweighting** | **after coupling** | **effective draws** |
| Tri(0, 0.05) | 0.033 | -0.033 | -0.001 | 2355 of 3000 |
| Tri(0, 0.20) | 0.134 | -0.103 | -0.006 | 191 of 3000 |
| Tri(0, 0.50) | 0.334 | -0.120 | -0.019 | 50 of 3000 |

Each entry is the mean of five runs, because none of these is a fixed
number. The sweep moves the fourth column by about 0.002 from one run to
the next in the first two rows and by 0.02 in the third, which is the
column's own message: what survives the coupling is of the size of the
noise of the sweep that produced it.

Only the first row is an answer. The other two rest on far fewer
effective draws than the floor of
[`drbayes_control()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_control.md),
so both are reported as not doubly robust and have their credible
intervals withheld, whether or not their coupling happens to meet the
moment condition. It did in two of the five runs each, which is the coin
flip the next section is about. They are tabulated to show what the
coupling does to a shift, not as results.

The spread of \\g()\\ survives the coupling but the location shift
largely does not, which is why the result reports both rows: `ate` after
the coupling and `ate_gcomp` before it. This is not an artefact of this
implementation. Section 6.3 of the paper reports the same thing: with
\\E\[\xi\] = 1/3\\ and with \\E\[\xi\] = -1/3\\ the average treatment
effect moves from 0.010 to 0.002 and to -0.004, about 0.01 in both cases
and in the same direction in both, which is the size of a shift that has
been given back rather than one of \\\pm 1/3\\. Read `ate` as the doubly
robust estimand under the assumed bias and `ate_gcomp` as the g-formula
one, and expect the sensitivity of the first to \\\xi\\ to be the milder
of the two.

The recoupling is a fresh sequential Monte Carlo sweep, so `xi = 0`
returns the original analysis only up to the Monte Carlo error of that
sweep, exactly as two calls to
[`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
on the same data differ. It is run with the `control` of the original
fit and against the tolerance that fit was judged against, so that the
same bar is applied to both couplings, and it costs about what the
coupling inside
[`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
cost. What is not inherited is the number of draws the mean being
compared with that bar is taken over; the next section is about what
follows from that.

**Reading the result.**

Importance reweighting can only report what the original draws already
covered. The reweighted posterior of the average treatment effect is the
original one shifted by about \\-E\[\xi\]\\, so the two overlap only
while \\\|E\[\xi\]\|\\ is within a few posterior standard deviations of
the effect being estimated. That ratio, and not the size of \\\xi\\ on
its own, is what decides whether the answer means anything: on a fit
whose posterior standard deviation is 0.11, the paper's own \\Tri(0,
0.5)\\ has \\E\[\xi\] = 1/3\\, three posterior standard deviations away,
and leaves an effective sample size of under 2 percent of the draws.

Two things have to hold before the analysis is still doubly robust, and
`smc$converged` is both of them. The coupling of the reweighted draws
has to meet the moment condition (3.4), which on its own is
`smc$converged_coupling`. And the reweighting that came before it has to
have left enough of the original draws in play, which is `ess` against
the floor `ess_frac` of
[`drbayes_control()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_control.md),
the same floor the coupling holds its own weights to. The second is not
something the coupling can see: it resamples the reweighted draws and
its kernel moves every particle, so the cloud looks fully populated to
every diagnostic it computes even where 28 of the 3000 draws are behind
it. Nor is the moment condition evidence at that point. The tolerance
was set from the fit's own 3000 draws, while the mean being compared
with it is now an average over an effective few dozen and carries
several times the tolerance as its own Monte Carlo error, so it is met
about as often as not: six runs of the identical call at \\Tri(0, 0.3)\\
on the gaussian example met it four times.

The credible interval is reported only where the effective sample size
can support it. A 2.5 percent point needs ten effective draws beyond it
to be a quantile rather than the most extreme handful that survived the
reweighting, so it needs 400 effective draws, and a median needs 20.
Below that the quantile is `NA` and the mean and the standard deviation
are what there is to read.

When the weights degenerate there are three things to do, in order of
how much they help: run
[`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
with more draws, since the effective sample size grows in proportion to
them; assume a sensitivity distribution with less mass far from zero,
which is a statement about the science and not a convenience; or trace
the analysis over a sequence of smaller \\\xi\\ and read off where the
conclusion changes, which is what a sensitivity analysis is for and
needs no single distribution to be believed.

**The residual scale of a gaussian outcome.**

The gaussian log likelihood ratio is a difference of squared residuals
over \\2\sigma^2\\ summed over the treated observations, so \\\sigma\\
is not a detail of the weight but the thing that sets its size. A point
estimate of it conditions the answer on a number the posterior is
uncertain about: on the gaussian example of this package's tests,
holding \\\sigma\\ one posterior standard deviation of its own either
side of 1, which is \\\sigma / \sqrt{2n} = 0.035\\, moves the shift the
reweighting produces by 13 percent and the effective sample size of the
weights by a factor of 2. By default \\\sigma\\ is therefore not fixed:
see the `sigma` argument.

## References

Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
Causal Inference via Posterior Coupling. arXiv preprint
arXiv:2506.04868.

Daniels, M. J., Linero, A., & Roy, J. (2023). *Bayesian Nonparametrics
for Causal Inference and Missing Data*. New York: Chapman and Hall/CRC.

## See also

[`xi_triangular()`](https://t-momozaki.github.io/DRBayes/reference/xi_triangular.md)
for the triangular sensitivity distributions of Section 6.3,
[`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
for the fit being reweighted, and
[`drbayes_control()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_control.md)
for `keep_particles`.

## Examples

``` r
# generate_dataset() leaves the caller's random number stream alone, so
# seed the coupling itself to make the fit below reproducible.
set.seed(1)
data <- generate_dataset(nn = 120, pp = 0, seed = 1)
fit  <- drbayes_pc(YY ~ AA + WW1 + WW2 + WW3 + WW4,
                   AA ~ WW1 + WW2 + WW3 + WW4,
                   data = data,
                   outcome.model = bayes_lm, ps.model = bayes_logit,
                   mc = 900, bn = 300, chains = 2L, verbose = FALSE,
                   control = drbayes_control(keep_particles = TRUE,
                                             n_steps = 200))
#> Warning: The posterior draws do not meet the convergence criteria (R-hat < 1.01, ESS > 400). outcome model: min ESS 355.7 (WW1). An effective sample size below 400 also makes R-hat itself unreliable, so a small R-hat here is not evidence of convergence. Raise mc, or thin less. Inspect $diagnostics.

# A bias somewhere in (0, 0.1) and most likely at its largest, which is the
# shape Section 6.3 of the paper assumes. What a sensitivity parameter has
# to be read against is not the size of the effect, 110 here, but the
# posterior standard deviation of it, which is the sd column of the summary
# below and about 0.16. The mean of xi, 0.067, is a little under half of
# that: a large bias to entertain, and one the original draws can still
# represent, since two thirds of them stay effective.
sens <- drbayes_sensitivity(fit, xi = xi_triangular(0, 0.1), M = 200)
sens
#> Sensitivity analysis for unmeasured confounding
#> 
#> Call:
#> drbayes_sensitivity(object = fit, xi = xi_triangular(0, 0.1), 
#>     M = 200)
#> 
#> Sensitivity parameter: one xi shared by every observation (Algorithm A.1)
#>                        200 draws from g(): mean 0.0669, range 0.00557 to 0.0999
#>                        integral resolved by a median of 197 of them, 123 at worst
#> Residual scale sigma:  drawn per posterior draw: mean 1.08, sd 0.0727
#> Posterior draws:       600
#> Effective sample size: 480 of 600 draws (79.9%)
#> Coupled again at:      lambda = 0.1
#> Moment condition:      |mean B_n| = 0.00917 against a tolerance of 0.0114, satisfied
#> Still doubly robust:   yes
#> 
#> Average treatment effect
#>                mean    2.5%   97.5%
#> sensitivity 110.029 109.610 110.422
#> reweighted  110.037 109.632 110.439
#> original    110.041 109.587 110.462
#> 
#> The sensitivity row is the posterior under the assumed unmeasured confounding;
#> the reweighted row is the same draws before the moment condition was imposed on
#> them again, and holds the whole of the shift; the original row is the fit that
#> was reweighted. See the "Reweight first, then couple" section of
#> ?drbayes_sensitivity for why the first two differ.
summary(sens)
#> Sensitivity analysis for unmeasured confounding
#> 
#> Call:
#> drbayes_sensitivity(object = fit, xi = xi_triangular(0, 0.1), 
#>     M = 200)
#> 
#> Sensitivity parameter: one xi shared by every observation (Algorithm A.1)
#> Draws from g():        200
#>                        integral resolved by a median of 197 of them, 123 at worst
#> Residual scale sigma:  drawn per posterior draw: mean 1.08, sd 0.0727
#> Posterior draws:       600
#> Effective sample size: 480 of 600 draws (79.9%)
#> Coupled again at:      lambda = 0.1
#> Moment condition:      |mean B_n| = 0.00917 against a tolerance of 0.0114, satisfied
#> Still doubly robust:   yes
#> 
#> Sensitivity parameter xi
#>    mean     sd     min    max
#>  0.0669 0.0235 0.00557 0.0999
#> 
#> Posterior summary of the average treatment effect
#>                mean    sd    2.5%     50%   97.5%
#> sensitivity 110.029 0.218 109.610 110.021 110.422
#> reweighted  110.037 0.212 109.632 110.032 110.439
#> original    110.041 0.227 109.587 110.036 110.462
#> 
#> The sensitivity row is the posterior under the assumed unmeasured confounding;
#> the reweighted row is the same draws before the moment condition was imposed on
#> them again, and holds the whole of the shift; the original row is the fit that
#> was reweighted. See the "Reweight first, then couple" section of
#> ?drbayes_sensitivity for why the first two differ.

# Whether the answer is still doubly robust: the coupling met the moment
# condition, on a reweighting that had draws left to spare.
sens$smc$converged
#> [1] TRUE
```

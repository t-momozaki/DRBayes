# Tuning Parameters for Posterior Coupling

Collects the settings that control the Markov chain Monte Carlo
sampling, the sequential Monte Carlo sweep over the tilting parameter,
and the convergence check performed by
[`drbayes_pc`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md).
Pass the result as the `control` argument, in the same way
[`glm.control`](https://rdrr.io/r/stats/glm.control.html) is passed to
[`glm`](https://rdrr.io/r/stats/glm.html). Only the settings you name
have to be given, everything else keeps its default.

## Usage

``` r
drbayes_control(
  mc = 5000,
  bn = 1000,
  thin = 2,
  chains = 4L,
  init = NULL,
  lambda_max = 10,
  n_steps = 1000,
  newton_steps = 100L,
  smoothing = 0.99,
  tol = NULL,
  ridge = 1e-05,
  ess_frac = 0.1,
  moment = c("ipw", "subclass"),
  n_subclass = 5L,
  pruning = FALSE,
  prune.rule = c("weight", "quantile"),
  prune.w = 0.3,
  prune.q = 0.1,
  rhat_max = 1.01,
  ess_min = 400,
  keep_draws = FALSE,
  keep_particles = FALSE
)
```

## Arguments

- mc:

  Integer. Number of MCMC iterations per chain, before burn-in and
  thinning. Ignored when posterior draws are supplied directly. Default
  5000.

- bn:

  Integer. Number of burn-in iterations to discard. Applied to
  internally generated and externally supplied draws alike, so that
  feeding a sampler's own output back in reproduces the internal path
  exactly. Default 1000.

- thin:

  Integer. Thinning interval applied after the burn-in. Default 2.

- chains:

  Integer. Number of Markov chains. Default 4. R-hat compares chains
  against each other, so at least four are recommended.

- init:

  Starting values for the samplers, or NULL (default) to place them
  automatically at a quick frequentist fit and disperse the remaining
  chains around it. A list of `chains` numeric vectors gives both models
  the same starting values, which is only meaningful when they have the
  same number of coefficients. A list with components `outcome` and
  `ps`, each itself a list of `chains` numeric vectors, gives each model
  its own. The length of each vector is checked when the sampler runs,
  since the number of coefficients is not known until the design
  matrices have been built.

- lambda_max:

  Numeric. Largest value of the tilting parameter the sweep will
  consider, written \\\bar{\lambda}\\ in the paper. Default 10.

- n_steps:

  Integer. Number of steps the sweep takes to travel from zero to
  `lambda_max`, written \\T\\ in the paper. The increment between
  successive values of the tilting parameter is `lambda_max / n_steps`.
  Default 1000, giving the increment of 0.01 used in the paper. This is
  the resolution of the grid of Algorithm 2 only; the iteration cap of
  Algorithm 1 is `newton_steps`.

- newton_steps:

  Integer. Largest number of Newton iterations Algorithm 1 may take when
  solving equation (3.9) for the tilting parameter, which is the work
  `drbayes_pc(method = "is")` does in place of the sweep. Default 100.
  The iteration is one dimensional and reaches the solution in a handful
  of steps whenever one exists, so the cap is there to end a search that
  is not converging rather than to be reached.

- smoothing:

  Numeric. Kernel smoothing coefficient of Liu and West (2001), written
  \\a\\ in the paper. A number in (0, 1\]. Default 0.99.

- tol:

  Numeric. Tolerance on the absolute posterior mean of the moment
  condition, which is the stopping rule for the sweep. NULL, the
  default, sets it to the Monte Carlo standard error of that mean over
  the retained draws: the mean is itself an average over a finite number
  of draws, so it cannot be resolved below its own sampling error, and
  asking for less than that is asking the sweep to chase noise. Give a
  positive number to fix the tolerance instead.

- ridge:

  Numeric. Small value added to the denominator of the Newton step for
  the tilting parameter, both the step Algorithm 1 iterates and the one
  that chooses the initial direction of the sweep. It keeps that step
  defined when every draw puts the moment condition at zero, which would
  otherwise be a division of zero by zero. Default 1e-5. Both
  denominators are an average of \\B_n^2\\ over the draws, the same
  average at \\\lambda = 0\\ and weighted by \\\exp(\lambda B_n)\\
  thereafter, so one value of `ridge` means the same thing in both. That
  average shrinks as the sample grows, since the moment condition itself
  does, and at a few hundred thousand observations the default is a
  noticeable share of it: that shortens the Newton steps and costs a few
  more iterations rather than moving the tilting parameter they converge
  on, and a smaller `ridge` removes the effect.

- ess_frac:

  Numeric in (0, 1\]. Smallest share of the posterior draws that has to
  be carrying the tilted posterior for the coupling to be reported as
  converged. Default 0.1. Meeting the moment condition is not on its own
  evidence that anything was learned: the importance weighted mean of
  \\B_n\\ tends to its smallest value as the tilting parameter runs off
  to minus infinity, so weights that have collapsed onto one draw can
  sit well inside `tol` while describing a posterior of a single point,
  whose standard deviation is zero and whose credible interval has zero
  width. A fit resting on fewer than `ess_frac` of the draws is
  therefore reported as not converged whatever the moment condition
  says. The quantity compared against it is the effective sample size of
  the importance weights under `method = "is"`, and the smallest number
  of distinct particles to survive a resampling step under the sweep,
  which resamples to equal weights and so has no importance weights of
  its own. A tenth is the usual rule of thumb for an importance sample
  that has been spent. Note that this is a different quantity from
  `ess_min`, which is about the Markov chains that produced the draws in
  the first place.

- moment:

  One of `"ipw"` (default) or `"subclass"`. Which form of the moment
  condition the tilting has to satisfy: the inverse probability weighted
  one of equation (3.4), or the propensity score subclassification one
  of equation (F.2). Subclassification replaces each individual weight
  by the treated fraction of the stratum the observation falls in, which
  is more stable when a few estimated propensity scores are close to 0
  or 1.

- n_subclass:

  Integer. Number of propensity score strata used when
  `moment = "subclass"`, written \\K\\ in Appendix F. Default 5,
  following Cochran (1968). The strata hold equal numbers of
  observations, so a large `n_subclass` leaves few in each and makes a
  stratum with no treated or no control units, which cannot be used,
  more likely.

- pruning:

  Logical. Whether to apply the sample pruning of Section 5.3.1, which
  discards the particles carrying the smallest weights at each step of
  the sweep so that the survivors concentrate where the moment condition
  holds. Default FALSE. It reduces the residual bias left when only the
  propensity score model is correctly specified, at the cost of
  discarding part of the particle cloud at every step, and it also ends
  the sweep much sooner: see the note on the pruning threshold below.

- prune.rule:

  Character. How pruning decides which particles to discard. `"weight"`,
  the default, discards every particle whose weight against the untilted
  posterior falls below `prune.w`, which is the rule the authors used.
  `"quantile"` instead discards a fixed fraction `prune.q` at every
  step, whatever the weights look like.

- prune.w:

  Numeric, positive. Weight below which a particle is discarded, as a
  multiple of \\1/S\\, the weight a particle would carry if all \\S\\ of
  them were equal. Default 0.3, so a particle is discarded once it
  carries less than three tenths of an equal share. Below about 0.12 the
  threshold is never reached before the sweep stops, which makes
  `pruning = TRUE` a no-op; above about 0.5 it discards so much that the
  sweep ends after one or two steps. Expressing the threshold relative
  to \\1/S\\ keeps its meaning when the number of draws changes. The
  weight tested is \\\exp(\lambda B_n)\\ at the tilting parameter
  reached so far, as in section 5.3.1, not the incremental weight of a
  single step: the increments are all within a fraction of a percent of
  \\1/S\\, so an absolute threshold on those would never discard
  anything. The paper does not state the value of the threshold, so the
  default is this package's choice; treat it as something to vary.

- prune.q:

  Numeric in \[0, 1). Fraction of the particles discarded at each step,
  smallest weight first, when `prune.rule = "quantile"`. Default 0.1.
  `prune.q = 0` leaves the sweep exactly as `pruning = FALSE` does. The
  particles it removes are those at one tail of the moment condition, so
  discarding a fixed fraction of them moves its posterior mean much
  further than an increment of the tilting parameter does, which ends
  the sweep after very few steps; the note below says how much further.

- rhat_max:

  Numeric. Largest acceptable rank-normalised split R-hat. Default 1.01,
  the threshold recommended by Vehtari et al. (2021).

- ess_min:

  Numeric. Smallest acceptable bulk and tail effective sample size.
  Default 400, again following Vehtari et al. (2021), who also note that
  below this the R-hat estimate is itself unreliable.

- keep_draws:

  Logical. Whether to return the per-chain posterior draws of the two
  models' coefficients, as the iterations by chains by parameters arrays
  the samplers produce. They are what
  [`rank_plot()`](https://t-momozaki.github.io/DRBayes/reference/rank_plot.md)
  and `plot(fit, type = "rank")` need, since a rank histogram is a
  statement about chains and the tilted particles have none. Default
  `FALSE`.

- keep_particles:

  Logical. Whether to return the tilted particles themselves. Default
  FALSE. See the note on memory below.

## Value

An object of class `"drbayes_control"`, a named list with one element
per argument above.

## Details

**Why these are settings rather than constants:**

`lambda_max`, `n_steps` and `smoothing` are the tuning constants of
Algorithm 2 of Orihara, Momozaki and Sugasawa (2025). The paper leaves
the choice of \\\bar{\lambda}\\, \\T\\ and the smoothing coefficient
\\a\\ to the user, so hard-coding them would hide exactly the knobs the
methodology asks you to tune. The defaults reproduce the settings used
in the paper and are a reasonable starting point, not a recommendation
to leave them alone.

`newton_steps` and `ess_frac` have no counterpart in the paper, which
asks only that Algorithm 1 be iterated "until convergence" and does not
say what to do when it does not converge. They are this package's answer
to that: a cap on the search, and a floor below which the sample the
search arrived at is too thin to be called a posterior. `lambda_max`
bounds Algorithm 1 as well as the sweep, for the same reason.

**Choosing the smoothing coefficient:**

The sweep resamples the particles at every step and then applies the Liu
and West (2001) kernel, which shrinks each particle towards the cloud
mean and adds normal noise carrying the variance that shrinking removed.
One application is close to harmless. Applying it `n_steps` times is
not: repeated shrink-and-jitter pulls the particle cloud towards a
multivariate normal with the right first two moments and the wrong
shape. That is a poor fit for the heavy-tailed posterior with a spike at
zero that the horseshoe priors produce, so a user fitting
[`bayes_lm_hs`](https://t-momozaki.github.io/DRBayes/reference/bayes_lm_hs.md),
[`bayes_logit_hs`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit_hs.md)
or
[`bayes_probit_hs`](https://t-momozaki.github.io/DRBayes/reference/bayes_probit_hs.md)
may want fewer and larger steps in the tilting parameter, for example
`n_steps = 200`, or a smoothing coefficient closer to 1, which shrinks
less at each application.

**The pruning threshold:**

Section 5.3.1 of the paper describes sample pruning as discarding the
samples with small sampling weights, and stops there: it gives no
threshold, and no rule for deriving one. `prune.q` is this package's
reading of that sentence, discarding a fixed fraction of the particles
at each step, and its default of 0.1 is a choice made here rather than a
value reported in the paper.

Pruning ends the sweep far sooner than its name suggests, and it is
worth being clear about why. The particles it discards are those at one
tail of the moment condition, the tail the tilting is pushing away from,
so deleting `prune.q` of the cloud moves the posterior mean of \\B_n\\
by a noticeable fraction of that mean's distance from zero, while one
increment of the tilting parameter moves it by far less: at the default
settings one pruning step is worth some tens of ordinary steps. Step 4
of Algorithm 2 evaluates the constraint on the particles it currently
holds, which are the pruned ones, so the sweep typically stops after one
or two steps at a tilting parameter close to zero. How many steps it
takes is then set by `prune.q` and hardly at all by `n_steps`.

The direction of that is what Section 5.3.1 asks for, since pruning
exists precisely so that \\\lambda\\ need not grow large enough to
violate condition (C.3) of Appendix A. What it costs is that the
constraint has been met largely by deleting particles rather than by
tilting them: the `lambda` reported in `$smc` is a poor measure of how
much tilting was done, and the returned cloud is a truncated one, so its
credible intervals are narrower than those of the tilted posterior by an
amount that reflects the deletion and not information gained. Fit both
ways. A result that changes materially between `pruning = FALSE` and
`pruning = TRUE`, or between `prune.q = 0.05` and `prune.q = 0.2`, is
telling you that the tilting is being carried by a small part of the
particle cloud, which is worth knowing.

**Keeping the particles:**

`keep_particles = TRUE` returns the tilted particles, one draws by
parameters matrix for the outcome model and one for the propensity score
model. They are useful for inspecting what the tilting did to the
coefficients, and expensive in memory: with the default settings and a
moderate number of covariates they are several times the size of
everything else the fit returns.

## References

Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
Causal Inference via Posterior Coupling. arXiv preprint
arXiv:2506.04868.

Cochran, W. G. (1968). The effectiveness of adjustment by
subclassification in removing bias in observational studies. Biometrics,
24(2), 295-313.

Liu, J., & West, M. (2001). Combined parameter and state estimation in
simulation-based filtering. In A. Doucet, N. de Freitas, & N. Gordon
(Eds.), Sequential Monte Carlo Methods in Practice (pp. 197-223).
Springer.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Buerkner, P.-C.
(2021). Rank-normalization, folding, and localization: an improved R-hat
for assessing convergence of MCMC. Bayesian Analysis, 16(2), 667-718.

## See also

[`drbayes_pc`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md),
which consumes the result.

## Examples

``` r
# The defaults, ready to be modified
drbayes_control()
#> Control settings for drbayes_pc()
#> 
#>   MCMC         mc = 5000, bn = 1000, thin = 2, chains = 4, init = automatic
#>   Tilting      lambda_max = 10, n_steps = 1000 (step 0.01), smoothing = 0.99
#>                tol = Monte Carlo standard error, ridge = 1e-05, ess_frac = 0.1
#>                moment = ipw, pruning = off
#>                newton_steps = 100 (importance sampling only)
#>   Diagnostics  rhat_max = 1.01, ess_min = 400
#>   Output       keep_draws = FALSE, keep_particles = FALSE

# Horseshoe priors give a heavy-tailed posterior that repeated kernel
# smoothing distorts, so take fewer and larger steps in lambda.
ctrl <- drbayes_control(n_steps = 200)
ctrl$n_steps
#> [1] 200

# A plain list of overrides is accepted wherever a control object is, so
# control = list(n_steps = 200) means the same thing.

set.seed(1)
n <- 200
dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
dat$A <- rbinom(n, 1, plogis(0.2 + 0.5 * dat$X1 - 0.3 * dat$X2))
dat$Y <- 1 + 1.5 * dat$A + 0.8 * dat$X1 + 0.6 * dat$X2 + rnorm(n)

fit <- drbayes_pc(
  outcome.formula = Y ~ A + X1 + X2,
  ps.formula      = A ~ X1 + X2,
  data            = dat,
  outcome.model   = bayes_lm_hs,
  ps.model        = bayes_logit,
  control         = drbayes_control(mc = 1000, bn = 200, n_steps = 60)
)
#> Formula processing complete:
#>   Observations used: 200 out of 200
#>   Treatment group: 107 | Control group: 93
#>   Outcome model variables: 3
#>   Propensity score variables: 2
#> Auto-detected continuous outcome: using family='gaussian'
#> Generating posterior samples using provided models...
#> Using 1600 posterior samples for analysis...
fit$smc$n_steps
#> [1] 5
```

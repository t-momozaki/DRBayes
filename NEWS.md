# DRBayes 0.1.1

First public release, of the reference implementation of Orihara, Momozaki and
Sugasawa (2025), "Bayesian Doubly Robust Causal Inference via Posterior
Coupling" ([doi:10.48550/arXiv.2506.04868](https://doi.org/10.48550/arXiv.2506.04868)).

## New features

* `drbayes_pc(Y ~ A + X1, A ~ X1, data)` fits with no further arguments.
  Given the family and the link the sampling algorithm follows, so the only
  choice left is the prior: `outcome.prior` and `ps.prior` select `"normal"` or
  `"horseshoe"` and the sampler follows. Supplying a sampler of your own still
  works and still takes precedence.
* `method` chooses between the two algorithms of the paper: `"smc"` (default)
  is Algorithm 2, `"is"` is the single importance sampling step of Algorithm 1.
* `drbayes_control()` collects the tuning constants of Algorithm 2, which were
  hard-coded before, together with the stopping tolerance, the pruning rule of
  Section 5.3.1, the diagnostic thresholds and `keep_particles`.
* `drbayes_sensitivity()` implements the sensitivity analysis of Algorithm 3
  for unmeasured confounding, and `drbayes_select()` the confounder selection
  of Algorithm 4.
* `convergence_diagnostics()` reports rank-normalised split R-hat, its folded
  counterpart, bulk and tail effective sample size and the Monte Carlo standard
  error, following Vehtari et al. (2021), and `rank_plot()` draws the rank
  histograms the paper proposes in place of trace plots. All six samplers now
  take `chains` and `init` and return an iterations by chains by parameters
  array. `drbayes_pc()` checks both posteriors before spending the sweep on
  draws that cannot support it; `diagnostics = "error"` refuses to continue and
  `"none"` skips the check.
  These statistics apply to Markov chain draws, so to the built-in samplers and
  to external MCMC output, and not to the sequential Monte Carlo particles in
  `$pc`, which are resampled and carry no time ordering.
* `print()`, `summary()` and `plot()` methods for the fitted object, giving a
  one screen summary of the two estimands, the state of the sweep and the worst
  convergence diagnostics.
* `bayes_stan()` fits the same six models with CmdStan through `cmdstanr` and
  returns draws in the same layout. It is an option, never the default: CRAN's
  machines have no CmdStan, so the binaries CRAN distributes would carry no
  compiled models. A missing CmdStan is an error naming both installation
  steps, never a silent substitution.
* `drbayes_bb()` and `drbayes_pc()` take the same `(outcome.formula,
  ps.formula, data)` interface and share one implementation of it, so the two
  estimators cannot disagree about how a formula is read. `drbayes_bb()` also
  gains `trim` to bound the propensity score.
* `convert_to_array()` is exported, for reshaping the output of a replicated
  simulation study.

## Renamed functions

The exported names were brought to one scheme and moved away from dots, since a
dot is how R marks an S3 method and `DRBayes.PC` could not be told apart from a
`print` method for class `DRBayes`. The old names still work and forward to
their replacements with a warning, and will be removed in the next release.

| Old | New |
| --- | --- |
| `DRBayes.PC()` | `drbayes_pc()` |
| `DRBayes.BB()`, `DRBayes.BB.nleqslv()` | `drbayes_bb()` |
| `B.LM()` | `bayes_lm()` |
| `B.Logit()` | `bayes_logit()` |
| `B.Probit()` | `bayes_probit()` |
| `HS.LM()` | `bayes_lm_hs()` |
| `HS.Logit()` | `bayes_logit_hs()` |
| `HS.Probit()` | `bayes_probit_hs()` |

`DRBayes.BB.nleqslv()` solved the same two weighted regressions with a general
nonlinear root finder and agreed with `DRBayes.BB()` to about 1e-08, the
difference being only the root finder's stopping tolerance, so it forwards to
`drbayes_bb()` rather than being kept as a second implementation. `nleqslv` and
`parallelly` are no longer dependencies.

## Speed and memory

The estimates are unchanged. Every rewrite below was checked against the
implementation it replaced on the same draws, and the two agree bit for bit,
not to a tolerance: at seventy-eight combinations of propensity score link,
outcome link and sample size for equation (3.4), and at twenty-one for equation
(F.2). The average treatment effect a whole fit returns is identical either
way.

* **The moment conditions are evaluated in compiled code.** They are the whole
  cost of the sequential Monte Carlo sweep, which evaluates one of them once a
  step. The matrix products go to the same BLAS R would have called; what
  changed is everything after them, which R performed by building about seven
  matrices of draws by observations per block and which now runs in one pass
  over two buffers. Whole fits, timed against the R implementation alternately
  so that a changing machine load cannot favour either and taking the minimum
  of three:

  | | steps | R | compiled | |
  | --- | ---: | ---: | ---: | ---: |
  | `moment = "ipw"`, n = 500 | 4 | 4.02 s | 3.77 s | 1.07x |
  | `moment = "ipw"`, n = 2000 | 99 | 37.52 s | 21.38 s | 1.76x |
  | `moment = "ipw"`, n = 5000 | 323 | 212.57 s | 105.61 s | 2.01x |
  | `moment = "subclass"`, n = 2000 | 205 | 25.96 s | 11.57 s | 2.24x |

  The gain is proportional to the number of steps the sweep takes, so it is
  largest exactly where fits are slowest. At n = 500 the sweep converges in
  four steps and there is almost nothing to save.

* **The propensity score posterior mean is accumulated block by block.** The
  strata of Appendix F are cut on it, and it is accumulated a block of draws at
  a time rather than by holding the linear predictors and the scores at full
  size: the working memory no longer grows with the number of draws.

* **Peak memory.** What the package holds, measured as R's own heap high-water
  mark, is 126 MB at n = 500 and 140 MB at both n = 2000 and n = 8000, so it no
  longer grows with the sample size. Measured instead as the peak resident set
  of the whole process, the compiled path is below the R one in every case
  tried, by 23 MB at the smallest and 286 MB at n = 8000 with 2000 draws. The
  resident set is larger than R's heap because the C library does not return
  freed pages to the system; that part is outside the package's control.

## Reproducibility of the numbers

* **A fit is reproducible from its seed on one machine, and reproducible to
  rounding across machines.** The compiled kernels return values identical bit
  for bit to the R implementations they replace, which is checked against them
  directly on every platform the package is tested on. What is not identical
  across machines is the matrix product underneath: R accumulates sums in
  `long double` where the build has it and in `double` where it does not, and a
  BLAS may reach a different kernel for a differently shaped product. Two
  machines will therefore agree to about the last two digits of a double rather
  than exactly. Within one machine, `seed` fixes the answer.
* **`drbayes_pc()` leaves the random number stream as it found it**, including
  the generator kind, so a fit can be dropped into an existing script without
  changing what follows it.

## Which average the treatment effect is over

* **A fit now says which estimand its interval belongs to.** `g.comp` and `pc`
  are posteriors for equation (3.7), the average of the fitted contrast over
  the covariate vectors that were observed. Section 2.1 of the paper defines
  the estimand as an average over the population those vectors came from. The
  two are different quantities, and `print()`, `summary()` and the help page
  now say which one a fit reports.

  The difference is in the spread, not the centre. Averaging over the observed
  covariates treats their distribution as known, so the interval does not
  carry the uncertainty of not knowing it. Weighting the observations by a
  Dirichlet, as a Bayesian bootstrap would, adds exactly `Var_i(c_i)/n` to the
  posterior variance, where `c_i` is the unit level contrast, and leaves the
  posterior mean alone. The interval widens by a factor of
  `sqrt(1 + Var_i(c_i)/(n * V))`, which is a quantity a reader can compute for
  their own fit; the vignette shows how.

  It matters only when the fitted contrast varies from unit to unit. On a unit
  level effect of `2 + 1.5 * X1` the factor is 1.11 at n = 500 and 1.21 at
  n = 2000. With no treatment interaction and an identity link it is exactly 1.
  With a logit link and no interaction it is 1.001: the link does make the
  contrast vary, and not by enough to notice.

  The default is unchanged, and deliberately. Equation (3.7) is what the paper
  computes, and a reference implementation should compute it. Theorem 1
  concerns the posterior mean, which the choice of weights does not move.

## Documented behaviour worth knowing

* `drbayes_bb()` targets the mixed estimand of Saarela et al. (2016) rather
  than the superpopulation ATE, which Section 4.3 of the paper is explicit
  about and the previous documentation did not mention.
* Progress reporting goes through `message()` behind a `verbose` argument, so
  it can be silenced. It used `cat()`, which cannot.
* A missing or non-finite value anywhere in the design made the two logistic
  samplers loop forever rather than fail, because the Polya-Gamma rejection
  sampler never returns on a `NaN`. Finite input is now required.
* An intercept-only design, which `ps.formula = A ~ 1` asks for, crashed five
  of the six samplers on a mismatch in the length of the dimnames. The three
  samplers under a normal prior now fit it; the three horseshoe samplers refuse
  it, there being nothing left for a shrinkage prior to shrink.

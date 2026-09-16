# Package index

## Estimating a treatment effect

The two estimators.
[`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
is the posterior coupling of the paper, which returns an explicit
posterior for the average treatment effect;
[`drbayes_bb()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_bb.md)
is the Bayesian bootstrap estimator of Saarela et al., kept for
comparison and targeting a different estimand.

- [`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
  : Bayesian Doubly Robust Causal Inference via Posterior Coupling
- [`drbayes_bb()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_bb.md)
  : Bayesian Bootstrap Doubly Robust Estimation of the Average Treatment
  Effect
- [`drbayes_control()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_control.md)
  : Tuning Parameters for Posterior Coupling

## Looking at a fit

What a fit carries and how to read it. The tilted and untilted
posteriors are both returned, so the two can always be compared.

- [`DRBayes-methods`](https://t-momozaki.github.io/DRBayes/reference/DRBayes-methods.md)
  : Methods for posterior coupling fits
- [`print(`*`<DRBayes>`*`)`](https://t-momozaki.github.io/DRBayes/reference/print.DRBayes.md)
  : Print a posterior coupling fit
- [`summary(`*`<DRBayes>`*`)`](https://t-momozaki.github.io/DRBayes/reference/summary.DRBayes.md)
  [`print(`*`<summary.DRBayes>`*`)`](https://t-momozaki.github.io/DRBayes/reference/summary.DRBayes.md)
  : Summarise a posterior coupling fit
- [`plot(`*`<DRBayes>`*`)`](https://t-momozaki.github.io/DRBayes/reference/plot.DRBayes.md)
  : Plot a posterior coupling fit
- [`print(`*`<drbayes_control>`*`)`](https://t-momozaki.github.io/DRBayes/reference/print.drbayes_control.md)
  : Print Posterior Coupling Control Settings

## Posterior samplers

The models whose posteriors are coupled. Each returns an iterations by
chains by parameters array and can be used on its own as a Bayesian
regression. The horseshoe variants take an `unshrunk` argument, which is
how the treatment effect is kept out of the shrinkage.

- [`bayes_lm()`](https://t-momozaki.github.io/DRBayes/reference/bayes_lm.md)
  : Bayesian Linear Regression using MCMC
- [`bayes_logit()`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit.md)
  : Bayesian Logistic Regression with Polya-Gamma Latent Variables
- [`bayes_probit()`](https://t-momozaki.github.io/DRBayes/reference/bayes_probit.md)
  : Bayesian Probit Regression with Latent Variables
- [`bayes_lm_hs()`](https://t-momozaki.github.io/DRBayes/reference/bayes_lm_hs.md)
  : Bayesian Linear Regression with Horseshoe Prior
- [`bayes_logit_hs()`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit_hs.md)
  : Bayesian Logistic Regression with Horseshoe Prior
- [`bayes_probit_hs()`](https://t-momozaki.github.io/DRBayes/reference/bayes_probit_hs.md)
  : Bayesian Probit Regression with Horseshoe Prior
- [`bayes_stan()`](https://t-momozaki.github.io/DRBayes/reference/bayes_stan.md)
  : Fit One of the Package Models with Stan Instead of the Gibbs Sampler

## Convergence diagnostics

Rank normalised R-hat and effective sample size, after Vehtari et
al. (2021). These are statements about Markov chain draws: they apply to
the samplers above, and not to the sequential Monte Carlo particles,
which are resampled and carry no order.

- [`convergence_diagnostics()`](https://t-momozaki.github.io/DRBayes/reference/convergence_diagnostics.md)
  : Convergence diagnostics for posterior draws
- [`rank_plot()`](https://t-momozaki.github.io/DRBayes/reference/rank_plot.md)
  : Rank plots of posterior draws

## Sensitivity to unmeasured confounding

Algorithm 3 of the paper, reweighting the fitted posterior under an
assumed distribution for the bias rather than refitting.

- [`drbayes_sensitivity()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_sensitivity.md)
  : Sensitivity analysis for unmeasured confounding
- [`xi_triangular()`](https://t-momozaki.github.io/DRBayes/reference/xi_triangular.md)
  : Triangular sensitivity distributions
- [`print(`*`<DRBayes_sensitivity>`*`)`](https://t-momozaki.github.io/DRBayes/reference/print.DRBayes_sensitivity.md)
  : Print a sensitivity analysis
- [`summary(`*`<DRBayes_sensitivity>`*`)`](https://t-momozaki.github.io/DRBayes/reference/summary.DRBayes_sensitivity.md)
  [`print(`*`<summary.DRBayes_sensitivity>`*`)`](https://t-momozaki.github.io/DRBayes/reference/summary.DRBayes_sensitivity.md)
  : Summarise a sensitivity analysis

## Confounder selection

Algorithm 4, selecting from a shrinkage posterior and updating only the
selected coefficients through a modified moment condition.

- [`drbayes_select()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_select.md)
  : Confounder Selection with Posterior Coupling
- [`print(`*`<DRBayes_select>`*`)`](https://t-momozaki.github.io/DRBayes/reference/print.DRBayes_select.md)
  : Print a confounder selection fit

## Simulation studies

The data generating process of Section 5 and the machinery for
replicating it, including the seeding that makes a study reproducible
across any number of workers.

- [`generate_dataset()`](https://t-momozaki.github.io/DRBayes/reference/generate_dataset.md)
  : Generate Synthetic Dataset for Causal Inference Simulation
- [`run_parallel_simulation()`](https://t-momozaki.github.io/DRBayes/reference/run_parallel_simulation.md)
  : Run Parallel Simulation Study
- [`convert_to_array()`](https://t-momozaki.github.io/DRBayes/reference/convert_to_array.md)
  : Convert Simulation Results to Array Format

## The package

- [`DRBayes`](https://t-momozaki.github.io/DRBayes/reference/DRBayes-package.md)
  [`DRBayes-package`](https://t-momozaki.github.io/DRBayes/reference/DRBayes-package.md)
  : DRBayes: Bayesian Doubly Robust Causal Inference Methods

## Superseded names

The names used before the release. Each forwards to its replacement and
warns once per session.

- [`DRBayes.PC()`](https://t-momozaki.github.io/DRBayes/reference/DRBayes-deprecated.md)
  [`DRBayes.BB()`](https://t-momozaki.github.io/DRBayes/reference/DRBayes-deprecated.md)
  [`DRBayes.BB.nleqslv()`](https://t-momozaki.github.io/DRBayes/reference/DRBayes-deprecated.md)
  [`B.LM()`](https://t-momozaki.github.io/DRBayes/reference/DRBayes-deprecated.md)
  [`B.Logit()`](https://t-momozaki.github.io/DRBayes/reference/DRBayes-deprecated.md)
  [`B.Probit()`](https://t-momozaki.github.io/DRBayes/reference/DRBayes-deprecated.md)
  [`HS.LM()`](https://t-momozaki.github.io/DRBayes/reference/DRBayes-deprecated.md)
  [`HS.Logit()`](https://t-momozaki.github.io/DRBayes/reference/DRBayes-deprecated.md)
  [`HS.Probit()`](https://t-momozaki.github.io/DRBayes/reference/DRBayes-deprecated.md)
  : Functions renamed in DRBayes 0.1.0

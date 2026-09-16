# DRBayes: Bayesian Doubly Robust Causal Inference via Posterior Coupling

Reference documentation and articles:
<https://t-momozaki.github.io/DRBayes/>

## Overview

DRBayes implements Bayesian doubly robust methods for causal inference
using the novel posterior coupling approach. This methodology combines
the robustness of doubly robust estimation with the uncertainty
quantification benefits of Bayesian inference, providing a principled
framework for estimating average treatment effects.

### Key Features

- **Formula Interface**: Intuitive R formula syntax compatible with
  [`lm()`](https://rdrr.io/r/stats/lm.html),
  [`glm()`](https://rdrr.io/r/stats/glm.html), and other standard
  functions
- **Doubly Robust Estimation**: Consistent estimation even when either
  the outcome model or propensity score model is misspecified
- **Posterior Coupling**: Novel methodology that avoids the feedback
  problem in traditional doubly robust methods
- **Bayesian Uncertainty Quantification**: Full posterior distributions
  for treatment effects
- **Multiple Outcome Types**: Support for continuous (linear regression)
  and binary outcomes (logistic/probit regression)
- **Automatic Detection**: Smart detection of outcome family and link
  functions based on data
- **Flexible Modeling**: Built-in functions with standard and horseshoe
  priors plus support for external Bayesian software
- **Missing Data Handling**: Comprehensive support for missing value
  treatment
- **Software Integration**: Seamless integration with Stan, brms, JAGS,
  and other Bayesian platforms

## Installation

You can install the development version of DRBayes from GitHub:

``` r

# Install devtools if you haven't already
if (!require(devtools)) install.packages("devtools")

# Install DRBayes with vignettes
devtools::install_github("t-momozaki/DRBayes", build_vignettes = TRUE)
```

### Dependencies

DRBayes installs the following R packages, all of them from CRAN:

- `Rcpp` and `RcppArmadillo` (the compiled moment conditions)
- `mvtnorm` (Multivariate normal sampling)
- `extraDistr` (Inverse gamma distribution)
- `RcppTN` (Truncated normal sampling for probit models)
- `pgdraw` (Polya-Gamma sampling for logistic models)
- `future`, `furrr` (Parallel replication in
  [`run_parallel_simulation()`](https://t-momozaki.github.io/DRBayes/reference/run_parallel_simulation.md))
- `stats`, `graphics`, `utils` (Standard R functions)

**Optional.** `cmdstanr` and `instantiate` are needed only by
[`bayes_stan()`](https://t-momozaki.github.io/DRBayes/reference/bayes_stan.md),
which fits the same six models with Stan instead of the built-in Gibbs
samplers; `ggplot2` and `dplyr` only by the vignettes. Posterior draws
from Stan, brms, JAGS or any other software can be passed to
[`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
as a matrix, which needs none of those packages installed.

## Quick Start

### Example 1: Continuous Outcome with Formula Interface

``` r

library(DRBayes)

# Generate example data with continuous outcome
set.seed(1234)
n <- 300
data <- data.frame(
  X1 = rnorm(n),
  X2 = rnorm(n),
  X3 = rnorm(n),
  X4 = rnorm(n)
)

# Generate treatment with confounding
data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$X1 - 0.3*data$X2))

# Generate continuous outcome with treatment effect and interactions
data$Y <- 1 + 1.5*data$A + 0.8*data$X1 + 0.6*data$X2 + 0.4*data$A*data$X1 + rnorm(n)

# Fit Bayesian doubly robust model using formula interface
result_continuous <- drbayes_pc(
  outcome.formula = Y ~ A + X1 + X2 + A:X1,      # Outcome model with interaction
  ps.formula = A ~ X1 + X2 + X3,                 # Propensity score model
  data = data,
  family = "gaussian",                           # Auto-detected if NULL
  outcome.model = bayes_lm,                          # Bayesian linear regression
  ps.model = bayes_logit,                            # Bayesian logistic regression
  mc = 3000, bn = 1000, thin = 2
)

# View results
cat("ATE estimate:", round(mean(result_continuous$pc), 3), "\n")
cat("Posterior SD:", round(sd(result_continuous$pc), 3), "\n")
cat("95% CI:", round(quantile(result_continuous$pc, c(0.025, 0.975)), 3), "\n")

# Compare with g-computation
cat("G-computation ATE:", round(mean(result_continuous$g.comp), 3), "\n")
```

### Example 2: Binary Outcome with Automatic Detection

``` r

# Generate binary outcome
data$Y_binary <- rbinom(n, 1, plogis(0.5 + 0.8*data$A + 0.4*data$X1 + 0.3*data$X2))

# Fit model with automatic family/link detection
result_binary <- drbayes_pc(
  outcome.formula = Y_binary ~ A + X1 + X2,
  ps.formula = A ~ X1 + X2 + X3,
  data = data,
  # family and link auto-detected as "binomial" and "logit"
  outcome.model = bayes_logit,
  ps.model = bayes_logit,
  mc = 3000, bn = 1000, thin = 2
)

# Results represent probability differences for binary outcomes
cat("Family used:", result_binary$family, "\n")
cat("Link function:", result_binary$link, "\n")
cat("ATE (probability difference):", round(mean(result_binary$pc), 3), "\n")
```

### Example 3: Using Horseshoe Priors for High-Dimensional Data

``` r

# Generate high-dimensional data
set.seed(1234)
n <- 200
p <- 50  # Many covariates
data_hd <- data.frame(matrix(rnorm(n * p), n, p))
names(data_hd) <- paste0("X", 1:p)

# Only first few covariates are truly relevant
data_hd$A <- rbinom(n, 1, plogis(0.3*data_hd$X1 + 0.2*data_hd$X2 - 0.1*data_hd$X3))
data_hd$Y <- 2 + 1.5*data_hd$A + 0.8*data_hd$X1 + 0.5*data_hd$X2 + rnorm(n)

# Use horseshoe priors for regularization. The treatment main effect and its
# interactions are exempt from shrinkage automatically.
result_horseshoe <- drbayes_pc(
  outcome.formula = Y ~ A + .,  # Include all covariates
  ps.formula = A ~ . - Y,       # All covariates except outcome
  data = data_hd,
  outcome.prior = "horseshoe",
  ps.prior = "horseshoe",
  mc = 2000, bn = 500, thin = 2
)

cat("Horseshoe ATE:", round(mean(result_horseshoe$pc), 3), "\n")
```

The horseshoe logistic sampler mixes slowly on 50 covariates and 200
observations, so this short run warns that the propensity score model
has not met `R-hat < 1.01`. That is the check doing its job;
`result_horseshoe$diagnostics` names the coefficients, and raising `mc`
is the remedy.

### Example 4: Advanced Formula Usage and Missing Data

``` r

# Create data with missing values and complex relationships
data_missing <- data
data_missing$X1[sample(n, 20)] <- NA
data_missing$X2[sample(n, 15)] <- NA

# Complex formulas with transformations and polynomials
result_advanced <- drbayes_pc(
  outcome.formula = Y ~ A + log(abs(X1) + 1) + I(X2^2) + poly(X3, 2) + A:X1,
  ps.formula = A ~ X1 + X2 + X3 + I(X1^2) + X2:X3,
  data = data_missing,
  na.action = "na.omit",                         # Handle missing values
  outcome.model = bayes_lm,
  ps.model = bayes_logit,
  mc = 2000, bn = 500
)

cat("Original sample size:", n, "\n")
cat("Analysis sample size:", result_advanced$data_info$n_observations, "\n")
cat("Missing observations:", result_advanced$data_info$missing_observations, "\n")
cat("Advanced formula ATE:", round(mean(result_advanced$pc), 3), "\n")
```

### Example 5: Bayesian Bootstrap Alternative

``` r

# Alternative: Bayesian Bootstrap approach (Saarela et al. 2016).
# Note that this targets Saarela's mixed estimand rather than the
# superpopulation ATE; see ?drbayes_bb.
result_bb <- drbayes_bb(
  outcome.formula = Y ~ A + X1 + X2 + X3 + A:X1,
  ps.formula      = A ~ X1 + X2 + X3,
  data            = data,
  num_iterations  = 1000,
  family          = "gaussian"
)

cat("Bayesian Bootstrap ATE:", round(mean(result_bb), 3), "\n")
```

## Main Functions

### Core Functions

- **[`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)**:
  Main function for Bayesian doubly robust estimation via posterior
  coupling with formula interface
- **[`drbayes_bb()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_bb.md)**:
  Bayesian bootstrap approach (Saarela et al. 2016)
- **[`drbayes_control()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_control.md)**:
  Tuning parameters for the sequential Monte Carlo sweep and the
  convergence check
- **[`drbayes_sensitivity()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_sensitivity.md)**:
  Sensitivity analysis for unmeasured confounding. Reweights the draws
  of an existing fit and couples them again; no model is refitted
- **[`drbayes_select()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_select.md)**:
  Confounder selection on the propensity score model under shrinkage
  priors, followed by posterior coupling on the selected set

### Built-in Model Functions

These are the samplers
[`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
chooses between. Call them directly only when you want the draws for
their own sake, or want to hand a sampler the main function does not
provide.

**Normal priors:** -
**[`bayes_lm()`](https://t-momozaki.github.io/DRBayes/reference/bayes_lm.md)**:
Bayesian linear regression with conjugate priors (for continuous
outcomes) -
**[`bayes_logit()`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit.md)**:
Bayesian logistic regression using Polya-Gamma data augmentation (for
binary outcomes/treatment) -
**[`bayes_probit()`](https://t-momozaki.github.io/DRBayes/reference/bayes_probit.md)**:
Bayesian probit regression using latent variable approach (for binary
outcomes/treatment)

**Horseshoe Priors (for high-dimensional data):** -
**[`bayes_lm_hs()`](https://t-momozaki.github.io/DRBayes/reference/bayes_lm_hs.md)**:
Horseshoe prior linear regression -
**[`bayes_logit_hs()`](https://t-momozaki.github.io/DRBayes/reference/bayes_logit_hs.md)**:
Horseshoe prior logistic regression -
**[`bayes_probit_hs()`](https://t-momozaki.github.io/DRBayes/reference/bayes_probit_hs.md)**:
Horseshoe prior probit regression

**Stan backend:** -
**[`bayes_stan()`](https://t-momozaki.github.io/DRBayes/reference/bayes_stan.md)**:
Fits any of the six models above with CmdStan and returns draws in the
same layout. Requires `cmdstanr` and CmdStan, so it is never the
default.

### Diagnostics

- **[`convergence_diagnostics()`](https://t-momozaki.github.io/DRBayes/reference/convergence_diagnostics.md)**:
  Rank-normalised split R-hat, bulk and tail effective sample size, and
  Monte Carlo standard error (Vehtari et al. 2021)
- **[`rank_plot()`](https://t-momozaki.github.io/DRBayes/reference/rank_plot.md)**:
  Rank histograms, the paper’s replacement for trace plots

### Simulation and Utility Functions

- **[`generate_dataset()`](https://t-momozaki.github.io/DRBayes/reference/generate_dataset.md)**:
  Generate synthetic datasets for simulation studies
- **[`run_parallel_simulation()`](https://t-momozaki.github.io/DRBayes/reference/run_parallel_simulation.md)**:
  Run parallel simulation studies
- **[`convert_to_array()`](https://t-momozaki.github.io/DRBayes/reference/convert_to_array.md)**:
  Collect the output of
  [`run_parallel_simulation()`](https://t-momozaki.github.io/DRBayes/reference/run_parallel_simulation.md)
  into an array

## Outcome Types and Model Families

DRBayes supports multiple outcome types with automatic detection:

| Family | Link Function | Outcome Type | Use Case | Auto-Detection |
|----|----|----|----|----|
| `"gaussian"` | `"identity"` | Continuous | Linear regression, mean differences | Non-binary Y values |
| `"binomial"` | `"logit"` | Binary (0/1) | Logistic regression, probability differences | All Y in {0,1} |
| `"binomial"` | `"probit"` | Binary (0/1) | Probit regression, probability differences | Manual specification |

**Automatic Detection**: When `family = NULL` (default), the function
automatically detects the appropriate family: - If all Y values are 0 or
1, `family = "binomial"` and `link = "logit"` - Otherwise,
`family = "gaussian"` and `link = "identity"`

## Formula Interface Features

The formula interface provides several advantages:

### Intuitive Syntax

``` r

# Similar to lm() and glm() - familiar to R users
drbayes_pc(
  outcome.formula = Y ~ A + X1 + X2 + A:X1,
  ps.formula = A ~ X1 + X2 + X3,
  data = data
)
```

### Automatic Handling

- **Design matrices**: Automatically created from formulas
- **Factor variables**: Properly handled with contrast coding
- **Interactions**: Simple syntax for interaction terms (`A:X1`, `A*X1`)
- **Transformations**: Support for
  [`I()`](https://rdrr.io/r/base/AsIs.html),
  [`log()`](https://rdrr.io/r/base/Log.html),
  [`poly()`](https://rdrr.io/r/stats/poly.html), etc.
- **Missing values**: Consistent handling across models

### Complex Modeling

``` r

# Advanced transformations and interactions
outcome.formula = Y ~ A + log(X1) + I(X2^2) + poly(X3, 3) + A:X1 + A:I(X2^2)
ps.formula = A ~ X1 + X2 + X3 + I(X1^2) + X1:X2
```

## External Software Integration

DRBayes supports posterior samples from popular Bayesian software with
the formula interface:

### Stan/cmdstanr Integration

``` r

# Extract samples from Stan fit
outcome_samples <- fit_stan$draws("beta", format = "matrix")
ps_samples <- fit_stan$draws("gamma", format = "matrix")

# Use with formula interface. Stan has already discarded its warm-up draws,
# so there is no burn-in left for DRBayes to remove: bn = 0.
result <- drbayes_pc(
  outcome.formula = Y ~ A + X1 + X2 + A:X1,
  ps.formula = A ~ X1 + X2 + X3,
  data = data,
  outcome.samples = outcome_samples,
  ps.samples = ps_samples,
  bn = 0, thin = 1
)
```

### brms Integration

``` r

# Extract fixed effects from brms
outcome_samples <- as.matrix(brms_fit)[, 1:ncol(model.matrix(outcome.formula, data))]
ps_samples <- as.matrix(brms_ps_fit)[, 1:ncol(model.matrix(ps.formula, data))]

result <- drbayes_pc(
  outcome.formula = Y ~ A + X1 + X2,
  ps.formula = A ~ X1 + X2 + X3,
  data = data,
  outcome.samples = outcome_samples,
  ps.samples = ps_samples
)
```

### JAGS/R2jags Integration

``` r

# Extract parameter samples from JAGS
outcome_samples <- jags_fit$BUGSoutput$sims.matrix[, grep("beta", colnames(...))]
ps_samples <- jags_fit$BUGSoutput$sims.matrix[, grep("gamma", colnames(...))]

result <- drbayes_pc(
  outcome.formula = Y ~ A + X1 + X2,
  ps.formula = A ~ X1 + X2 + X3,
  data = data,
  outcome.samples = outcome_samples,
  ps.samples = ps_samples
)
```

## Prior Specification

Given the family and the link, the sampling algorithm follows, so the
only choice left is the prior. `outcome.prior` and `ps.prior` make it,
and no sampler has to be named:

``` r

# Normal priors (the default) on both models
result <- drbayes_pc(Y ~ A + X1 + X2, A ~ X1 + X2 + X3, data = data)

# Horseshoe on both models, for many candidate confounders. The treatment
# effect and its interactions are exempt from shrinkage automatically.
result_hs <- drbayes_pc(
  outcome.formula = Y ~ A + .,
  ps.formula = A ~ . - Y,
  data = data_hd,
  outcome.prior = "horseshoe",
  ps.prior = "horseshoe"
)
```

The hyperparameters of whichever sampler is chosen are set through
`outcome.priors` and `ps.priors`:

``` r

# Custom priors for internal sampling
result <- drbayes_pc(
  outcome.formula = Y ~ A + X1 + X2,
  ps.formula = A ~ X1 + X2 + X3,
  data = data,
  outcome.model = bayes_lm,
  ps.model = bayes_logit,
  outcome.priors = list(
    theta_prior = 1/50,        # Prior precision for coefficients
    sigma_prior = c(2, 1)      # Inverse gamma prior for error variance
  ),
  ps.priors = list(
    theta_prior = 1/25         # Prior precision for PS coefficients
  )
)

# Horseshoe priors with custom tau
result_hs <- drbayes_pc(
  outcome.formula = Y ~ A + .,
  ps.formula = A ~ . - Y,
  data = data_hd,
  outcome.model = bayes_lm_hs,
  ps.model = bayes_logit_hs,
  outcome.priors = list(
    tau_prior = 1/sqrt(ncol(data_hd))  # Custom global shrinkage
  ),
  ps.priors = list(
    tau_prior = 1/sqrt(ncol(data_hd))
  )
)
```

## Documentation

### Vignettes

DRBayes includes comprehensive vignettes:

``` r

# Basic usage tutorial
vignette("getting-started", package = "DRBayes")

# External software integration guide
vignette("external-software-integration", package = "DRBayes")

# Browse all vignettes
browseVignettes("DRBayes")
```

### Help Documentation

``` r

# Main function help
?drbayes_pc

# Model functions
?bayes_lm
?bayes_logit
?bayes_lm_hs

# Package overview
help(package = "DRBayes")
```

## Methodology

The posterior coupling approach implemented in DRBayes addresses key
limitations of traditional doubly robust methods:

1.  **Separate Model Fitting**: Outcome and propensity score models are
    fit independently using the formula interface
2.  **Moment Condition Enforcement**: The joint posterior is tilted by
    `exp(lambda * B_n)` until the posterior mean of the moment condition
    is zero. `method = "smc"` (the default) walks a grid of tilting
    parameters and rejuvenates the particles at each step, which is
    Algorithm 2 of the paper; `method = "is"` takes a single importance
    sampling step, which is Algorithm 1
3.  **Coupled Posterior**: Results in a posterior distribution that
    satisfies the doubly robust property

The method enforces the moment condition:

\\E\left\[\frac{(A - \pi(X)) \cdot (Y -
m_A(X))}{\pi(X)(1-\pi(X))}\right\] = 0\\

where \\\pi(X)\\ is the propensity score and \\m_A(X)\\ is the outcome
model evaluated at the observed treatment.

**For binary outcomes**: \\m_A(X)\\ represents the predicted probability
using the specified link function, and the ATE represents the average
probability difference between treatment and control groups.

## Citation

If you use DRBayes in your research, please cite:

``` bibtex
@misc{orihara2025bayesian,
  title={Bayesian Doubly Robust Causal Inference via Posterior Coupling},
  author={Orihara, Shunichiro and Momozaki, Tomotaka and Sugasawa, Shonosuke},
  year={2025},
  eprint={2506.04868},
  archivePrefix={arXiv},
  primaryClass={stat.ME}
}
```

**Reference**: Orihara, S., Momozaki, T., & Sugasawa, S. (2025).
Bayesian Doubly Robust Causal Inference via Posterior Coupling. *arXiv
preprint arXiv:2506.04868*.

## Examples and Use Cases

### Observational Studies

- **Healthcare outcomes research**: Both continuous (blood pressure
  reduction) and binary (treatment success) outcomes
- **Educational intervention evaluation**: Test scores (continuous) or
  graduation rates (binary)
- **Economic policy analysis**: Income effects (continuous) or
  employment status (binary)
- **Social science applications**: Attitude scales (continuous) or
  behavioral choices (binary)

## Troubleshooting

### Common Issues

1.  **Convergence problems**: Increase `mc` or adjust priors
2.  **Extreme propensity scores**: Check for positivity violations
3.  **High-dimensional data**: Use `outcome.prior = "horseshoe"` and
    `ps.prior = "horseshoe"`
4.  **Missing data**: Use `na.action = "na.omit"` or preprocess data

### Diagnostic Tools

``` r

# Both estimands, their credible intervals, the tilting parameter, whether
# the moment condition was met, and the worst R-hat and ESS of the two models
print(result)
summary(result)

# The two posterior densities
plot(result)

# The full convergence table, one row per coefficient of each model
result$diagnostics

# Data diagnostics
print(result$data_info)
```

R-hat and effective sample size are reported for the Markov chain draws
of the outcome and propensity score models, which is where they apply.
They are not reported for `result$pc`: those are sequential Monte Carlo
particles, which are resampled and carry no time ordering, so a trace
plot or an autocorrelation based effective sample size computed on them
would be estimating nothing.

## Reporting Issues

If you encounter bugs or have feature requests, please file an issue on
our [GitHub Issues page](https://github.com/t-momozaki/DRBayes/issues).

## Related Packages

- [`AIPW`](https://CRAN.R-project.org/package=AIPW): Augmented inverse
  probability weighting
- [`tmle`](https://CRAN.R-project.org/package=tmle): Targeted maximum
  likelihood estimation
- [`bartCause`](https://CRAN.R-project.org/package=bartCause): Bayesian
  causal inference with BART
- [`MatchIt`](https://CRAN.R-project.org/package=MatchIt): Matching
  methods for causal inference

## Authors

- **Tomotaka Momozaki** - *Package Developer* - Tokyo University of
  Science
- **Shunichiro Orihara** - *Methodology Development* - Tokyo Medical
  University
- **Shonosuke Sugasawa** - *Methodology Development* - Keio University

## License

This project is licensed under the MIT License - see the
[LICENSE](https://t-momozaki.github.io/DRBayes/LICENSE.md) file for
details.

## Acknowledgments

- Thanks to the R community for providing excellent tools for Bayesian
  computation
- Special thanks to the Stan, brms, and JAGS development teams for
  creating compatible software

------------------------------------------------------------------------

For more information, visit our [GitHub
repository](https://github.com/t-momozaki/DRBayes) or contact the
maintainer at <t_momozaki@rs.tus.ac.jp>.

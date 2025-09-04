---
title: 'DRBayes: Bayesian Doubly Robust Causal Inference via Posterior Coupling'
---

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%E2%89%A5%203.5.0-blue.svg)](https://www.r-project.org/)
[![Development](https://img.shields.io/badge/Status-Development-orange.svg)](https://github.com/t-momozaki/DRBayes)
<!-- badges: end -->

## Overview

DRBayes implements Bayesian doubly robust methods for causal inference using the novel posterior coupling approach. This methodology combines the robustness of doubly robust estimation with the uncertainty quantification benefits of Bayesian inference, providing a principled framework for estimating average treatment effects.

### Key Features

- **Formula Interface**: Intuitive R formula syntax compatible with `lm()`, `glm()`, and other standard functions
- **Doubly Robust Estimation**: Consistent estimation even when either the outcome model or propensity score model is misspecified
- **Posterior Coupling**: Novel methodology that avoids the feedback problem in traditional doubly robust methods
- **Bayesian Uncertainty Quantification**: Full posterior distributions for treatment effects
- **Multiple Outcome Types**: Support for continuous (linear regression) and binary outcomes (logistic/probit regression) 
- **Automatic Detection**: Smart detection of outcome family and link functions based on data
- **Flexible Modeling**: Built-in functions with standard and horseshoe priors plus support for external Bayesian software
- **Missing Data Handling**: Comprehensive support for missing value treatment
- **Software Integration**: Seamless integration with Stan, brms, JAGS, and other Bayesian platforms

## Installation

You can install the development version of DRBayes from GitHub:

```r
# Install devtools if you haven't already
if (!require(devtools)) install.packages("devtools")

# Install DRBayes with vignettes
devtools::install_github("t-momozaki/DRBayes", build_vignettes = TRUE)
```

### Dependencies

DRBayes requires the following R packages:

**Core dependencies:**
- `MCMCpack` (Dirichlet random numbers)
- `mvtnorm` (Multivariate normal sampling)
- `stats`, `base` (Standard R functions)

**Model-specific dependencies:**
- `extraDistr` (Inverse gamma distribution)
- `RcppTN` (Truncated normal sampling for probit models)
- `pgdraw` (Pólya-Gamma sampling for logistic models)

**Optional for enhanced functionality:**
- `future`, `furrr`, `parallelly` (Parallel processing for simulations)
- `nleqslv` (Equation solving)
- `cmdstanr` (Stan integration)
- `brms` (High-level Bayesian modeling)
- `R2jags` (JAGS integration)

## Quick Start

### Example 1: Continuous Outcome with Formula Interface

```r
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
result_continuous <- DRBayes.PC(
  outcome.formula = Y ~ A + X1 + X2 + A:X1,      # Outcome model with interaction
  ps.formula = A ~ X1 + X2 + X3,                 # Propensity score model
  data = data,
  family = "gaussian",                           # Auto-detected if NULL
  outcome.model = B.LM,                          # Bayesian linear regression
  ps.model = B.Logit,                            # Bayesian logistic regression
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

```r
# Generate binary outcome
data$Y_binary <- rbinom(n, 1, plogis(0.5 + 0.8*data$A + 0.4*data$X1 + 0.3*data$X2))

# Fit model with automatic family/link detection
result_binary <- DRBayes.PC(
  outcome.formula = Y_binary ~ A + X1 + X2,
  ps.formula = A ~ X1 + X2 + X3,
  data = data,
  # family and link auto-detected as "binomial" and "logit"
  outcome.model = B.Logit,
  ps.model = B.Logit,
  mc = 3000, bn = 1000, thin = 2
)

# Results represent probability differences for binary outcomes
cat("Family used:", result_binary$family, "\n")
cat("Link function:", result_binary$link, "\n")
cat("ATE (probability difference):", round(mean(result_binary$pc), 3), "\n")
```

### Example 3: Using Horseshoe Priors for High-Dimensional Data

```r
# Generate high-dimensional data
set.seed(1234)
n <- 200
p <- 50  # Many covariates
data_hd <- data.frame(matrix(rnorm(n * p), n, p))
names(data_hd) <- paste0("X", 1:p)

# Only first few covariates are truly relevant
data_hd$A <- rbinom(n, 1, plogis(0.3*data_hd$X1 + 0.2*data_hd$X2 - 0.1*data_hd$X3))
data_hd$Y <- 2 + 1.5*data_hd$A + 0.8*data_hd$X1 + 0.5*data_hd$X2 + rnorm(n)

# Use horseshoe priors for regularization
result_horseshoe <- DRBayes.PC(
  outcome.formula = Y ~ A + .,  # Include all covariates
  ps.formula = A ~ . - Y,       # All covariates except outcome
  data = data_hd,
  outcome.model = HS.LM,        # Horseshoe linear regression
  ps.model = HS.Logit,          # Horseshoe logistic regression
  mc = 2000, bn = 500, thin = 2
)

cat("Horseshoe ATE:", round(mean(result_horseshoe$pc), 3), "\n")
```

### Example 4: Advanced Formula Usage and Missing Data

```r
# Create data with missing values and complex relationships
data_missing <- data
data_missing$X1[sample(n, 20)] <- NA
data_missing$X2[sample(n, 15)] <- NA

# Complex formulas with transformations and polynomials
result_advanced <- DRBayes.PC(
  outcome.formula = Y ~ A + log(abs(X1) + 1) + I(X2^2) + poly(X3, 2) + A:X1,
  ps.formula = A ~ X1 + X2 + X3 + I(X1^2) + X2:X3,
  data = data_missing,
  na.action = "na.omit",                         # Handle missing values
  outcome.model = B.LM,
  ps.model = B.Logit,
  mc = 2000, bn = 500
)

cat("Original sample size:", n, "\n")
cat("Effective sample size:", result_advanced$data_info$n_observations, "\n")
cat("Missing observations:", result_advanced$data_info$missing_observations, "\n")
cat("Advanced formula ATE:", round(mean(result_advanced$pc), 3), "\n")
```

### Example 5: Bayesian Bootstrap Alternative

```r
# Alternative: Bayesian Bootstrap approach (Saarela et al. 2016)
X.lm <- data.frame(
  A = data$A,
  X1 = data$X1, X2 = data$X2, X3 = data$X3,
  A_X1 = data$A * data$X1  # Treatment-covariate interaction
)
X.ps <- data.frame(X1 = data$X1, X2 = data$X2, X3 = data$X3)

result_bb <- DRBayes.BB(
  Y = data$Y,
  A = data$A,
  X.lm = X.lm,
  X.ps = X.ps,
  num_iterations = 1000,
  family = "gaussian"
)

cat("Bayesian Bootstrap ATE:", round(mean(result_bb), 3), "\n")
```

## Main Functions

### Core Functions

- **`DRBayes.PC()`**: Main function for Bayesian doubly robust estimation via posterior coupling with formula interface
- **`DRBayes.BB()`**: Bayesian bootstrap approach (Saarela et al. 2016)

### Built-in Model Functions

**Standard Priors:**
- **`B.LM()`**: Bayesian linear regression with conjugate priors (for continuous outcomes)
- **`B.Logit()`**: Bayesian logistic regression using Pólya-Gamma data augmentation (for binary outcomes/treatment)
- **`B.Probit()`**: Bayesian probit regression using latent variable approach (for binary outcomes/treatment)

**Horseshoe Priors (for high-dimensional data):**
- **`HS.LM()`**: Horseshoe prior linear regression
- **`HS.Logit()`**: Horseshoe prior logistic regression
- **`HS.Probit()`**: Horseshoe prior probit regression

### Simulation and Utility Functions

- **`generate_dataset()`**: Generate synthetic datasets for simulation studies
- **`run_parallel_simulation()`**: Run parallel simulation studies

## Outcome Types and Model Families 

DRBayes supports multiple outcome types with automatic detection:

| Family | Link Function | Outcome Type | Use Case | Auto-Detection |
|--------|---------------|--------------|----------|----------------|
| `"gaussian"` | `"identity"` | Continuous | Linear regression, mean differences | Non-binary Y values |
| `"binomial"` | `"logit"` | Binary (0/1) | Logistic regression, probability differences | All Y ∈ {0,1} |
| `"binomial"` | `"probit"` | Binary (0/1) | Probit regression, probability differences | Manual specification |

**Automatic Detection**: When `family = NULL` (default), the function automatically detects the appropriate family:
- If all Y values are 0 or 1 → `family = "binomial"`, `link = "logit"`
- Otherwise → `family = "gaussian"`, `link = "identity"`

## Formula Interface Features

The formula interface provides several advantages:

### Intuitive Syntax
```r
# Similar to lm() and glm() - familiar to R users
DRBayes.PC(
  outcome.formula = Y ~ A + X1 + X2 + A:X1,
  ps.formula = A ~ X1 + X2 + X3,
  data = data
)
```

### Automatic Handling
- **Design matrices**: Automatically created from formulas
- **Factor variables**: Properly handled with contrast coding
- **Interactions**: Simple syntax for interaction terms (`A:X1`, `A*X1`)
- **Transformations**: Support for `I()`, `log()`, `poly()`, etc.
- **Missing values**: Consistent handling across models

### Complex Modeling
```r
# Advanced transformations and interactions
outcome.formula = Y ~ A + log(X1) + I(X2^2) + poly(X3, 3) + A:X1 + A:I(X2^2)
ps.formula = A ~ X1 + X2 + X3 + I(X1^2) + X1:X2
```

## External Software Integration

DRBayes supports posterior samples from popular Bayesian software with the formula interface:

### Stan/cmdstanr Integration

```r
# Extract samples from Stan fit
outcome_samples <- fit_stan$draws("beta", format = "matrix")
ps_samples <- fit_stan$draws("gamma", format = "matrix")

# Use with formula interface
result <- DRBayes.PC(
  outcome.formula = Y ~ A + X1 + X2 + A:X1,
  ps.formula = A ~ X1 + X2 + X3,
  data = data,
  outcome.samples = outcome_samples,
  ps.samples = ps_samples,
  bn = 100, thin = 2
)
```

### brms Integration

```r
# Extract fixed effects from brms
outcome_samples <- as.matrix(brms_fit)[, 1:ncol(model.matrix(outcome.formula, data))]
ps_samples <- as.matrix(brms_ps_fit)[, 1:ncol(model.matrix(ps.formula, data))]

result <- DRBayes.PC(
  outcome.formula = Y ~ A + X1 + X2,
  ps.formula = A ~ X1 + X2 + X3,
  data = data,
  outcome.samples = outcome_samples,
  ps.samples = ps_samples
)
```

### JAGS/R2jags Integration

```r
# Extract parameter samples from JAGS
outcome_samples <- jags_fit$BUGSoutput$sims.matrix[, grep("beta", colnames(...))]
ps_samples <- jags_fit$BUGSoutput$sims.matrix[, grep("gamma", colnames(...))]

result <- DRBayes.PC(
  outcome.formula = Y ~ A + X1 + X2,
  ps.formula = A ~ X1 + X2 + X3,
  data = data,
  outcome.samples = outcome_samples,
  ps.samples = ps_samples
)
```

## Prior Specification

Customize priors through the built-in model functions:

```r
# Custom priors for internal sampling
result <- DRBayes.PC(
  outcome.formula = Y ~ A + X1 + X2,
  ps.formula = A ~ X1 + X2 + X3,
  data = data,
  outcome.model = B.LM,
  ps.model = B.Logit,
  outcome.priors = list(
    theta_prior = 1/50,        # Prior precision for coefficients
    sigma_prior = c(2, 1)      # Inverse gamma prior for error variance
  ),
  ps.priors = list(
    theta_prior = 1/25         # Prior precision for PS coefficients
  )
)

# Horseshoe priors with custom tau
result_hs <- DRBayes.PC(
  outcome.formula = Y ~ A + .,
  ps.formula = A ~ . - Y,
  data = data_hd,
  outcome.model = HS.LM,
  ps.model = HS.Logit,
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

```r
# Basic usage tutorial
vignette("getting-started", package = "DRBayes")

# External software integration guide
vignette("external-software-integration", package = "DRBayes")

# Browse all vignettes
browseVignettes("DRBayes")
```

### Help Documentation

```r
# Main function help
?DRBayes.PC

# Model functions
?B.LM
?B.Logit
?HS.LM

# Package overview
help(package = "DRBayes")
```

## Methodology

The posterior coupling approach implemented in DRBayes addresses key limitations of traditional doubly robust methods:

1. **Separate Model Fitting**: Outcome and propensity score models are fit independently using the formula interface
2. **Moment Condition Enforcement**: Sequential Monte Carlo enforces the doubly robust moment condition
3. **Coupled Posterior**: Results in a posterior distribution that satisfies the doubly robust property

The method enforces the moment condition:

$$E\left[\frac{(A - \pi(X)) \cdot (Y - \mu(X))}{\pi(X)(1-\pi(X))}\right] = 0$$

where $\pi(X)$ is the propensity score and $\mu(X)$ is the outcome model.

**For binary outcomes**: $\mu(X)$ represents the predicted probability using the specified link function, and the ATE represents the average probability difference between treatment and control groups.

## Citation

If you use DRBayes in your research, please cite:

```bibtex
@misc{orihara2025bayesian,
  title={Bayesian Doubly Robust Causal Inference via Posterior Coupling},
  author={Orihara, Shunichiro and Momozaki, Tomotaka and Sugasawa, Shonosuke},
  year={2025},
  eprint={2506.04868},
  archivePrefix={arXiv},
  primaryClass={stat.ME}
}
```

**Reference**: Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust Causal Inference via Posterior Coupling. *arXiv preprint arXiv:2506.04868*.

## Examples and Use Cases

### Observational Studies

- **Healthcare outcomes research**: Both continuous (blood pressure reduction) and binary (treatment success) outcomes
- **Educational intervention evaluation**: Test scores (continuous) or graduation rates (binary)
- **Economic policy analysis**: Income effects (continuous) or employment status (binary)
- **Social science applications**: Attitude scales (continuous) or behavioral choices (binary)

## Troubleshooting

### Common Issues

1. **Convergence problems**: Increase `mc` or adjust priors
2. **Extreme propensity scores**: Check for positivity violations
3. **High-dimensional data**: Use horseshoe priors (`HS.*` functions)
4. **Missing data**: Use `na.action = "na.omit"` or preprocess data

### Diagnostic Tools

```r
# Check posterior samples
plot(result$pc, type = "l")  # Trace plot
hist(result$pc)              # Posterior distribution

# Compare methods
cat("PC ATE:", mean(result$pc), "\n")
cat("G-comp ATE:", mean(result$g.comp), "\n")

# Data diagnostics
print(result$data_info)
```

## Reporting Issues

If you encounter bugs or have feature requests, please file an issue on our [GitHub Issues page](https://github.com/t-momozaki/DRBayes/issues).

## Related Packages

- [`AIPW`](https://cran.r-project.org/package=AIPW): Augmented inverse probability weighting
- [`tmle`](https://cran.r-project.org/package=tmle): Targeted maximum likelihood estimation
- [`CausalInference`](https://cran.r-project.org/package=CausalInference): Various causal inference methods
- [`MatchIt`](https://cran.r-project.org/package=MatchIt): Matching methods for causal inference

## Authors

- **Tomotaka Momozaki** - *Package Developer* - Tokyo University of Science
- **Shunichiro Orihara** - *Methodology Development* - Tokyo Medical University
- **Shonosuke Sugasawa** - *Methodology Development* - Keio University

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Thanks to the R community for providing excellent tools for Bayesian computation
- Special thanks to the Stan, brms, and JAGS development teams for creating compatible software

---

For more information, visit our [GitHub repository](https://github.com/t-momozaki/DRBayes) or contact the maintainer at t_momozaki@rs.tus.ac.jp.

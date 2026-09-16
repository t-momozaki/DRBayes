# Integrating DRBayes with External Bayesian Software

## Introduction

The **DRBayes** package implements Bayesian doubly robust causal
inference methods via posterior coupling. One of the key features of
this package is its flexibility in working with posterior samples from
various Bayesian software packages, including Stan, JAGS, and brms.

This vignette demonstrates how to integrate
[`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
with external Bayesian software for both **continuous and binary
outcomes**, allowing researchers to leverage their preferred modeling
frameworks while benefiting from the doubly robust posterior coupling
methodology.

### Key Features

- **Formula Interface**: Intuitive R formula syntax compatible with
  [`lm()`](https://rdrr.io/r/stats/lm.html),
  [`glm()`](https://rdrr.io/r/stats/glm.html), and other standard
  functions
- **Dual operation modes**: Internal sampling with built-in functions or
  external posterior samples
- **Multiple outcome types**: Support for both continuous and binary
  outcomes
- **Flexible model families**: Works with Gaussian (linear) and binomial
  (logistic/probit) models
- **Automatic detection**: Smart detection of outcome family and link
  functions based on data
- **Wide compatibility**: Works with Stan/cmdstanr, brms, JAGS, and
  other Bayesian software
- **Seamless integration**: Uses standard formula interface for
  consistent syntax

## Installation and Setup

``` r
# Install from GitHub (when available)
# devtools::install_github("t-momozaki/DRBayes")

# For Stan integration
# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
```

``` r
# Load the DRBayes package
library(DRBayes)

# Load additional packages for external software integration
# library(cmdstanr)  # For Stan integration
# ggplot2 and dplyr are optional (Suggests); they are attached only if present
if (has_ggplot2) library(ggplot2)   # For visualization
if (has_dplyr) library(dplyr)       # For data manipulation and pipe operator (%>%)
```

## Data Setup

Let’s generate both continuous and binary outcome examples that we’ll
use throughout this vignette:

``` r
# Set seed for reproducibility
set.seed(1234)

# Generate synthetic data
n <- 500  # Sample size

# Create a comprehensive data frame
data <- data.frame(
  X1 = rnorm(n),
  X2 = rnorm(n),
  X3 = rnorm(n),
  X4 = rbinom(n, 1, 0.4)  # Binary covariate
)

# True treatment assignment mechanism
ps_logits <- 0.5 + 0.4*data$X1 + 0.3*data$X2 - 0.2*data$X3 + 0.5*data$X4
data$A <- rbinom(n, 1, plogis(ps_logits))

# Continuous outcome (ATE = 1.5)
data$Y_continuous <- 2 + 1.5*data$A + 0.8*data$X1 + 0.6*data$X2 + 0.4*data$X3 + 0.3*data$A*data$X1 + rnorm(n, 0, 1)

# Binary outcome; the true ATE in probability difference is computed below
linear_pred <- -0.5 + 0.8*data$A + 0.4*data$X1 + 0.3*data$X2 + 0.2*data$X3
data$Y_binary <- rbinom(n, 1, plogis(linear_pred))

# Calculate true ATE for binary outcome using counterfactual approach
prob_treated <- plogis(-0.5 + 0.8*1 + 0.4*data$X1 + 0.3*data$X2 + 0.2*data$X3)
prob_control <- plogis(-0.5 + 0.8*0 + 0.4*data$X1 + 0.3*data$X2 + 0.2*data$X3)
true_ate_binary <- mean(prob_treated - prob_control)

cat("Dataset created:\n")
#> Dataset created:
cat("Sample size:", n, "\n")
#> Sample size: 500
cat("Treatment prevalence:", round(mean(data$A), 3), "\n")
#> Treatment prevalence: 0.618
cat("True ATE (continuous) =", 1.5, "\n")
#> True ATE (continuous) = 1.5
cat("True ATE (binary, prob diff) =", round(true_ate_binary, 3), "\n")
#> True ATE (binary, prob diff) = 0.184

# Display data structure
cat("\nData structure:\n")
#> 
#> Data structure:
str(data)
#> 'data.frame':    500 obs. of  7 variables:
#>  $ X1          : num  -1.207 0.277 1.084 -2.346 0.429 ...
#>  $ X2          : num  0.985 -1.225 0.71 -0.109 1.783 ...
#>  $ X3          : num  -1.205 0.301 -1.539 0.635 0.703 ...
#>  $ X4          : int  1 0 0 0 0 0 0 1 0 1 ...
#>  $ A           : int  0 1 1 1 1 1 0 1 0 1 ...
#>  $ Y_continuous: num  0.169 3.091 4.392 2.301 3.667 ...
#>  $ Y_binary    : int  0 1 0 1 1 0 0 1 0 0 ...
```

## Basic Usage: Internal Sampling

Before diving into external software integration, let’s see how
[`drbayes_pc()`](https://t-momozaki.github.io/DRBayes/reference/drbayes_pc.md)
works with its built-in model functions for both outcome types:

``` r
# Continuous outcome using built-in models
result_internal_continuous <- drbayes_pc(
  outcome.formula = Y_continuous ~ A + X1 + X2 + X3 + A:X1,  # With interaction
  ps.formula = A ~ X1 + X2 + X3 + X4,                        # Include more predictors
  data = data,
  family = "gaussian",              # For continuous outcomes
  outcome.model = bayes_lm,             # Bayesian linear regression
  ps.model = bayes_logit,               # Bayesian logistic regression
  mc = 1200, bn = 300, thin = 2,
  outcome.priors = list(theta_prior = 1/100, sigma_prior = c(1, 1)),
  ps.priors = list(theta_prior = 1/100)
)

# Binary outcome using built-in models
result_internal_binary <- drbayes_pc(
  outcome.formula = Y_binary ~ A + X1 + X2 + X3,
  ps.formula = A ~ X1 + X2 + X3 + X4,
  data = data,
  # family automatically detected as "binomial" since Y_binary is 0/1
  outcome.model = bayes_logit,          # Bayesian logistic regression
  ps.model = bayes_logit,               # Bayesian logistic regression
  mc = 1200, bn = 300, thin = 2,
  outcome.priors = list(theta_prior = 1/100),
  ps.priors = list(theta_prior = 1/100)
)

# Summarize results
cat("=== INTERNAL SAMPLING RESULTS ===\n")
#> === INTERNAL SAMPLING RESULTS ===
cat("Continuous Outcome:\n")
#> Continuous Outcome:
cat("  Family used:", result_internal_continuous$family, "with", result_internal_continuous$link, "link\n")
#>   Family used: gaussian with identity link
cat("  ATE (PC): Mean =", round(mean(result_internal_continuous$pc), 3), 
    ", SD =", round(sd(result_internal_continuous$pc), 3), "\n")
#>   ATE (PC): Mean = 1.563 , SD = 0.023

cat("Binary Outcome:\n")
#> Binary Outcome:
cat("  Family used:", result_internal_binary$family, "with", result_internal_binary$link, "link\n")
#>   Family used: binomial with logit link
cat("  ATE (PC): Mean =", round(mean(result_internal_binary$pc), 3), 
    ", SD =", round(sd(result_internal_binary$pc), 3), "\n")
#>   ATE (PC): Mean = 0.144 , SD = 0.045
```

## External Software Integration

### Understanding Design Matrix Dimensions

When using external samples, it’s crucial to understand how the formula
interface creates design matrices:

``` r
# Check design matrix dimensions for our formulas
outcome_matrix_continuous <- model.matrix(Y_continuous ~ A + X1 + X2 + X3 + A:X1, data)
ps_matrix <- model.matrix(A ~ X1 + X2 + X3 + X4, data)
outcome_matrix_binary <- model.matrix(Y_binary ~ A + X1 + X2 + X3, data)

cat("Design Matrix Dimensions:\n")
#> Design Matrix Dimensions:
cat("Continuous outcome model:", ncol(outcome_matrix_continuous), "parameters\n")
#> Continuous outcome model: 6 parameters
cat("Binary outcome model:", ncol(outcome_matrix_binary), "parameters\n")
#> Binary outcome model: 5 parameters
cat("Propensity score model:", ncol(ps_matrix), "parameters\n")
#> Propensity score model: 5 parameters

# Show parameter names
cat("\nParameter names:\n")
#> 
#> Parameter names:
cat("Continuous outcome:", colnames(outcome_matrix_continuous), "\n")
#> Continuous outcome: (Intercept) A X1 X2 X3 A:X1
cat("Binary outcome:", colnames(outcome_matrix_binary), "\n")
#> Binary outcome: (Intercept) A X1 X2 X3
cat("PS model:", colnames(ps_matrix), "\n")
#> PS model: (Intercept) X1 X2 X3 X4
```

### Stan/cmdstanr Integration

Stan is a popular platform for statistical modeling. Here’s how to
integrate Stan models with DRBayes for both continuous and binary
outcomes:

#### Step 1: Define Stan Models

``` r
# Stan model for linear regression (continuous outcome)
stan_lm <- "
data {
  int<lower=0> N;              // number of observations
  int<lower=0> K;              // number of parameters
  matrix[N, K] X;              // design matrix (including intercept)
  vector[N] y;                 // response variable
}
parameters {
  vector[K] beta;              // regression coefficients
  real<lower=0> sigma;         // error standard deviation
}
model {
  // Priors
  beta ~ normal(0, sqrt(100)); // weakly informative priors
  sigma ~ inv_gamma(1, 1);     // prior for error variance
  
  // Likelihood
  y ~ normal(X * beta, sigma);
}
"

# Stan model for logistic regression (binary outcome & propensity score)
stan_logit <- "
data {
  int<lower=0> N;                     // number of observations
  int<lower=0> K;                     // number of parameters
  matrix[N, K] X;                     // design matrix (including intercept)
  array[N] int<lower=0, upper=1> y;   // binary response
}
parameters {
  vector[K] beta;              // regression coefficients
}
model {
  // Priors
  beta ~ normal(0, sqrt(100)); // weakly informative priors
  
  // Likelihood
  y ~ bernoulli_logit(X * beta);
}
"

# Stan model for probit regression (alternative for binary outcomes)
stan_probit <- "
data {
  int<lower=0> N;                     // number of observations
  int<lower=0> K;                     // number of parameters
  matrix[N, K] X;                     // design matrix (including intercept)
  array[N] int<lower=0, upper=1> y;   // binary response
}
parameters {
  vector[K] beta;              // regression coefficients
}
model {
  // Priors
  beta ~ normal(0, sqrt(100)); // weakly informative priors
  
  // Likelihood
  y ~ bernoulli(Phi(X * beta)); // Probit model
}
"
```

#### Step 2: Compile and Fit Stan Models

``` r
# Compile models
mod_lm     <- cmdstanr::cmdstan_model(stan_file = cmdstanr::write_stan_file(stan_lm))
mod_logit  <- cmdstanr::cmdstan_model(stan_file = cmdstanr::write_stan_file(stan_logit))
mod_probit <- cmdstanr::cmdstan_model(stan_file = cmdstanr::write_stan_file(stan_probit))

cat("Fitting Stan models...\n")
#> Fitting Stan models...

# === CONTINUOUS OUTCOME MODEL ===
stan_data_continuous <- list(
  N = nrow(data),
  K = ncol(outcome_matrix_continuous),
  X = outcome_matrix_continuous,
  y = data$Y_continuous
)

fit_continuous_stan <- mod_lm$sample(
  data = stan_data_continuous,
  chains = 2,
  parallel_chains = 2,
  iter_sampling = 1000,
  iter_warmup = 1000,
  refresh = 0,
  show_messages = FALSE, 
  show_exceptions = FALSE
)

# === BINARY OUTCOME MODELS ===
stan_data_binary <- list(
  N = nrow(data),
  K = ncol(outcome_matrix_binary),
  X = outcome_matrix_binary,
  y = data$Y_binary
)

# Logistic regression
fit_binary_logit_stan <- mod_logit$sample(
  data = stan_data_binary,
  chains = 2,
  parallel_chains = 2,
  iter_sampling = 1000,
  iter_warmup = 1000,
  refresh = 0,
  show_messages = FALSE, 
  show_exceptions = FALSE
)

# Probit regression
fit_binary_probit_stan <- mod_probit$sample(
  data = stan_data_binary,
  chains = 2,
  parallel_chains = 2,
  iter_sampling = 1000,
  iter_warmup = 1000,
  refresh = 0,
  show_messages = FALSE, 
  show_exceptions = FALSE
)

# === PROPENSITY SCORE MODEL ===
stan_data_ps <- list(
  N = nrow(data),
  K = ncol(ps_matrix),
  X = ps_matrix,
  y = data$A
)

fit_ps_stan <- mod_logit$sample(
  data = stan_data_ps,
  chains = 2,
  parallel_chains = 2,
  iter_sampling = 1000,
  iter_warmup = 1000,
  refresh = 0,
  show_messages = FALSE, 
  show_exceptions = FALSE
)

cat("Stan model fitting completed!\n")
#> Stan model fitting completed!
```

#### Step 3: Extract Samples and Integrate with DRBayes

``` r
# Extract samples from Stan fits
ps_samples_stan <- fit_ps_stan$draws("beta", format = "matrix")
outcome_samples_continuous <- fit_continuous_stan$draws("beta", format = "matrix")
outcome_samples_binary_logit <- fit_binary_logit_stan$draws("beta", format = "matrix")
outcome_samples_binary_probit <- fit_binary_probit_stan$draws("beta", format = "matrix")

# Verify sample dimensions
cat("Sample dimensions:\n")
#> Sample dimensions:
cat("PS samples:", dim(ps_samples_stan), "\n")
#> PS samples: 2000 5
cat("Continuous outcome samples:", dim(outcome_samples_continuous), "\n")
#> Continuous outcome samples: 2000 6
cat("Binary logit samples:", dim(outcome_samples_binary_logit), "\n")
#> Binary logit samples: 2000 5
cat("Binary probit samples:", dim(outcome_samples_binary_probit), "\n")
#> Binary probit samples: 2000 5

# === CONTINUOUS OUTCOME INTEGRATION ===
result_stan_continuous <- drbayes_pc(
  outcome.formula = Y_continuous ~ A + X1 + X2 + X3 + A:X1,
  ps.formula = A ~ X1 + X2 + X3 + X4,
  data = data,
  family = "gaussian",                    # Must specify for external samples
  outcome.samples = outcome_samples_continuous,
  ps.samples = ps_samples_stan,
  bn = 0, thin = 2  
)

# === BINARY OUTCOME INTEGRATION (LOGIT) ===
result_stan_binary_logit <- drbayes_pc(
  outcome.formula = Y_binary ~ A + X1 + X2 + X3,
  ps.formula = A ~ X1 + X2 + X3 + X4,
  data = data,
  family = "binomial", 
  link = "logit",
  outcome.samples = outcome_samples_binary_logit,
  ps.samples = ps_samples_stan,
  bn = 0, thin = 2
)

# === BINARY OUTCOME INTEGRATION (PROBIT) ===
result_stan_binary_probit <- drbayes_pc(
  outcome.formula = Y_binary ~ A + X1 + X2 + X3,
  ps.formula = A ~ X1 + X2 + X3 + X4,
  data = data,
  family = "binomial", 
  link = "probit",
  outcome.samples = outcome_samples_binary_probit,
  ps.samples = ps_samples_stan,
  bn = 0, thin = 2
)

# Summarize Stan integration results
cat("\n=== STAN INTEGRATION RESULTS ===\n")
#> 
#> === STAN INTEGRATION RESULTS ===
cat("Continuous (Linear):\n")
#> Continuous (Linear):
cat("  ATE (PC): Mean =", round(mean(result_stan_continuous$pc), 3), 
    ", SD =", round(sd(result_stan_continuous$pc), 3), "\n")
#>   ATE (PC): Mean = 1.563 , SD = 0.097

cat("Binary (Logit):\n")
#> Binary (Logit):
cat("  ATE (PC): Mean =", round(mean(result_stan_binary_logit$pc), 3), 
    ", SD =", round(sd(result_stan_binary_logit$pc), 3), "\n")
#>   ATE (PC): Mean = 0.147 , SD = 0.044

cat("Binary (Probit):\n")
#> Binary (Probit):
cat("  ATE (PC): Mean =", round(mean(result_stan_binary_probit$pc), 3), 
    ", SD =", round(sd(result_stan_binary_probit$pc), 3), "\n")
#>   ATE (PC): Mean = 0.145 , SD = 0.047
```

## Comprehensive Comparison

``` r
# Create a comprehensive comparison function
compare_all_results <- function() {
  
  results_list <- list(
    "Built-in Continuous" = result_internal_continuous,
    "Built-in Binary" = result_internal_binary,
    "Stan Continuous" = result_stan_continuous,
    "Stan Binary (Logit)" = result_stan_binary_logit,
    "Stan Binary (Probit)" = result_stan_binary_probit
  )
  
  outcome_types <- c("Continuous", "Binary", "Continuous", "Binary", "Binary")
  true_ates <- c(1.5, true_ate_binary, 1.5, true_ate_binary, true_ate_binary)
  
  comparison_df <- data.frame(
    Method = names(results_list),
    Outcome_Type = outcome_types,
    True_ATE = true_ates,
    PC_Mean = sapply(results_list, function(x) mean(x$pc)),
    PC_SD = sapply(results_list, function(x) sd(x$pc)),
    PC_CI_Lower = sapply(results_list, function(x) quantile(x$pc, 0.025)),
    PC_CI_Upper = sapply(results_list, function(x) quantile(x$pc, 0.975)),
    Family = sapply(results_list, function(x) x$family),
    Link = sapply(results_list, function(x) x$link)
  )
  
  # Calculate bias and coverage
  comparison_df$Bias <- comparison_df$PC_Mean - comparison_df$True_ATE
  comparison_df$Coverage <- mapply(function(true_val, ci_low, ci_high) {
    (true_val >= ci_low) & (true_val <= ci_high)
  }, comparison_df$True_ATE, comparison_df$PC_CI_Lower, comparison_df$PC_CI_Upper)
  
  # Round numeric columns
  numeric_cols <- c("True_ATE", "PC_Mean", "PC_SD", "PC_CI_Lower", "PC_CI_Upper", "Bias")
  comparison_df[, numeric_cols] <- round(comparison_df[, numeric_cols], 3)
  
  return(comparison_df)
}

comparison_table <- compare_all_results()
print(comparison_table)
#>                                    Method Outcome_Type True_ATE PC_Mean PC_SD
#> Built-in Continuous   Built-in Continuous   Continuous    1.500   1.563 0.023
#> Built-in Binary           Built-in Binary       Binary    0.184   0.144 0.045
#> Stan Continuous           Stan Continuous   Continuous    1.500   1.563 0.097
#> Stan Binary (Logit)   Stan Binary (Logit)       Binary    0.184   0.147 0.044
#> Stan Binary (Probit) Stan Binary (Probit)       Binary    0.184   0.145 0.047
#>                      PC_CI_Lower PC_CI_Upper   Family     Link   Bias Coverage
#> Built-in Continuous        1.520       1.606 gaussian identity  0.063    FALSE
#> Built-in Binary            0.055       0.230 binomial    logit -0.039     TRUE
#> Stan Continuous            1.355       1.726 gaussian identity  0.063     TRUE
#> Stan Binary (Logit)        0.060       0.227 binomial    logit -0.037     TRUE
#> Stan Binary (Probit)       0.056       0.241 binomial   probit -0.039     TRUE
```

``` r
# Create comprehensive visualization
plot_data <- comparison_table
plot_data$Method_Clean <- gsub("Built-in ", "", plot_data$Method)
plot_data$Method_Clean <- gsub("Stan ", "", plot_data$Method_Clean)
plot_data$CI_Width <- plot_data$PC_CI_Upper - plot_data$PC_CI_Lower
plot_data$Software <- ifelse(grepl("Built-in", plot_data$Method), "Built-in", "Stan")

# Alternative using dplyr (if available)
if (require(dplyr, quietly = TRUE)) {
  plot_data <- comparison_table %>%
    mutate(
      Method_Clean = gsub("Built-in ", "", Method),
      Method_Clean = gsub("Stan ", "", Method_Clean),
      CI_Width = PC_CI_Upper - PC_CI_Lower,
      Software = ifelse(grepl("Built-in", Method), "Built-in", "Stan")
    )
}

# ATE comparison plot
p1 <- ggplot(plot_data, aes(x = Method_Clean, y = PC_Mean, color = Software)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = PC_CI_Lower, ymax = PC_CI_Upper), 
                width = 0.2, linewidth = 1, position = position_dodge(width = 0.3)) +
  geom_point(aes(y = True_ATE), shape = 4, size = 5, color = "red") +
  facet_wrap(~ Outcome_Type, scales = "free_y", ncol = 2) +
  labs(
    title = "ATE Estimates: Built-in vs External Software Integration",
    subtitle = "Points show posterior means, bars show 95% CIs, red X shows true ATE",
    x = "Method", 
    y = "ATE Estimate",
    color = "Software"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold")
  )

print(p1)
```

## Posterior Distribution Comparisons

``` r
# Create data frame for all posterior distributions
create_posterior_data <- function() {
  data.frame(
    ATE = c(
      result_internal_continuous$pc,
      result_stan_continuous$pc,
      result_internal_binary$pc,
      result_stan_binary_logit$pc,
      result_stan_binary_probit$pc
    ),
    Method = rep(c("Built-in", "Stan", "Built-in", "Stan (Logit)", "Stan (Probit)"), 
                 times = c(
                   length(result_internal_continuous$pc),
                   length(result_stan_continuous$pc),
                   length(result_internal_binary$pc),
                   length(result_stan_binary_logit$pc),
                   length(result_stan_binary_probit$pc)
                 )),
    Outcome_Type = rep(c("Continuous", "Continuous", "Binary", "Binary", "Binary"),
                      times = c(
                        length(result_internal_continuous$pc),
                        length(result_stan_continuous$pc),
                        length(result_internal_binary$pc),
                        length(result_stan_binary_logit$pc),
                        length(result_stan_binary_probit$pc)
                      ))
  )
}

posterior_data <- create_posterior_data()

# Create density plots
p2 <- ggplot(posterior_data, aes(x = ATE, fill = Method)) +
  geom_density(alpha = 0.7) +
  geom_vline(data = data.frame(
    Outcome_Type = c("Continuous", "Binary"),
    true_ate = c(1.5, true_ate_binary)
  ), aes(xintercept = true_ate), linetype = "dashed", 
             color = "red", linewidth = 1) +
  facet_wrap(~ Outcome_Type, scales = "free", ncol = 1) +
  labs(
    title = "Posterior Distributions: Built-in vs External Software",
    subtitle = "Red line shows true ATE",
    x = "Average Treatment Effect",
    y = "Density",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold")
  )

print(p2)
```

## Best Practices for External Integration

### Formula Interface Guidelines

#### 1. Design Matrix Consistency

``` r
# Always check design matrix dimensions before fitting external models
outcome_matrix <- model.matrix(your_outcome_formula, data)
ps_matrix <- model.matrix(your_ps_formula, data)

cat("Expected parameter counts:\n")
cat("Outcome model:", ncol(outcome_matrix), "\n")
cat("PS model:", ncol(ps_matrix), "\n")

# Your external software must produce samples with these exact dimensions
```

#### 2. Family and Link Specification

`family` and `link` are auto-detected from the outcome, but the data
cannot say which link the external software used: a binary outcome could
have come from either a logit or a probit fit, and the two put the
coefficients on different scales. `ps.link` matters for the same reason.
Name all three when the samples come from elsewhere.

``` r
result <- drbayes_pc(
  outcome.formula = Y ~ A + X1 + X2,
  ps.formula = A ~ X1 + X2 + X3,
  data = data,
  family = "binomial",       # auto-detected, but say it anyway
  link = "logit",            # the scale the outcome draws are on
  ps.link = "logit",         # the scale the propensity score draws are on
  outcome.samples = external_outcome_samples,
  ps.samples = external_ps_samples
)
```

#### 3. Sample Dimension Validation

``` r
# Create validation function
validate_external_samples <- function(outcome_samples, ps_samples, 
                                     outcome.formula, ps.formula, data) {
  outcome_matrix <- model.matrix(outcome.formula, data)
  ps_matrix <- model.matrix(ps.formula, data)
  
  errors <- c()
  
  if (ncol(outcome_samples) != ncol(outcome_matrix)) {
    errors <- c(errors, sprintf("Outcome samples: expected %d cols, got %d", 
                               ncol(outcome_matrix), ncol(outcome_samples)))
  }
  
  if (ncol(ps_samples) != ncol(ps_matrix)) {
    errors <- c(errors, sprintf("PS samples: expected %d cols, got %d", 
                               ncol(ps_matrix), ncol(ps_samples)))
  }
  
  if (nrow(outcome_samples) != nrow(ps_samples)) {
    errors <- c(errors, "Sample counts don't match between outcome and PS models")
  }
  
  if (length(errors) > 0) {
    stop(paste(errors, collapse = "; "))
  } else {
    cat("Sample dimensions validated successfully!\n")
  }
}

# Use validation
validate_external_samples(outcome_samples_continuous, ps_samples_stan,
                         Y_continuous ~ A + X1 + X2 + X3 + A:X1,
                         A ~ X1 + X2 + X3 + X4, data)
```

### Software-Specific Integration Tips

#### Stan/cmdstanr

- Use `refresh = 0` for silent fitting in production
- Consider `parallel_chains = 4` for faster computation
- For binary outcomes, ensure correct likelihood specification

``` r
# Best practices for Stan integration
# 1. Use consistent parameterization
# 2. Include intercepts in design matrices
# 3. Match DRBayes prior specifications when possible
# 4. Use 'format = "matrix"' for draws extraction

# Example Stan workflow
fit <- cmdstan_model(stan_file)$sample(
  data = list(N = nrow(data), K = ncol(design_matrix), 
              X = design_matrix, y = response),
  chains = 4,
  iter_sampling = 1000,
  refresh = 0
)

samples <- fit$draws("beta", format = "matrix")
```

### Troubleshooting Guide

``` r
# 1. The design matrix and the samples disagree about how many coefficients
#    there are. This is the mistake the formula interface can catch, and it
#    names the two counts.
tryCatch({
  wrong_samples <- outcome_samples_continuous[, 1:3]  # the formula needs 6
  drbayes_pc(
    outcome.formula = Y_continuous ~ A + X1 + X2 + X3 + A:X1,
    ps.formula = A ~ X1 + X2 + X3 + X4,
    data = data,
    family = "gaussian",
    outcome.samples = wrong_samples,
    ps.samples = ps_samples_stan
  )
}, error = function(e) {
  cat("Expected error - dimension mismatch:\n")
  cat(conditionMessage(e), "\n\n")
})
#> Expected error - dimension mismatch:
#> outcome.samples must have 6 parameters (intercept + 5 covariates), but has 3

# 2. The two sets of samples must have the same number of iterations, since a
#    tilted draw pairs one outcome draw with one propensity score draw.
tryCatch({
  drbayes_pc(
    outcome.formula = Y_continuous ~ A + X1 + X2 + X3 + A:X1,
    ps.formula = A ~ X1 + X2 + X3 + X4,
    data = data,
    family = "gaussian",
    outcome.samples = outcome_samples_continuous,
    ps.samples = ps_samples_stan[1:500, ],
    bn = 0
  )
}, error = function(e) {
  cat("Expected error - unequal iteration counts:\n")
  cat(conditionMessage(e), "\n\n")
})
#> Expected error - unequal iteration counts:
#> outcome.samples and ps.samples must have the same number of iterations

# 3. What cannot be caught. Neither the formula nor the data records which
#    link the external software used, so probit draws handed over as logit
#    draws produce a number rather than an error. Name family, link and
#    ps.link yourself.
cat("Successful integration requires:\n")
#> Successful integration requires:
cat("1. Correct sample dimensions\n")
#> 1. Correct sample dimensions
cat("2. The same number of iterations in both sets of samples\n")
#> 2. The same number of iterations in both sets of samples
cat("3. family, link and ps.link matching the models that were fitted\n")
#> 3. family, link and ps.link matching the models that were fitted
```

## Conclusion

The **DRBayes** package provides seamless integration with popular
Bayesian software for both continuous and binary outcomes, allowing
researchers to:

1.  **Leverage existing expertise** with their preferred Bayesian
    software
2.  **Handle multiple outcome types** with appropriate model families
    and link functions
3.  **Maintain modeling flexibility** while benefiting from doubly
    robust inference
4.  **Compare results** across different software implementations and
    link functions
5.  **Scale to complex models** that may be difficult to implement in
    built-in functions

The posterior coupling approach ensures that regardless of which
software is used for initial model fitting, the final causal inference
benefits from the doubly robust properties and uncertainty
quantification that the method provides.

For more information about the theoretical foundations of the posterior
coupling method, see:

> Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly
> Robust Causal Inference via Posterior Coupling. arXiv preprint
> arXiv:2506.04868.

## Session Information

``` r
sessionInfo()
#> R version 4.5.1 (2025-06-13)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Tahoe 26.5.2
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8
#> 
#> time zone: Asia/Tokyo
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices datasets  utils     methods   base     
#> 
#> other attached packages:
#> [1] DRBayes_0.1.1
#> 
#> loaded via a namespace (and not attached):
#>  [1] tensorA_0.36.2.1     sass_0.4.10          future_1.67.0       
#>  [4] generics_0.1.4       renv_1.1.4           listenv_0.9.1       
#>  [7] digest_0.6.37        magrittr_2.0.3       evaluate_1.0.5      
#> [10] instantiate_0.2.3    mvtnorm_1.3-3        fastmap_1.2.0       
#> [13] jsonlite_2.0.0       processx_3.9.0       backports_1.5.1     
#> [16] ps_1.9.3             extraDistr_1.10.0    purrr_1.1.0         
#> [19] codetools_0.2-20     textshaping_1.0.5    jquerylib_0.1.4     
#> [22] abind_1.4-8          cli_3.6.5            rlang_1.3.0         
#> [25] cmdstanr_0.9.0       parallelly_1.45.1    withr_3.0.3         
#> [28] cachem_1.1.0         yaml_2.3.12          otel_0.2.0          
#> [31] tools_4.5.1          parallel_4.5.1       checkmate_2.3.4     
#> [34] pgdraw_1.1           globals_0.18.0       vctrs_0.6.5         
#> [37] posterior_1.7.1      R6_2.6.1             lifecycle_1.0.4     
#> [40] fs_2.1.0             ragg_1.5.2           furrr_0.3.1         
#> [43] pkgconfig_2.0.3      desc_1.4.3           callr_3.8.0         
#> [46] pkgdown_2.2.1        pillar_1.11.1        bslib_0.12.0        
#> [49] data.table_1.18.4    glue_1.8.0           Rcpp_1.1.2          
#> [52] RcppTN_0.2-2         systemfonts_1.3.2    xfun_0.60           
#> [55] tibble_3.3.1         knitr_1.51           htmltools_0.5.9     
#> [58] rmarkdown_2.31       compiler_4.5.1       distributional_0.8.1
```

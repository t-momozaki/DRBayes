// Logistic regression, the Stan counterpart of bayes_logit().
//
// The prior is the one that sampler documents, so the two backends target the
// same posterior:
//   theta ~ N(0, prior_precision^-1)
//   Y     ~ Bernoulli(logit^-1(Z theta))
// Z carries its intercept in column one, matching cbind(1, X) on the R side.

data {
  int<lower=0> N;
  int<lower=1> P;
  matrix[N, P] Z;
  array[N] int<lower=0, upper=1> Y;
  matrix[P, P] prior_precision;
}

transformed data {
  // The prior is data, so its factor is computed once rather than at every
  // gradient evaluation, which is what multi_normal_prec() would cost.
  matrix[P, P] L_prior = cholesky_decompose(prior_precision);
}

parameters {
  vector[P] theta;
}

model {
  // theta' Q theta = || L' theta ||^2 with Q = L L'.
  target += -0.5 * dot_self(L_prior' * theta);
  Y ~ bernoulli_logit_glm(Z, 0.0, theta);
}

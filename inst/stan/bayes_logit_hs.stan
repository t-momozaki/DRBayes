// Logistic regression with a horseshoe prior, the Stan counterpart of
// bayes_logit_hs().
//
// The priors are the ones that sampler documents, so the two backends target
// the same posterior:
//   theta[1]            ~ N(0, 1 / beta0_prior)
//   theta[unshrunk_idx] ~ N(0, 1 ./ unshrunk_prior)
//   theta[shrunk_idx]   ~ N(0, tau^2 lambda^2)
//   lambda              ~ Cauchy+(0, 1)
//   tau                 ~ Cauchy+(0, tau_prior)
//   Y                   ~ Bernoulli(logit^-1(Z theta))
//
// The coefficients listed in unshrunk_idx keep a normal prior instead of the
// horseshoe. In a causal model those are the treatment main effect and every
// effect modification term, because shrinking them attenuates the estimand.
//
// Everything is written non-centred: the coefficients are standard normal
// draws multiplied by their scales. This is the reason to offer Stan for the
// horseshoe at all, since the centred version has a funnel that no gradient
// based sampler steps through cleanly.

data {
  int<lower=0> N;
  int<lower=2> P;
  matrix[N, P] Z;
  array[N] int<lower=0, upper=1> Y;
  int<lower=0> K;
  int<lower=1> M;
  array[K] int<lower=2, upper=P> unshrunk_idx;
  array[M] int<lower=2, upper=P> shrunk_idx;
  real<lower=0> beta0_prior;
  vector<lower=0>[K] unshrunk_prior;
  real<lower=0> tau_prior;
}

transformed data {
  real beta0_sd = inv_sqrt(beta0_prior);
  vector[K] unshrunk_sd = inv_sqrt(unshrunk_prior);
}

parameters {
  real beta0_raw;
  vector[K] beta_unshrunk_raw;
  vector[M] beta_shrunk_raw;
  vector<lower=0>[M] lambda;
  real<lower=0> tau;
}

transformed parameters {
  vector[P] theta;
  theta[1] = beta0_raw * beta0_sd;
  theta[unshrunk_idx] = beta_unshrunk_raw .* unshrunk_sd;
  theta[shrunk_idx] = beta_shrunk_raw .* lambda * tau;
}

model {
  beta0_raw ~ std_normal();
  beta_unshrunk_raw ~ std_normal();
  beta_shrunk_raw ~ std_normal();
  lambda ~ cauchy(0, 1);
  tau ~ cauchy(0, tau_prior);
  Y ~ bernoulli_logit_glm(Z, 0.0, theta);
}

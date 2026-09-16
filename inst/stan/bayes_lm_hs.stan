// Gaussian linear model with a horseshoe prior, the Stan counterpart of
// bayes_lm_hs().
//
// The priors are the ones that sampler documents, so the two backends target
// the same posterior:
//   theta[1]            ~ N(0, 1 / beta0_prior)
//   theta[unshrunk_idx] ~ N(0, 1 ./ unshrunk_prior)
//   theta[shrunk_idx]   ~ N(0, sigma2 tau^2 lambda^2)
//   lambda              ~ Cauchy+(0, 1)
//   tau                 ~ Cauchy+(0, tau_prior)
//   sigma2              ~ InvGamma(sigma_a, sigma_b)
// The horseshoe scales carry a factor of sigma, as in Carvalho, Polson and
// Scott (2010), so their shrinkage does not move with the units of Y and tau
// is free of them. The two fixed priors do not, because a prior variance in
// units of sigma2 also adds its quadratic form to the residual sum of squares,
// and for an intercept and a treatment effect much larger than sigma that
// inflates the error variance. bayes_lm_hs() keeps those two weakly
// informative by scaling their precisions to the data instead.
//
// The coefficients listed in unshrunk_idx keep a normal prior instead of the
// horseshoe. In a causal model those are the treatment main effect and every
// effect modification term, because shrinking them attenuates the estimand.
//
// Everything is written non-centred: the coefficients are standard normal
// draws multiplied by their scales, sigma among them where it belongs, so no
// Jacobian is needed. This is the reason to offer Stan for the horseshoe at
// all, since the centred version has a funnel that no gradient based sampler
// steps through cleanly.

data {
  int<lower=0> N;
  int<lower=2> P;
  matrix[N, P] Z;
  vector[N] Y;
  int<lower=0> K;
  int<lower=1> M;
  array[K] int<lower=2, upper=P> unshrunk_idx;
  array[M] int<lower=2, upper=P> shrunk_idx;
  real<lower=0> beta0_prior;
  vector<lower=0>[K] unshrunk_prior;
  real<lower=0> tau_prior;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
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
  real<lower=0> sigma2;
}

transformed parameters {
  vector[P] theta;
  real sigma = sqrt(sigma2);
  theta[1] = beta0_raw * beta0_sd;
  theta[unshrunk_idx] = beta_unshrunk_raw .* unshrunk_sd;
  theta[shrunk_idx] = beta_shrunk_raw .* lambda * tau * sigma;
}

model {
  beta0_raw ~ std_normal();
  beta_unshrunk_raw ~ std_normal();
  beta_shrunk_raw ~ std_normal();
  lambda ~ cauchy(0, 1);
  tau ~ cauchy(0, tau_prior);
  sigma2 ~ inv_gamma(sigma_a, sigma_b);
  Y ~ normal_id_glm(Z, 0.0, theta, sigma);
}

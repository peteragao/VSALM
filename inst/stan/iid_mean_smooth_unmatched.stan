
functions {
  real pcprec_lpdf(real prec, real u, real alpha) {
    real lambda = -log(alpha)/u;
    real s = 1/sqrt(prec);
    real d = -log(2) - 1.5 * log(prec) + log(lambda) - lambda * s;
    return d;
  }
}

data {
  // DIMENSIONS
  int<lower=0> N;               // total number of areas
  int<lower=0> N_data;          // total number of areas with valid ests
  int<lower=0> K; // num covariates
  int ind_data[N_data];         // indices of non-missing data

  // DATA
  real Yhat[N_data]; // direct estimates of small area means
  real<lower=0> Vhat[N_data];  // variance estimates of direct estimates
  matrix[N, K] X; // design matrix

  // HYPERPARAMETERS
  real<lower=0> pc_u_v;
  real<lower=0> pc_u_alpha;
}

parameters {
// FIXED EFFECTS
  real mu;
  vector[K] betas; // covariates

  // RANDOM EFFECTS
  vector[N_data] u; // (scaled) area effects for sampled areas
  real<lower=0> sigma_u; // overall standard deviation

}
transformed parameters {
  vector[N] theta_obs;
  // variance of each component should be approximately equal to 1
  for (i in 1:N) {
    theta_obs[i] = mu;
  }
  if (K > 0) {
    theta_obs = theta_obs + X * betas;
  }
  for (i in 1:N_data) {
    theta_obs[ind_data[i]] = theta_obs[ind_data[i]] + sigma_u * u[i];
  }
  theta_obs = inv_logit(theta_obs);
}

model {
  int theta_k; // indexing thetas
  if (K > 0) {
    # precision ~ 1000
    betas ~ normal(0.0, 31.62278);
  }
  target += normal_lpdf(mu | 0, 31.62278);
  for (i in 1:N_data) {
    theta_k = ind_data[i];
    target += normal_lpdf(Yhat[i] | theta_obs[theta_k], sqrt(Vhat[i]));
  }
  target += pcprec_lpdf(1 / pow(sigma_u, 2) | pc_u_v,  pc_u_alpha);
  target += normal_lpdf(u | 0, 1);
}
generated quantities {
  vector[N] theta;
  for (i in 1:N) {
    theta[i] = inv_logit(normal_rng(logit(theta_obs[i]), sigma_u));
  }
  for (i in ind_data) {
    theta[i] = theta_obs[i];
  }
}


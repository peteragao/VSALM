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
  int<lower=0> na[N_data];      // number of sampled units in area
  int<lower=0> df[N_data];      // number of degrees of freedom

  // DATA
  real Yhat[N_data]; // direct estimates of small area means
  real<lower=0> Vhat[N_data];  // variance estimates of direct estimates
  matrix[N, K] X; // design matrix

  // HYPERPARAMETERS
  real<lower=0> pc_u_v;
  real<lower=0> pc_u_alpha;
  real<lower=0> pc_tau_v;
  real<lower=0> pc_tau_alpha;
}
parameters {
// FIXED EFFECTS
  real mu;
  vector[K] betas; // covariates

  // RANDOM EFFECTS
  vector[N_data] u; // (scaled) area effects for sampled areas
  real<lower=0> sigma_u; // overall standard deviation

  // GENERALIZED VARIANCE FUNCTION
  real g0;
  real g1;
  real g2;
  vector[N_data] tau; // gvf errors;
  real<lower=0> sigma_tau; // error in gvf

}
transformed parameters {
  vector[N] theta_obs;
  real<lower=0> V[N_data]; // sampling errors
  real<lower=0> sigma_e[N_data]; // sampling errors
  real<lower=0> ss_cl[N_data]; // cluster sum of squares

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
  for (i in 1:N_data) {
    sigma_e[i] = sqrt(exp(g0 + g1 * log(theta_obs[ind_data[i]] * (1 - theta_obs[ind_data[i]])) + g2 * log(na[i]) + sigma_tau * tau[i]));
    V[i] = pow(sigma_e[i], 2);
    ss_cl[i] = (df[i]) * Vhat[i] / V[i];
  }
}

model {
  int theta_k; // indexing thetas
  if (K > 0) {
    betas ~ normal(0.0, 1.0);
  }

  target += normal_lpdf(g0 | 0, 1);
  target += normal_lpdf(g1 | -1, .5);
  target += normal_lpdf(g2 | -1, .5);

  target += pcprec_lpdf(1 / pow(sigma_u, 2) | pc_u_v, pc_u_alpha);
  target += pcprec_lpdf(1 / pow(sigma_tau, 2) | pc_tau_v, pc_tau_alpha);
  target += normal_lpdf(mu | 0, 5);
  target += normal_lpdf(u | 0, 1);
  target += normal_lpdf(tau | 0, 1);
  for (i in 1:N_data) {
    theta_k = ind_data[i];
    target += chi_square_lpdf(ss_cl[i] | df[i]);
    target += normal_lpdf(Yhat[i] | theta_obs[theta_k], sigma_e[i]);
  }
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

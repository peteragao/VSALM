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
  int<lower=0> N_edges; // number of edges
  int<lower=0> K; // num covariates
  int ind_data[N_data];         // indices of non-missing data
  int<lower=0> df[N_data];      // number of degrees of freedom
  int<lower=1, upper=N> n1[N_edges]; // n1[i] adjacent to node2[i]
  int<lower=1, upper=N> n2[N_edges]; //

  // DATA
  real Yhat[N_data]; // direct estimates of small area means
  real<lower=0> Vhat[N_data];  // variance estimates of direct estimates
  matrix[N, K] X; // design matrix

  // HYPERPARAMETERS
  real<lower=0> scaling_factor; // scales the variance of the spatial effects
  real<lower=0> pc_u_v;
  real<lower=0> pc_u_alpha;
  real<lower=0> pc_tau_v;
  real<lower=0> pc_tau_alpha;
}
transformed data {
  real logit_y[N_data] = logit(Yhat);
}
parameters {
// FIXED EFFECTS
  real mu;
  vector[K] betas; // covariates

  // RANDOM EFFECTS
  vector[N] u_ns; // (scaled) area effects for sampled areas

  real<lower=0> sigma_u; // overall standard deviation
  real<lower=0, upper=1> phi; // proportion unstructured vs. spatially structured variance
  vector[N] u_sp; // spatial effects

  // GENERALIZED VARIANCE FUNCTION
  real g0;
  real g1;
  real g2;
  vector[N_data] tau; // gvf errors;
  real<lower=0> sigma_tau; // error in gvf

}
transformed parameters {
  vector[N] theta;
  real<lower=0> V[N_data]; // sampling errors
  real<lower=0> sigma_e[N_data]; // sampling errors
  real<lower=0> ss_cl[N_data]; // cluster sum of squares
  vector[N] u;
  // variance of each component should be approximately equal to 1
  u = sqrt(1 - phi) * u_ns + sqrt(phi / scaling_factor) * u_sp;
  theta = mu + sigma_u * u;
  if (K > 0) {
    theta = theta + X * betas;
  }
  theta = inv_logit(theta);
  for (i in 1:N_data) {
    sigma_e[i] = sqrt(exp(g0 + g1 * log(theta[ind_data[i]] * (1 - theta[ind_data[i]])) + g2 * log(df[i]) + sigma_tau * tau[i]));
    V[i] = pow(sigma_e[i], 2);
    ss_cl[i] = (df[i] - 1) * Vhat[i] / V[i];
  }
}

model {
  int theta_k; // indexing thetas
  if (K > 0) {
    betas ~ normal(0.0, 1.0);
  }
  phi ~ beta(0.5, 0.5);

  target += normal_lpdf(g0 | 0, 1);
  target += normal_lpdf(g1 | -1, .5);
  target += normal_lpdf(g2 | -1, .5);
  // This is the prior for u_sp! (up to proportionality)
  target += -0.5 * dot_self(u_sp[n1] - u_sp[n2]);
  // soft sum-to-zero constraint on u_sp)
  sum(u_sp) ~ normal(0, 0.001 * N); // equivalent to mean(u_sp) ~ normal(0,0.001)

  target += pcprec_lpdf(1 / pow(sigma_u, 2) | pc_u_v, pc_u_alpha);
  target += pcprec_lpdf(1 / pow(sigma_tau, 2) | pc_tau_v, pc_tau_alpha);
  target += normal_lpdf(mu | 0, 5);
  target += normal_lpdf(u_ns | 0, 1);
  target += normal_lpdf(tau | 0, 1);
  for (i in 1:N_data) {
    theta_k = ind_data[i];
    target += chi_square_lpdf(ss_cl[i] | df[i] - 1);
    target += normal_lpdf(Yhat[i] | theta[theta_k], sigma_e[i]);
  }
}

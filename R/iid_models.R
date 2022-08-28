#' Mean smoothing only model with iid area effects
#'
#' @param Yhat vector of direct weighted estimates
#' @param Vhat vector of estimated variances of direct weighted estimators
#' @param domain vector of domain labels
#' @param X matrix of area level covariates
#' @param pc_u hyperparameters for prior c(t, alpha) such that P(sigma_u > t) = alpha
#' @param var_tol tolerance parameter; all direct estimates with Vhat below var_tol are ignored when smoothing
#' @param initf optional stan parameter with starting parameter values
#' @param seed random seed
#' @param detailed_output if T, return stan object as well as estimates
#' @param chains number of chains for Stan
#' @param warmup number of warm up samples per chain for Stan
#' @param iter total number of samples per chain for Stan
#' @param ... additional parameters
#'
#' @return A data frame of summaries of the draws from the posterior, including mean, median, variance, and prediction interval
#' @export
#'
#' @examples
iidMeanSmooth <- function(Yhat, Vhat,
                          domain, X = NULL,
                          pc_u = c(1, .01),
                          var_tol = 1e-5,
                          initf = NULL, seed = NULL,
                          detailed_output = F,
                          chains = 4,
                          warmup = 2000,
                          iter = 4000,
                          ...) {
  if(is.null(seed)) {
    seed = 20220504
  }
  direct_est <- data.frame(
    domain = domain,
    Yhat = Yhat,
    Vhat = Vhat
  )
  to_smooth <- direct_est$Vhat > var_tol & !is.na(direct_est$Vhat)
  if (is.null(X)) {
    X <-  matrix(0, nrow = length(Yhat), ncol = 0)
  }

  dat <- list(N = length(to_smooth),
              N_data = sum(to_smooth),
              ind_data = which(to_smooth),
              K = ncol(X),
              Yhat = Yhat[to_smooth],
              Vhat = Vhat[to_smooth],
              X = X,
              pc_u_v = pc_u[1],
              pc_u_alpha = pc_u[2])
  if (is.null(initf)) {
    initf <- function() {
      list(theta_obs = rep(.5, length(to_smooth)),
           prec_u = runif(1, .001, .3),
           u = rnorm(sum(to_smooth), 0, .1))
    }
  }

  out_stan <- sampling(
    stanmodels$iid_mean_smooth,
    data = dat,         # named list of data
    chains = chains,         # number of Markov chains
    warmup = warmup,      # number of warmup iterations per chain
    iter = iter,        # total number of iterations per chain
    cores = 1,          # number of cores
    init = initf,
    seed = seed,
    ...
  )
  out_est <- summarize_estimates(domain, out_stan, "iidMeanSmooth", logit = F)
  if (detailed_output) {
    return(
      list(
        stan = out_stan,
        est = out_est
      )
    )
  }
  out_est
}
#' Logit-mean smoothing only model with iid area effects
#'
#' @param Yhat vector of direct weighted estimates
#' @param Vhat vector of estimated variances of direct weighted estimators
#' @param domain vector of domain labels
#' @param X matrix of area level covariates
#' @param pc_u hyperparameters for prior c(t, alpha) such that P(sigma_u > t) = alpha
#' @param var_tol tolerance parameter; all direct estimates with Vhat below var_tol are ignored when smoothing
#' @param initf optional stan parameter with starting parameter values
#' @param seed random seed
#' @param detailed_output if T, return stan object as well as estimates
#' @param chains number of chains for Stan
#' @param warmup number of warm up samples per chain for Stan
#' @param iter total number of samples per chain for Stan
#' @param ... additional parameters
#'
#' @return A data frame of summaries of the draws from the posterior, including mean, median, variance, and prediction interval
#' @export
#'
#' @examples
iidMeanSmoothLogit <- function(Yhat, Vhat,
                               domain, X = NULL,
                               pc_u = c(1, .01),
                               var_tol = 1e-5,
                               initf = NULL, seed = NULL,
                               detailed_output = F,
                               chains = 4,
                               warmup = 2000,
                               iter = 4000,
                               ...) {
  if(is.null(seed)) {
    seed = 20220504
  }
  direct_est <- data.frame(
    domain = domain,
    Yhat = SUMMER::logit(Yhat),
    Vhat = Vhat / Yhat^2 / (1-Yhat)^2
  )
  to_smooth <- direct_est$Vhat > var_tol & !is.na(direct_est$Vhat)
  if (is.null(X)) {
    X <-  matrix(0, nrow = length(Yhat), ncol = 0)
  }

  dat <- list(N = length(to_smooth),
              N_data = sum(to_smooth),
              ind_data = which(to_smooth),
              K = ncol(X),
              Yhat = direct_est$Yhat[to_smooth],
              Vhat = direct_est$Vhat[to_smooth],
              X = X,
              pc_u_v = pc_u[1],
              pc_u_alpha = pc_u[2])
  if (is.null(initf)) {
    initf <- function() {
      list(theta_obs = rep(.5, length(to_smooth)),
           prec_u = runif(1, .001, .3),
           u = rnorm(sum(to_smooth), 0, .1))
    }
  }

  out_stan <- sampling(
    stanmodels$iid_mean_smooth,
    data = dat,         # named list of data
    chains = chains,         # number of Markov chains
    warmup = warmup,      # number of warmup iterations per chain
    iter = iter,        # total number of iterations per chain
    cores = 1,          # number of cores
    init = initf,
    seed = seed,
    ...
  )
  out_est <- summarize_estimates(domain, out_stan, "iidMeanSmoothLogit", logit = T)
  if (detailed_output) {
    return(
      list(
        stan = out_stan,
        est = out_est
      )
    )
  }
  out_est
}
#' Unmatched mean smoothing only model with iid area effects
#'
#' @param Yhat vector of direct weighted estimates
#' @param Vhat vector of estimated variances of direct weighted estimators
#' @param domain vector of domain labels
#' @param X matrix of area level covariates
#' @param pc_u hyperparameters for prior c(t, alpha) such that P(sigma_u > t) = alpha
#' @param var_tol tolerance parameter; all direct estimates with Vhat below var_tol are ignored when smoothing
#' @param initf optional stan parameter with starting parameter values
#' @param seed random seed
#' @param detailed_output if T, return stan object as well as estimates
#' @param chains number of chains for Stan
#' @param warmup number of warm up samples per chain for Stan
#' @param iter total number of samples per chain for Stan
#' @param ... additional parameters
#'
#' @return A data frame of summaries of the draws from the posterior, including mean, median, variance, and prediction interval
#' @export
#'
#' @examples
iidMeanSmoothUnmatched <- function(Yhat, Vhat,
                                   domain, X = NULL,
                                   pc_u = c(1, .01),
                                   var_tol = 1e-5,
                                   initf = NULL, seed = NULL,
                                   detailed_output = F,
                                   chains = 4,
                                   warmup = 2000,
                                   iter = 4000,
                                   ...) {
  if(is.null(seed)) {
    seed = 20220504
  }
  direct_est <- data.frame(
    domain = domain,
    Yhat = Yhat,
    Vhat = Vhat
  )
  to_smooth <- direct_est$Vhat > var_tol & !is.na(direct_est$Vhat)
  if (is.null(X)) {
    X <-  matrix(0, nrow = length(Yhat), ncol = 0)
  }

  dat <- list(N = length(to_smooth),
              N_data = sum(to_smooth),
              ind_data = which(to_smooth),
              K = ncol(X),
              Yhat = Yhat[to_smooth],
              Vhat = Vhat[to_smooth],
              X = X,
              pc_u_v = pc_u[1],
              pc_u_alpha = pc_u[2])
  if (is.null(initf)) {
    initf <- function() {
      list(theta = rep(.5, length(to_smooth)),
           prec_u = runif(1, .001, .3),
           u = rnorm(sum(to_smooth), 0, .1))
    }
  }

  out_stan <- sampling(
    stanmodels$iid_mean_smooth_unmatched,
    data = dat,         # named list of data
    chains = chains,         # number of Markov chains
    warmup = warmup,      # number of warmup iterations per chain
    iter = iter,        # total number of iterations per chain
    cores = 1,          # number of cores
    init = initf,
    seed = seed,
    ...
  )
  out_est <- summarize_estimates(domain, out_stan, "iidMeanSmoothUnmatched", logit = F)
  if (detailed_output) {
    return(
      list(
        stan = out_stan,
        est = out_est
      )
    )
  }
  out_est
}
#' Joint smoothing model with iid area effects
#'
#' @param Yhat vector of direct weighted estimates
#' @param Vhat vector of estimated variances of direct weighted estimators
#' @param domain vector of domain labels
#' @param na vector of sample sizes for areas
#' @param df vector of degrees of freedom for areas
#' @param X matrix of area level covariates
#' @param pc_u hyperparameters for prior c(t, alpha) such that P(sigma_u > t) = alpha
#' @param pc_tau hyperparameters for prior c(t, alpha) such that P(prec_tau > t) = alpha
#' @param var_tol tolerance parameter; all direct estimates with Vhat below var_tol are ignored when smoothing
#' @param initf optional stan parameter with starting parameter values
#' @param seed random seed
#' @param detailed_output if T, return stan object as well as estimates
#' @param chains number of chains for Stan
#' @param warmup number of warm up samples per chain for Stan
#' @param iter total number of samples per chain for Stan
#' @param ... additional parameters
#'
#' @return A data frame of summaries of the draws from the posterior, including mean, median, variance, and prediction interval
#' @import Matrix
#' @export
#'
#' @examples
iidJointSmooth <- function(Yhat, Vhat,
                           domain, na, df,
                           X = NULL,
                           pc_u = c(1, .01),
                           pc_tau = c(1, .01),
                           var_tol = 1e-5,
                           initf = NULL, seed = NULL,
                           detailed_output = F,
                           chains = 4,
                           warmup = 2000,
                           iter = 4000,
                           ...) {
  if(is.null(seed)) {
    seed = 20220504
  }
  direct_est <- data.frame(
    domain = domain,
    Yhat = Yhat,
    Vhat = Vhat
  )
  to_smooth <- direct_est$Vhat > var_tol & !is.na(direct_est$Vhat)
  if (is.null(X)) {
    X <-  matrix(0, nrow = length(Yhat), ncol = 0)
  }
  dat <- list(N = length(to_smooth),
              N_data = sum(to_smooth),
              ind_data = which(to_smooth),
              na = na[to_smooth],
              df = df[to_smooth],
              K = ncol(X),
              Yhat = Yhat[to_smooth],
              Vhat = Vhat[to_smooth],
              X = X,
              pc_u_v = pc_u[1],
              pc_u_alpha = pc_u[2],
              pc_tau_v = pc_tau[1],
              pc_tau_alpha = pc_tau[2])
  if (is.null(initf)) {
    initf <- function() {
      list(theta = rep(SUMMER::logit(.05), length(to_smooth)),
           prec_u = runif(1, .001, .3),
           u = rnorm(sum(to_smooth), 0, .1),
           g0 = rnorm(1, 0, .5),
           g1 = rnorm(1, 1, .5),
           g2 = rnorm(1, 1, .5),
           prec_tau = runif(1, .001, .3),
           tau = rnorm(sum(to_smooth), 0, .1))
    }
  }
  out_stan <- sampling(
    stanmodels$iid_joint_smooth,
    data = dat,         # named list of data
    chains = chains,         # number of Markov chains
    warmup = warmup,      # number of warmup iterations per chain
    iter = iter,        # total number of iterations per chain
    cores = 1,          # number of cores
    init = initf,
    seed = seed,
    ...
  )
  out_est <- summarize_estimates(domain, out_stan, "iidJointSmooth", logit = F)
  if (detailed_output) {
    return(
      list(
        stan = out_stan,
        est = out_est
      )
    )
  }
  out_est
}
#' Logit joint smoothing model with iid area effects
#'
#' @param Yhat vector of direct weighted estimates
#' @param Vhat vector of estimated variances of direct weighted estimators
#' @param domain vector of domain labels
#' @param na vector of sample sizes for areas
#' @param df vector of degrees of freedom for areas
#' @param X matrix of area level covariates
#' @param pc_u hyperparameters for prior c(t, alpha) such that P(sigma_u > t) = alpha
#' @param pc_tau hyperparameters for prior c(t, alpha) such that P(sigma_tau > t) = alpha
#' @param var_tol tolerance parameter; all direct estimates with Vhat below var_tol are ignored when smoothing
#' @param initf optional stan parameter with starting parameter values
#' @param seed random seed
#' @param detailed_output if T, return stan object as well as estimates
#' @param chains number of chains for Stan
#' @param warmup number of warm up samples per chain for Stan
#' @param iter total number of samples per chain for Stan
#' @param ... additional parameters
#'
#' @return A data frame of summaries of the draws from the posterior, including mean, median, variance, and prediction interval
#' @export
#'
#' @examples
iidJointSmoothLogit <- function(Yhat, Vhat,
                                domain, na, df,
                                X = NULL,
                                pc_u = c(1, .01),
                                pc_tau = c(1, .01),
                                var_tol = 1e-5,
                                initf = NULL, seed = NULL,
                                detailed_output = F,
                                chains = 4,
                                warmup = 2000,
                                iter = 4000,
                                ...) {
  if(is.null(seed)) {
    seed = 20220504
  }
  direct_est <- data.frame(
    domain = domain,
    Yhat = Yhat,
    Vhat = Vhat
  )
  to_smooth <- direct_est$Vhat > var_tol & !is.na(direct_est$Vhat)
  if (is.null(X)) {
    X <-  matrix(0, nrow = length(Yhat), ncol = 0)
  }

  dat <- list(N = length(to_smooth),
              N_data = sum(to_smooth),
              ind_data = which(to_smooth),
              na = na[to_smooth],
              df = df[to_smooth],
              K = ncol(X),
              Yhat = Yhat[to_smooth],
              Vhat = Vhat[to_smooth],
              X = X,
              pc_u_v = pc_u[1],
              pc_u_alpha = pc_u[2],
              pc_tau_v = pc_tau[1],
              pc_tau_alpha = pc_tau[2])
  if (is.null(initf)) {
    initf <- function() {
      list(theta = rep(SUMMER::logit(.05), length(to_smooth)),
           prec_u = runif(1, .001, .3),
           u = rnorm(sum(to_smooth), 0, .1),
           g0 = rnorm(1, 0, .5),
           g1 = rnorm(1, 1, .5),
           g2 = rnorm(1, 1, .5),
           prec_tau = runif(1, .001, .3),
           tau = rnorm(sum(to_smooth), 0, .1))
    }
  }
  out_stan <- sampling(
    stanmodels$iid_joint_smooth_logit,
    data = dat,         # named list of data
    chains = chains,         # number of Markov chains
    warmup = warmup,      # number of warmup iterations per chain
    iter = iter,        # total number of iterations per chain
    cores = 1,          # number of cores
    init = initf,
    seed = seed,
    ...
  )
  out_est <- summarize_estimates(domain, out_stan, "iidJointSmoothLogit", logit = T)
  if (detailed_output) {
    return(
      list(
        stan = out_stan,
        est = out_est
      )
    )
  }
  out_est
}

#' Unmatched joint smoothing model with iid area effects
#'
#' @param Yhat vector of direct weighted estimates
#' @param Vhat vector of estimated variances of direct weighted estimators
#' @param domain vector of domain labels
#' @param na vector of sample sizes for areas
#' @param df vector of degrees of freedom for areas
#' @param X matrix of area level covariates
#' @param pc_u hyperparameters for prior c(t, alpha) such that P(sigma_u > t) = alpha
#' @param pc_tau hyperparameters for prior c(t, alpha) such that P(sigma_tau > t) = alpha
#' @param var_tol tolerance parameter; all direct estimates with Vhat below var_tol are ignored when smoothing
#' @param initf optional stan parameter with starting parameter values
#' @param seed random seed
#' @param detailed_output if T, return stan object as well as estimates
#' @param chains number of chains for Stan
#' @param warmup number of warm up samples per chain for Stan
#' @param iter total number of samples per chain for Stan
#' @param ... additional parameters
#'
#' @return A data frame of summaries of the draws from the posterior, including mean, median, variance, and prediction interval
#' @import Matrix
#' @export
#'
#' @examples
iidJointSmoothUnmatched <- function(Yhat, Vhat,
                                    domain, na, df,
                                    X = NULL,
                                    pc_u = c(1, .01),
                                    pc_tau = c(1, .01),
                                    var_tol = 1e-5,
                                    initf = NULL, seed = NULL,
                                    detailed_output = F,
                                    chains = 4,
                                    warmup = 2000,
                                    iter = 4000,
                                    ...) {
  if(is.null(seed)) {
    seed = 20220504
  }
  direct_est <- data.frame(
    domain = domain,
    Yhat = Yhat,
    Vhat = Vhat
  )
  to_smooth <- direct_est$Vhat > var_tol & !is.na(direct_est$Vhat)
  if (is.null(X)) {
    X <-  matrix(0, nrow = length(Yhat), ncol = 0)
  }

  dat <- list(N = length(to_smooth),
              N_data = sum(to_smooth),
              ind_data = which(to_smooth),
              na = na[to_smooth],
              df = df[to_smooth],
              K = ncol(X),
              Yhat = Yhat[to_smooth],
              Vhat = Vhat[to_smooth],
              X = X,
              pc_u_v = pc_u[1],
              pc_u_alpha = pc_u[2],
              pc_tau_v = pc_tau[1],
              pc_tau_alpha = pc_tau[2])
  if (is.null(initf)) {
    initf <- function() {
      list(theta = rep(SUMMER::logit(.05), length(to_smooth)),
           prec_u = runif(1, .001, .3),
           u = rnorm(sum(to_smooth), 0, .1),
           g0 = rnorm(1, 0, .5),
           g1 = rnorm(1, 1, .5),
           g2 = rnorm(1, 1, .5),
           prec_tau = runif(1, .001, .3),
           tau = rnorm(sum(to_smooth), 0, .1))
    }
  }
  out_stan <- sampling(
    stanmodels$iid_joint_smooth_unmatched,
    data = dat,         # named list of data
    chains = chains,         # number of Markov chains
    warmup = warmup,      # number of warmup iterations per chain
    iter = iter,        # total number of iterations per chain
    cores = 1,          # number of cores
    init = initf,
    seed = seed,
    ...
  )
  out_est <- summarize_estimates(domain, out_stan, "iidJointSmoothUnmatched", logit = F)
  if (detailed_output) {
    return(
      list(
        stan = out_stan,
        est = out_est
      )
    )
  }
  out_est
}

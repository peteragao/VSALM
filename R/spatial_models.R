
#' Summarize draws from posterior
#'
#' @param domain vector of domain labels
#' @param out_stan stan object from fit
#' @param method name of model
#' @param logit if T, apply logit transform to draws
#'
#' @return A data frame of summaries of the draws from the posterior, including mean, median, variance, and prediction interval
#' @import SUMMER
#'
#' @examples
summarize_estimates <- function(domain, out_stan, method, logit = F) {
  draws <- as.matrix(out_stan)
  out_summary <- summary(out_stan, probs = c(0.05, 0.25, 0.50, 0.75, 0.95))$summary
  combined_draws <- t(draws[, grep("theta\\[", colnames(draws))])
  out_theta <-
    out_summary[grep("theta\\[", rownames(out_summary)), ]
  if (logit) {
    combined_draws <- SUMMER::expit(combined_draws)
    out_theta <- SUMMER::expit(out_theta)
  }

  out_est <- data.frame(
    domain = domain,
    mean = NA, median = NA, var = NA,
    lower = NA, upper = NA
  )
  out_est$mean <- out_theta[, 1]
  out_est$median <- out_theta[, 6]
  out_est$var <- apply(combined_draws, 1, var, na.rm = T)
  out_est$lower <- out_theta[, 4]
  out_est$upper <- out_theta[, 8]
  out_est$method <-  method
  return(out_est)
}
#' Create list of nodes from adjacency matrix for spatial modeling
#'
#' @param adj_mat adjacency matrix
#'
#' @return list of two vectors of nodes that represent edges in graph
#' @import Matrix
#'
#' @examples
nb2mat_to_nodes <- function(adj_mat) {
  N_A <- nrow(adj_mat)
  N_edges <- sum(adj_mat != 0) / 2
  n1 <- vector(mode="numeric", length = N_edges)
  n2 <- vector(mode="numeric", length = N_edges)
  k <- 1
  for (i in 1:N_A) {
    for (j in i:N_A) {
      if (adj_mat[i, j] != 0) {
        n1[k] <- i
        n2[k] <- j
        k <- k + 1
      }
    }
  }
  return(list(n1 = n1, n2 = n2))
}

#' Prepare objects necessary for BYM2 model
#'
#' @param adj_mat adjacency matrix
#'
#' @return list of two vectors of nodes that represent edges in graph and scaling factor necessary for precision matrix
#' @import Matrix
#'
#' @examples
prepare_bym2 <- function(adj_mat) {
  nodes <- nb2mat_to_nodes(adj_mat)
  inla_adj <- sparseMatrix(i = nodes$n1, j = nodes$n2,
                           x = 1, symmetric = T)
  # ICAR precision matrix
  Q <- Diagonal(nrow(adj_mat), Matrix::rowSums(inla_adj)) - inla_adj
  Q_jit = Q + Diagonal(nrow(adj_mat)) * max(diag(Q)) * sqrt(.Machine$double.eps)

  Q_inv = inla.qinv(Q_jit, constr=list(A = matrix(1, 1, nrow(adj_mat)), e=0))

  #Compute the geometric mean of the variances, which are on the diagonal of Q.inv
  scl = exp(mean(log(diag(Q_inv))))
  return(list(n1 = nodes$n1, n2 = nodes$n2, scaling_factor = scl))
}



#' Joint smoothing model with BYM2 area effects
#'
#' @param Yhat vector of direct weighted estimates
#' @param Vhat vector of estimated variances of direct weighted estimators
#' @param domain vector of domain labels
#' @param na vector of sample sizes for areas
#' @param df vector of degrees of freedom for areas
#' @param adj_mat adjacency matrix
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
spatialJointSmooth <- function(Yhat, Vhat,
                               domain, na, df,
                               adj_mat, X = NULL,
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

  prepped <- prepare_bym2(adj_mat)

  dat <- list(N = length(to_smooth),
              N_data = sum(to_smooth),
              ind_data = which(to_smooth),
              na = na[to_smooth],
              df = df[to_smooth],
              N_edges = length(prepped$n1),
              n1 = prepped$n1,
              n2 = prepped$n2,
              K = ncol(X),
              Yhat = Yhat[to_smooth],
              Vhat = Vhat[to_smooth],
              X = X,
              scaling_factor = prepped$scaling_factor,
              pc_u_v = pc_u[1],
              pc_u_alpha = pc_u[2],
              pc_tau_v = pc_tau[1],
              pc_tau_alpha = pc_tau[2])
  if (is.null(initf)) {
    initf <- function() {
      list(theta = rep(SUMMER::logit(.05), length(to_smooth)),
           prec_u = runif(1, .001, .3),
           phi = runif(1, .001, .999),
           u_ns = rnorm(length(to_smooth), 0, .1),
           u_sp = rnorm(length(to_smooth), 0, .1),
           g0 = rnorm(1, 0, .5),
           g1 = rnorm(1, 1, .5),
           g2 = rnorm(1, 1, .5),
           prec_tau = runif(1, .001, .3),
           tau = rnorm(sum(to_smooth), 0, .1))
    }
  }

  out_stan <- sampling(
    stanmodels$spatial_joint_smooth,
    data = dat,         # named list of data
    chains = chains,         # number of Markov chains
    warmup = warmup,      # number of warmup iterations per chain
    iter = iter,        # total number of iterations per chain
    cores = 1,          # number of cores
    init = initf,
    seed = seed,
    ...
  )
  out_est <- summarize_estimates(domain, out_stan, "spatialJointSmooth", logit = F)
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
#' Unmatched joint smoothing model with BYM2 area effects
#'
#' @param Yhat vector of direct weighted estimates
#' @param Vhat vector of estimated variances of direct weighted estimators
#' @param domain vector of domain labels
#' @param na vector of sample sizes for areas
#' @param df vector of degrees of freedom for areas
#' @param adj_mat adjacency matrix
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
spatialJointSmoothUnmatched <- function(Yhat, Vhat,
                                        domain, na, df,
                                        adj_mat, X = NULL,
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

  prepped <- prepare_bym2(adj_mat)

  dat <- list(N = length(to_smooth),
              N_data = sum(to_smooth),
              ind_data = which(to_smooth),
              na = na[to_smooth],
              df = df[to_smooth],
              N_edges = length(prepped$n1),
              n1 = prepped$n1,
              n2 = prepped$n2,
              K = ncol(X),
              Yhat = Yhat[to_smooth],
              Vhat = Vhat[to_smooth],
              X = X,
              scaling_factor = prepped$scaling_factor,
              pc_u_v = pc_u[1],
              pc_u_alpha = pc_u[2],
              pc_tau_v = pc_tau[1],
              pc_tau_alpha = pc_tau[2])
  if (is.null(initf)) {
    initf <- function() {
      list(theta = rep(SUMMER::logit(.05), length(to_smooth)),
           prec_u = runif(1, .001, .3),
           phi = runif(1, .001, .999),
           u_ns = rnorm(length(to_smooth), 0, .1),
           u_sp = rnorm(length(to_smooth), 0, .1),
           g0 = rnorm(1, 0, .5),
           g1 = rnorm(1, 1, .5),
           g2 = rnorm(1, 1, .5),
           prec_tau = runif(1, .001, .3),
           tau = rnorm(sum(to_smooth), 0, .1))
    }
  }

  out_stan <- sampling(
    stanmodels$spatial_joint_smooth_unmatched,
    data = dat,         # named list of data
    chains = chains,         # number of Markov chains
    warmup = warmup,      # number of warmup iterations per chain
    iter = iter,        # total number of iterations per chain
    cores = 1,          # number of cores
    init = initf,
    seed = seed,
    ...
  )
  out_est <- summarize_estimates(domain, out_stan, "spatialJointSmoothUnmatched", logit = F)
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
#' Logit joint smoothing model with BYM2 area effects
#'
#' @param Yhat vector of direct weighted estimates
#' @param Vhat vector of estimated variances of direct weighted estimators
#' @param domain vector of domain labels
#' @param na vector of sample sizes for areas
#' @param df vector of degrees of freedom for areas
#' @param adj_mat adjacency matrix
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
spatialJointSmoothLogit <- function(Yhat, Vhat,
                                    domain, na, df,
                                    adj_mat, X = NULL,
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

  prepped <- prepare_bym2(adj_mat)

  dat <- list(N = length(to_smooth),
              N_data = sum(to_smooth),
              ind_data = which(to_smooth),
              na = na[to_smooth],
              df = df[to_smooth],
              N_edges = length(prepped$n1),
              n1 = prepped$n1,
              n2 = prepped$n2,
              K = ncol(X),
              Yhat = Yhat[to_smooth],
              Vhat = Vhat[to_smooth],
              X = X,
              scaling_factor = prepped$scaling_factor,
              pc_u_v = pc_u[1],
              pc_u_alpha = pc_u[2],
              pc_tau_v = pc_tau[1],
              pc_tau_alpha = pc_tau[2])
  if (is.null(initf)) {
    initf <- function() {
      list(theta = rep(SUMMER::logit(.05), length(to_smooth)),
           prec_u = runif(1, .001, .3),
           phi = runif(1, .001, .999),
           u_ns = rnorm(length(to_smooth), 0, .1),
           u_sp = rnorm(length(to_smooth), 0, .1),
           g0 = rnorm(1, 0, .5),
           g1 = rnorm(1, 1, .5),
           g2 = rnorm(1, 1, .5),
           prec_tau = runif(1, .001, .3),
           tau = rnorm(sum(to_smooth), 0, .1))
    }
  }

  out_stan <- sampling(
    stanmodels$spatial_joint_smooth_logit,
    data = dat,         # named list of data
    chains = chains,         # number of Markov chains
    warmup = warmup,      # number of warmup iterations per chain
    iter = iter,        # total number of iterations per chain
    cores = 1,          # number of cores
    init = initf,
    seed = seed,
    ...
  )
  out_est <- summarize_estimates(domain, out_stan, "spatialJointSmoothLogit", logit = T)
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
#' Unmatched mean smoothing only model with BYM2 area effects
#'
#' @param Yhat vector of direct weighted estimates
#' @param Vhat vector of estimated variances of direct weighted estimators
#' @param domain vector of domain labels
#' @param adj_mat adjacency matrix
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
spatialMeanSmooth <- function(Yhat, Vhat,
                              domain,
                              adj_mat, X = NULL,
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

  prepped <- prepare_bym2(adj_mat)

  dat <- list(N = length(to_smooth),
              N_data = sum(to_smooth),
              ind_data = which(to_smooth),
              N_edges = length(prepped$n1),
              n1 = prepped$n1,
              n2 = prepped$n2,
              K = ncol(X),
              Yhat = Yhat[to_smooth],
              Vhat = Vhat[to_smooth],
              X = X,
              scaling_factor = prepped$scaling_factor,
              pc_u_v = pc_u[1],
              pc_u_alpha = pc_u[2])
  if (is.null(initf)) {
    initf <- function() {
      list(theta = rep(.5, length(to_smooth)),
           prec_u = runif(1, .001, .3),
           phi = runif(1, .001, .999),
           u_ns = rnorm(length(to_smooth), 0, .1),
           u_sp = rnorm(length(to_smooth), 0, .1))
    }
  }

  out_stan <- sampling(
    stanmodels$spatial_mean_smooth,
    data = dat,         # named list of data
    chains = chains,         # number of Markov chains
    warmup = warmup,      # number of warmup iterations per chain
    iter = iter,        # total number of iterations per chain
    cores = 1,          # number of cores
    init = initf,
    seed = seed,
    ...
  )
  out_est <- summarize_estimates(domain, out_stan, "spatialMeanSmooth", logit = F)
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
#' Logit-mean smoothing only model with BYM2 area effects
#'
#' @param Yhat vector of direct weighted estimates
#' @param Vhat vector of estimated variances of direct weighted estimators
#' @param domain vector of domain labels
#' @param adj_mat adjacency matrix
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
spatialMeanSmoothLogit <- function(Yhat, Vhat,
                                   domain,
                                   adj_mat, X = NULL,
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

  prepped <- prepare_bym2(adj_mat)

  dat <- list(N = length(to_smooth),
              N_data = sum(to_smooth),
              ind_data = which(to_smooth),
              N_edges = length(prepped$n1),
              n1 = prepped$n1,
              n2 = prepped$n2,
              K = ncol(X),
              Yhat = direct_est$Yhat[to_smooth],
              Vhat = direct_est$Vhat[to_smooth],
              X = X,
              scaling_factor = prepped$scaling_factor,
              pc_u_v = pc_u[1],
              pc_u_alpha = pc_u[2])
  if (is.null(initf)) {
    initf <- function() {
      list(theta = rep(.5, length(to_smooth)),
           prec_u = runif(1, .001, .3),
           phi = runif(1, .001, .999),
           u_ns = rnorm(length(to_smooth), 0, .1),
           u_sp = rnorm(length(to_smooth), 0, .1))
    }
  }

  out_stan <- sampling(
    stanmodels$spatial_mean_smooth,
    data = dat,         # named list of data
    chains = chains,         # number of Markov chains
    warmup = warmup,      # number of warmup iterations per chain
    iter = iter,        # total number of iterations per chain
    cores = 1,          # number of cores
    init = initf,
    seed = seed,
    ...
  )
  out_est <- summarize_estimates(domain, out_stan, "spatialMeanSmoothLogit", logit = T)
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
#' Mean smoothing only model with BYM2 area effects
#'
#' @param Yhat vector of direct weighted estimates
#' @param Vhat vector of estimated variances of direct weighted estimators
#' @param domain vector of domain labels
#' @param adj_mat adjacency matrix
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
spatialMeanSmoothUnmatched <- function(Yhat, Vhat,
                                       domain,
                                       adj_mat, X = NULL,
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

  prepped <- prepare_bym2(adj_mat)

  dat <- list(N = length(to_smooth),
              N_data = sum(to_smooth),
              ind_data = which(to_smooth),
              N_edges = length(prepped$n1),
              n1 = prepped$n1,
              n2 = prepped$n2,
              K = ncol(X),
              Yhat = Yhat[to_smooth],
              Vhat = Vhat[to_smooth],
              X = X,
              scaling_factor = prepped$scaling_factor,
              pc_u_v = pc_u[1],
              pc_u_alpha = pc_u[2])
  if (is.null(initf)) {
    initf <- function() {
      list(theta = rep(.5, length(to_smooth)),
           prec_u = runif(1, .001, .3),
           phi = runif(1, .001, .999),
           u_ns = rnorm(length(to_smooth), 0, .1),
           u_sp = rnorm(length(to_smooth), 0, .1))
    }
  }

  out_stan <- sampling(
    stanmodels$spatial_mean_smooth_unmatched,
    data = dat,         # named list of data
    chains = chains,         # number of Markov chains
    warmup = warmup,      # number of warmup iterations per chain
    iter = iter,        # total number of iterations per chain
    cores = 1,          # number of cores
    init = initf,
    seed = seed,
    ...
  )
  out_est <- summarize_estimates(domain, out_stan, "spatialMeanSmoothUnmatched", logit = F)
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

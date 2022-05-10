#' Title
#'
#' @param Yhat
#' @param Vhat
#' @param domain
#' @param X
#' @param pc_u
#' @param var_tol
#' @param initf
#' @param seed
#' @param detailed_output
#' @param chains
#' @param warmup
#' @param iter
#' @param ...
#'
#' @return
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
           sigma_u = runif(1, .001, .3),
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
#' Title
#'
#' @param Yhat
#' @param Vhat
#' @param domain
#' @param X
#' @param pc_u
#' @param var_tol
#' @param initf
#' @param seed
#' @param detailed_output
#' @param chains
#' @param warmup
#' @param iter
#' @param ...
#'
#' @return
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
           sigma_u = runif(1, .001, .3),
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
#' Title
#'
#' @param Yhat
#' @param Vhat
#' @param domain
#' @param X
#' @param pc_u
#' @param var_tol
#' @param initf
#' @param seed
#' @param detailed_output
#' @param chains
#' @param warmup
#' @param iter
#' @param ...
#'
#' @return
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
           sigma_u = runif(1, .001, .3),
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
#' Title
#'
#' @param Yhat
#' @param Vhat
#' @param domain
#' @param df
#' @param X
#' @param pc_u
#' @param pc_tau
#' @param var_tol
#' @param initf
#' @param seed
#' @param ...
#'
#' @return
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
           sigma_u = runif(1, .001, .3),
           u = rnorm(sum(to_smooth), 0, .1),
           g0 = rnorm(1, 0, .5),
           g1 = rnorm(1, 1, .5),
           g2 = rnorm(1, 1, .5),
           sigma_tau = runif(1, .001, .3),
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
#' Title
#'
#' @param Yhat
#' @param Vhat
#' @param domain
#' @param df
#' @param X
#' @param pc_u
#' @param pc_tau
#' @param var_tol
#' @param initf
#' @param seed
#' @param detailed_output
#' @param chains
#' @param warmup
#' @param iter
#' @param ...
#'
#' @return
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
           sigma_u = runif(1, .001, .3),
           u = rnorm(sum(to_smooth), 0, .1),
           g0 = rnorm(1, 0, .5),
           g1 = rnorm(1, 1, .5),
           g2 = rnorm(1, 1, .5),
           sigma_tau = runif(1, .001, .3),
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

#' Title
#'
#' @param Yhat
#' @param Vhat
#' @param domain
#' @param df
#' @param X
#' @param pc_u
#' @param pc_tau
#' @param var_tol
#' @param initf
#' @param seed
#' @param ...
#'
#' @return
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
           sigma_u = runif(1, .001, .3),
           u = rnorm(sum(to_smooth), 0, .1),
           g0 = rnorm(1, 0, .5),
           g1 = rnorm(1, 1, .5),
           g2 = rnorm(1, 1, .5),
           sigma_tau = runif(1, .001, .3),
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

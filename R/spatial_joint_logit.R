#' Title
#'
#' @param area_vec
#' @param data_list
#' @param initf hello
#' @param method_name
#' @param detailed_output
#'
#' @return
#' @export
#'
#' @examples
spatial_joint_logit <- function(area_vec, data_list,
                                initf, method_name, detailed_output = F, ...) {
  if (is.null(model)) {
    model <- stan_model(filename)
  }
  out_stan <- sampling(
    stanmodels$spatial_joint_logit,
    data = data_list,    # named list of data
    chains = 4,                # number of Markov chains
    warmup = 3000,             # number of warmup iterations per chain
    iter = 6000,               # total number of iterations per chain
    cores = 1,                 # number of cores
    init=initf,
    seed = 20220301,
    ...
  )
  out_summary <- summary(out_stan, probs = c(0.05, 0.25, 0.50, 0.75, 0.95))$summary
  out_theta <-
    out_summary[grep("theta_pred", rownames(out_summary)), ]
  out_est <- data.frame(
    domain = area_vec,
    mean = NA, median = NA, var = NA,
    lower = NA, upper = NA
  )
  out_est$mean <- out_theta[, 1]
  out_est$median <- out_theta[, 6]
  out_est$var <- out_theta[,3]^2
  out_est$lower <- out_theta[, 4]
  out_est$upper <- out_theta[, 8]
  out_est$method <-  method_name
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

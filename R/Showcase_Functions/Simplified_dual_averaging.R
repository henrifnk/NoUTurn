#' update stepsize
#'
#' updates stpesize parameters by performing dual averaging
#'
#' @inheritParams sample_noUturn
#' @param dap dual_avereging_params parameters necessary for dual averaging
#' @param iter actual iteration cycle
#' @param average_acceptance average achieved acceptance in iteration m
#'
update_stepsize <- function(dap, iter, target_accpentance, average_acceptance) {
  weight <- (iter + dap$stability)^(-1)
  dap$mcmc_behavior[iter + 1] = (1 - weight) * dap$mcmc_behavior[iter] + weight * (target_accpentance - average_acceptance)
  dap$stepsize[iter + 1] <- exp(dap$level - ((sqrt(iter) / dap$shrinkage) * dap$mcmc_behavior[iter + 1]))
  dap$stepsize_weight[iter + 1] <- exp(iter^(-dap$adaption) * log(dap$stepsize[iter + 1]) +
                                         (1 - iter^(-dap$adaption)) * log(dap$stepsize_weight[iter]))

  dap
}

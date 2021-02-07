find_initial_stepsize <- function(position) {
  stepsize <- 1
  momentum <- rnorm(length(position))
  proposal <- leapfrog(position, momentum, stepsize)
  exp <- 2 * (acceptance_rate(position, momentum, proposal) > 0.5) - 1
  while(acceptance_rate(position, momentum, proposal)^exp > 2^(-exp)) {
    stepsize <- stepsize * 2^exp
    proposal <- leapfrog(position, momentum, stepsize)
  }
  stepsize
}

acceptance_rate <- function(position, momentum, proposal) {
  proposal_dens <- do.call(joint_log_density, proposal)
  initial_dens <- joint_log_density(position, momentum)
  exp(proposal_dens - initial_dens)
}

#' Initialize parameters for dual averaging
#'
#' @param stepsize_weight weighted stepsize from previous iteration
#' @param mcmc_behavoir behavior of algorithm at previous iterations
#' @param shrinkage controls the amount of shrinkage
#' @param stability controls stability at early iterations
#' @param adaption controls the speed of converges and so the influence of more recent weights
#' @example dap <- init_parmeters(0.25)
#'
init_parmeters <- function(stepsize, stepsize_weight = 1,
                           mcmc_behavior = 0, shrinkage = 0.05,
                           stability = 10, adaption = 0.75) {
  level <- log(10 * stepsize) # where stepsize is shrunken towards (\mhy)
  list(
    "stepsize" = stepsize,
    "level" = level, "stepsize_weight" = stepsize_weight,
    "mcmc_behavior" = mcmc_behavior, "shrinkage" = shrinkage,
    "stability" = stability, "adaption" = adaption
  )
}


#' update stepsize
#'
#' updates stpesize parameters by performing dual averaging
#'
#' @inheritParams sample_noUturn
#' @param dap dual averaging parameters (seee init_parameters)
#' @param iter actual iteration cycle
#' @param average_acceptance average achieved acceptance in iteration m
#' @param target_acceptance aimed average acceptance level per iteration
#' @example dap <- update_stepsize(dap, 1, 0.65, 0.5)
#'
update_stepsize <- function(dap, iter, target_accpentance, average_acceptance) {
  list2env(dap, environment())
  dif <- target_accpentance - average_acceptance
  weight <- (iter + stability)^(-1)
  dap$mcmc_behavior[iter + 1] <- (1 - weight) * mcmc_behavior[iter] + weight * dif
  dap$stepsize[iter + 1] <- exp(level - ((sqrt(iter) / shrinkage) * dap$mcmc_behavior[iter + 1]))
  dap$stepsize_weight[iter + 1] <- exp(iter^(-adaption) * log(dap$stepsize[iter + 1]) +
                                         (1 - iter^(-adaption)) * log(stepsize_weight[iter]))
  dap
}

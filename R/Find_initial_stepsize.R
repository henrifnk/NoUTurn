#' Find A Reasonable initial value for stepsize parameter
#'
#' supplement to dual averaging to speed up convergence by initializing a reasonable stepsize via heuristic
#'
#' @inheritParams nouturn_sampler
find_initial_stepsize <- function(position_init) {
  stepsize <- 1
  momentum <- rnorm(length(position_init))
  proposal <- leapfrog(position_init, momentum, stepsize)
  exponent <- 2 * (acceptance_rate(position_init, momentum, proposal) > 0.5) -1
  while(acceptance_rate(position_init, momentum, proposal)^exponent > 2^(-exponent)) {
    stepsize <- stepsize * 2^exponent
    proposal <- leapfrog(position_init, momentum, stepsize)
  }
  stepsize
}
# ______________________________________________________________________________
acceptance_rate <- function(position_init, momentum, proposal) {
  do.call(joint_probability, proposal) / joint_probability(position_init, momentum)
}

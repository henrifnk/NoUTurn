#' Find A Reasonable initial value for stepsize parameter
#'
#' supplement to dual averaging to speed up convergence by initializing a reasonable stepsize via heuristic
#'
#' @inheritParams nouturn_sampler
find_initial_stepsize <- function(position_init, design = NULL, target = NULL) {
  stepsize <- 1
  momentum <- rnorm(length(position_init))
  args <- if(is.null(design) && is.null(target)) {
    list(position_init)
  } else list(position_init, design, target)
  gradient_step <- do.call(gradient, args)
  proposal <- leapfrog(position_init, momentum, stepsize, gradient_step)
  exponent <- 2 * (acceptance_rate(position_init, momentum, proposal) > 0.5) -1
  while(acceptance_rate(position_init, momentum, proposal)^exponent > 2^(-exponent)) {
    stepsize <- stepsize * 2^exponent
    proposal <- leapfrog(position_init, momentum, stepsize, gradient_step)
  }
  stepsize
}
# ______________________________________________________________________________
acceptance_rate <- function(position_init, momentum, proposal) {
  exp(do.call(joint_probability, proposal)) / exp(joint_probability(position_init, momentum))
}

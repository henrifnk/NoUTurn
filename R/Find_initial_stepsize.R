#' Find A Reasonable initial value for stepsize parameter
#'
#' supplement to dual averaging to speed up convergence by initializing a reasonable stepsize via heuristic
#'
#' @inheritParams nouturn_sampler
find_initial_stepsize <- function(position_init, design = NULL, target = NULL, is_log) {
  stepsize <- 1
  momentum <- rnorm(length(position_init))
  args <- if(is.null(design) && is.null(target)) {
    list("pos" = position_init)
  } else list("pos" = position_init, "des" = design, "tar" = target)
  gradient_step <- do.call(gradient, args)
  proposal <- leapfrog(position_init, momentum, stepsize, gradient_step)
  exponent <- 2 * (acceptance_rate(args, momentum, proposal, is_log) > 0.5) - 1
  while(acceptance_rate(args, momentum, proposal, is_log)^exponent > 2^(-exponent)) {
    stepsize <- stepsize * 2^exponent
    proposal <- leapfrog(position_init, momentum, stepsize, gradient_step)
  }
  stepsize
}
# ______________________________________________________________________________
acceptance_rate <- function(args, momentum, proposal, is_log) {
  proposal_dens <- joint_probability(proposal[[1]], proposal[[2]], args$des, args$tar, is_log)
  initial_dens <- joint_probability(args$pos, momentum, args$des, args$tar, is_log)
  exp(proposal_dens - initial_dens)
}

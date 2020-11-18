#' Naive U-Turn-Sampler
#'
#' @param theta_init initial value for theta
#' @param joint_prob log of joint density probability
#' @param iteration number of iterations to sample from
#' @param seed set a random seed to guarantee reproducability
#' @example posterior_density = mvtnorm::dmvnorm
#' @return
#'
#TODO: calculate posterior from data up to normalization constant
naive_nouturn_sampler <- function(position_init, stepsize, iteration, seed = 123L) {
  position <- position_init
  set.seed(seed)
  for(m in 1:(iteration + 1)) {
    momentum <- rnorm(length(position))
    slice <- runif(1L, max = joint_probability(position, momentum))
    valid_states <- structure(vector("list", length = 2L), names= c("position", "momentum"))
    # Initialize state to call on Build Tree
    state <- initialize_state(position, momentum, slice, stepsize)

    while(state$run) {
      # 1 means forward, -1 means backward doubling
      state$direction <- sample(c(-1, 1), 1L)
      state1 <- build_trees(state)

   }
  }
}

#'  Joint Probability of position and momentum
#'  @inheritParams leapfrog
#'
joint_probability <- function(position, momentum) {
  as.numeric(exp(log(posterior_density(position)) - 0.5 * t(momentum) %*% momentum))
}

#' Initialize state
#'
#' Helperfunction to intialize state
#' @inheritParams leapfrog
#' @inheritParams build_tree
initialize_state <- function(position, momentum, slice, stepsize) {
  list(
    "run" = 1L, "tree_depth" = 0L, "slice_sample" = slice, "direction" = 0L, "stepsize" = stepsize,
    "state_ritghmost"= list("position" = position, "momentum" = momentum),
    "state_leftmost" = list("position" = position, "momentum" = momentum),
    "valid_state" = valid_states
  )
}
